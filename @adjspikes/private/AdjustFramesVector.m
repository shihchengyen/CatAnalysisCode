function [spikes,frame_vector,removed_spikes] = AdjustFramesVector(sp,stimInfo,frame_vector);
%%% This function will take the original spiketrain and the original frames_vector
%%% and adjust both to account for extra syncs and create a new
%%% spiketrain and new frames_vector to be used for generating the raster plot.
%%% Both the spikes and frames_vector will be in milliseconds.

extra_syncs_frames = stimInfo.data.extrasyncs.frames;
if isfield(stimInfo.data.iniInfo,'total_number_frames_displayed')
    if ~isempty(extra_syncs_frames)
        if extra_syncs_frames(end)>stimInfo.data.iniInfo.total_number_frames_displayed
            fprintf(['Extra Sync on End Trigger','\n'])
            extra_syncs_frames(end)=[];
            num_extra_syncs = stimInfo.data.extrasyncs.numExtraS;
            num_extra_syncs(end)=[];
        else
            num_extra_syncs = stimInfo.data.extrasyncs.numExtraS;
        end
    else
        num_extra_syncs = stimInfo.data.extrasyncs.numExtraS;
    end
else
    num_extra_syncs = stimInfo.data.extrasyncs.numExtraS;
end
% Some old files have a -1 extra sync on the end frame, this is an error
% and will trip up the program, so we remove this error.
if ~isempty(num_extra_syncs)
    if num_extra_syncs(end)==-1
        extra_syncs_frames(end)=[];
        num_extra_syncs(end)=[];
    end
end
sync_duration = stimInfo.data.catInfo.sync_duration;
syncs_per_frame = stimInfo.data.iniInfo.refreshes_per_frame;
frame_duration = sync_duration*syncs_per_frame;

% Since the new UEI data is broadband and broken up into chunks, the
% ispikes file is trial based, equivalent to the chunks, so I concatenate
% the chunks into one long spiketrain as in the experiment
if isfield(stimInfo.data.iniInfo,'DAQ')
    if strcmp(stimInfo.data.iniInfo.DAQ,'UEI')
        numTrials=sp.data.numTrials;
        duration=sp.data.duration*1000; % Convert to Milliseconds
        time_add=0;
        spikes=[];
        for b = 1:numTrials
            spikes=[spikes (sp.data.trial(b).cluster.spikes+time_add)];
            time_add=time_add+duration;
        end
    else
        spikes = sp.data.trial.cluster.spikes;
    end
else
    spikes = sp.data.trial.cluster.spikes;
end

%%%%% Confirm the frames numbers where extra syncs occurred
if isfield(stimInfo.data.iniInfo,'obj_type')
    extra_syncs_ind = find(diff(frame_vector) > frame_duration+(sync_duration/2));
    missing_syncs_ind = find(diff(frame_vector) < frame_duration-(sync_duration/2));

else % Some ini info does not contain an obj_type field
    extra_syncs_ind = find(diff(frame_vector) > frame_duration+(sync_duration/2)); % Half the sync duration.
    missing_syncs_ind = find(diff(frame_vector) < frame_duration-(sync_duration/2));
end
numSpikes = length(spikes);

c = setdiff(extra_syncs_frames,extra_syncs_ind);
if ~isempty(c)
    error('EXTRA SYNC ERROR')
    return
end

% Flag for missing syncs
if ~isempty(missing_syncs_ind)
    error('MISSING SYNC DETECTED')
    return
end
removed_spikes=[];
%%%%% adjust the frames_vector along with the spiketrain.
if ~isempty(extra_syncs_ind)
    %% Find the amount of temporal shift due to extra syncs
    esd = num_extra_syncs*sync_duration;
    % Remove spikes that occured during extra syncs, otherwise sub-millisecond ISI's.
    for e = 1:length(extra_syncs_ind)
        ex_spikes = find(spikes > (frame_vector(extra_syncs_ind(e))+frame_duration) & spikes <= frame_vector(extra_syncs_ind(e)+1));
        if ~isempty(ex_spikes)
            removed_spikes = [removed_spikes spikes(ex_spikes)];
            spikes(ex_spikes) = [];
        end
    end
    if numSpikes~=(length(spikes)+length(removed_spikes))
        error('MISSING SPIKES')
        return
    end
    %%%%%%%%% Now adjust both the original spiketrain and the frame_vector
    %%%%%%%%% to account for the extrasyncs durations %%%%%%%%%
    for e = 1:length(extra_syncs_ind)
        adj_spikes = find(spikes > frame_vector(extra_syncs_ind(e))+frame_duration);
        if ~isempty(adj_spikes)
            if spikes(adj_spikes(1)) == spikes(1)
                spikes = spikes(adj_spikes)-esd(e);
            else
                spikes = [spikes(1:adj_spikes(1)-1) spikes(adj_spikes)-esd(e)];
            end
        end
        frame_vector = [frame_vector(1:extra_syncs_ind(e)) frame_vector(extra_syncs_ind(e)+1:end)-esd(e)];
    end
end
