function [raster,new_frame_vector] = CreateRaster(adjusted_frames_vector,adjusted_spiketrain,stimInfo)
% Generates the raster of spiketimes for repeated stimuli.

if strcmpi(stimInfo.data.iniInfo.type,'Movie')
    stimulus = 'Movie';
elseif strcmpi(stimInfo.data.iniInfo.type,'sparse_noise')
    if strcmpi(stimInfo.data.iniInfo.obj_type,'Grating')
        stimulus = 'Grating';
    elseif strcmpi(stimInfo.data.iniInfo.obj_type,'Square Grating')
        stimulus = 'Square Grating';
    elseif strcmpi(stimInfo.data.iniInfo.obj_type,'Bar')
        stimulus = 'Bar';
    elseif strcmpi(stimInfo.data.iniInfo.obj_type,'Square')
        stimulus = 'Square';
    end
elseif strcmpi(stimInfo.data.iniInfo.type,'m_sequence')
    stimulus = 'MSequence';
end

new_frame_vector=[];
if strcmp(stimulus,'Movie')
    fv=adjusted_frames_vector;
    st=adjusted_spiketrain;
    num_reps=stimInfo.data.catInfo.num_repetitions;
    num_frames=stimInfo.data.catInfo.num_frames;
    %%% Creation of the raster %%%%%%
    for ii = 1:num_reps
        %%%% Frames for the repetition %%%%
        if isfield(stimInfo.data.iniInfo,'spontaneous_activity_duration_in_frames')
            num_spon=stimInfo.data.iniInfo.spontaneous_activity_duration_in_frames;
            frames = fv((ii*(num_frames+num_spon)-(num_frames+num_spon))+1:ii*(num_frames+num_spon)+1);
        else
            frames = fv(ii*num_frames-num_frames+1:(ii*num_frames)+1);
        end
        %%%% Spiketimes for the given repetition %%%%%
        rep_spikes = st(find(st > frames(1) & st <= frames(end)));
        if ii > 1
            raster = concatenate(raster,rep_spikes-time_adjust);
            time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
        else
            raster = rep_spikes-frames(1);
            time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
        end
    end
elseif strcmp(stimulus,'Grating')
    num_stimulus_repetitions = stimInfo.data.catInfo.num_stimulus_repetitions;
    isi=1;
    %%% Creation of the raster %%%%%%
    for ii = 1:num_stimulus_repetitions
        %%%% Frames for the repetition %%%%
        if ii == 1
            if stimInfo.data.catInfo.num_stimulus_repetitions>stimInfo.data.catInfo.num_blank_repetitions
                num_stimulus_frames = stimInfo.data.catInfo.num_stimulus_frames-1+isi; % Due to A3 and A4 missing first frame
            else
                num_stimulus_frames = stimInfo.data.catInfo.num_stimulus_frames+isi;
            end
            frames = adjusted_frames_vector(1:num_stimulus_frames+1);
            count = num_stimulus_frames+1;
        elseif ii == num_stimulus_repetitions
            if stimInfo.data.catInfo.num_stimulus_repetitions>stimInfo.data.catInfo.num_blank_repetitions
                num_stimulus_frames = stimInfo.data.catInfo.num_stimulus_frames; % Due to A3 and A4 missing first frame
            else
                num_stimulus_frames = stimInfo.data.catInfo.num_stimulus_frames+isi;
            end
            frames = adjusted_frames_vector(count:count+num_stimulus_frames);
            count = count+num_stimulus_frames;
        else
            num_stimulus_frames = stimInfo.data.catInfo.num_stimulus_frames+isi;
            frames = adjusted_frames_vector(count:count+num_stimulus_frames);
            count = count+num_stimulus_frames;
        end
        %%%% Spiketimes for the given repetition %%%%%
        rep_spikes = adjusted_spiketrain(find(adjusted_spiketrain > frames(1) & adjusted_spiketrain <= frames(end)));
        if ii > 1
            raster = concatenate(raster,rep_spikes-time_adjust);
            time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
            new_frame_vector{ii} = frames(1:end-1); %take out the ISI Period
        else
            if stimInfo.data.catInfo.num_stimulus_repetitions>stimInfo.data.catInfo.num_blank_repetitions
                raster = rep_spikes+stimInfo.data.catInfo.frame_duration;
                time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
                new_frame_vector{ii} = frames(1:end-1); %take out the ISI Period
            else
                raster = rep_spikes-frames(1);
                time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
                new_frame_vector{ii} = frames(1:end-1); %take out the ISI Period
            end
        end
    end

    %%%% Obtain the rasters with the highest firing rate %%%%%%
    stim_order = stimInfo.data.catInfo.seq_file(:,3);
    stim_num = unique(stim_order);
    for r = 1:length(stim_num)
        index = find(stim_order==stim_num(r));
        spike_counts(r) = length(find(raster(index,:)>0));
    end
    [xs,ind] = max(spike_counts);
    stim_num = stim_num(ind);
    max_ind = find(stim_order==stim_num);
    raster = raster(max_ind,:);

    new_frame_vector = cell2array(new_frame_vector)';
    new_frame_vector = new_frame_vector(max_ind,:);
    %figure
    %ii = transpose([1:size(raster,1)])*ones(1,size(raster,2));
    %plot(raster,ii,'LineStyle','none','Marker','.','Color','b');

elseif strcmp(stimulus,'Square Grating')

    num_reps = stimInfo.data.catInfo.num_repetitions;
    num_grating_stimuli = stimInfo.data.catInfo.num_grating_stimuli;
    num_frames_per_grating = stimInfo.data.iniInfo.individual_stimulus_duration_in_frames;
    if isfield(stimInfo.data.iniInfo,'spontaneous_activity_duration_in_frames')
        num_spon=stimInfo.data.iniInfo.spontaneous_activity_duration_in_frames;
    end
    % This includes the spontaneous activity, which is displayed after each
    % repeated block of stimuli.
    total_num_stimuli=num_grating_stimuli*num_reps;
    spac_frames = cumsum(zeros(1,num_reps)+num_grating_stimuli);
    fv=adjusted_frames_vector;
    st=adjusted_spiketrain;

    %%% Creation of the raster %%%%%%
    for a = 1:total_num_stimuli
        if a==1
            ff=1; lf=num_frames_per_grating+1;
        elseif a==spac_frames(1)
            spac_frames(1)=[]; % Add Spontaneous Activity Frames
            ff=lf; lf=ff+num_frames_per_grating+num_spon;
        else
            ff=lf; lf=ff+num_frames_per_grating;
        end
        %%%%% Spiketimes for the given repetition %%%%%
        rep_spikes = st(find(st > fv(ff) & st <= fv(lf)));
        if a > 1
            raster = concatenate(raster,rep_spikes-time_adjust);
            time_adjust = fv(lf); % This adjusts the spiketimes for each repetition of the stimulus
        else
            raster = rep_spikes;
            time_adjust = fv(lf); % This adjusts the spiketimes for each repetition of the stimulus
        end
    end
   
elseif strcmp(stimulus,'MSequence')
    %%% For the M-sequence
    raster = adjusted_spiketrain;
elseif strcmp(stimulus,'Bar');
    %%% For Sparse Bar Stimulus
    raster = adjusted_spiketrain;
elseif strcmp(stimulus,'Square')
    %%% For Sparse Bar Stimulus
    raster = adjusted_spiketrain;
else
    fprintf('Stimulus Type is Not Found')
    raster=[];
end