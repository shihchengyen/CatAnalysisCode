function [S,Values] = Sparseness(obj,stimInfo,Args)
% The function will calculate the sparseness of the firing based on the
% psth at the frame resolution, using the vinje and gallant sparseness
% equation.

FramePoints = obj.data.adjFramePoints;
numFrames = stimInfo.data.catInfo.num_frames;

if(Args.IgnoreFrames)
    % find end time of last frame
    endtime = FramePoints(numFrames+1);
    % set up histogram bins
    histbins = 0:Args.WindowSize:endtime;
    frameCounts = histcie(obj.data.raster,histbins,'DropLast');
    % average over repetitions to get mean count for each window
    Values = mean(frameCounts,2);
    % get number of bins
    nbins = length(Values);
    % A is from Vinje and Gallant, 2000: with numerator being (sum(Values)/numFrames)^2
    % and denominator being sum( (Values^2) / numFrames)
    % but since there is 1/numFrames^2 in the numerator and 1/numFrames in
    % the denominator, A simplifies to ( using Values' * Values to simplify
    % calculation: of sum(Values^2) )
    A = ( sum(Values)^2 ) / ( (Values' * Values) * nbins );
    S = (1 - A) / (1 - (1/nbins)) * 100;
else
    Spiketrain = obj.data.adjSpiketrain;
    % Flag for multiunit data, specifically for the long movie responses.
    % Removes most of the double triggers by looking for double trigger
    % waveforms, loads the ispikes file and uses the unadjusted FramesVector
    unit_type=obj.data.setNames{1}(end);
    if strcmpi(unit_type,'m')
        iswd=pwd;
        is=ispikes('auto');
        Spiketrain=is.data.trial.cluster.spikes;
        FramePoints=stimInfo.data.framePoints/(stimInfo.data.catInfo.samplingrate/1000);
        cd ../..; cd sort/FD
        gname=[is.data.sessionname 'g' is.data.groupname];
        load([gname 'WaveForms']);
        [DoubleTriggerTimes]=RemoveDoubleTriggers(WF,Spiketrain,.4);
        cd(iswd)
        Spiketrain(DoubleTriggerTimes)=[];
    end

    %%% Create the spike count matrix %%%%
    frameCounts = histcie(Spiketrain,FramePoints,'DropLast');

    %%% Create Sparseness Value %%%%
    if Args.Repetitions
        spike_matrix = reshape(frameCounts,numFrames,[])';
        if length(frameCounts) == numFrames % Some Sessions only have 1 Repetition
            mscframes = spike_matrix;
        else
            mscframes = mean(spike_matrix);
        end
        if(Args.Rate)
            % convert counts to spike rate or spikes per second
            bin=(FramePoints(2)-FramePoints(1))/1000;
            Values = mscframes/bin;
        else
            Values = mscframes;
        end
        S = (1-((sum(mscframes)/length(mscframes))^2)/sum((mscframes.*mscframes)/length(mscframes)))/(1-(1/length(mscframes)))*100;
    else
        if Args.Rate
            bin=(FramePoints(2)-FramePoints(1))/1000;
            frameCounts = frameCounts/bin;
        end
        Values = frameCounts;
        S = (1-((sum(frameCounts)/length(frameCounts))^2)/sum((frameCounts.*frameCounts)/length(frameCounts)))/(1-(1/length(frameCounts)))*100;
    end
end
