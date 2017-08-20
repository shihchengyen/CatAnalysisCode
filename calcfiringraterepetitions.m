function [FR] = calcfiringraterepetitions(obj,stimInfo,Args);
% The function will calculate the firing rate values based on the input
% arguments for each repetition.

Spiketrain = obj.data.adjSpiketrain;
FramePoints = obj.data.adjFramePoints;
Raster = obj.data.raster;
numFrames = stimInfo.data.catInfo.num_frames;
numRepetitions = stimInfo.data.catInfo.num_repetitions;
FrameDuration = stimInfo.data.catInfo.frame_duration;

if strcmpi(Args.Binsize,'frame')
    Args.Binsize=FrameDuration;
end

if Args.Overlap % If a specified overlap is desired, the output will be a single vector
    
    m = round(Args.Binsize/Args.Overlap);
    shape = ones(1,m);
    for r = 1:numRepetitions
        points = FramePoints(r*numFrames-numFrames+1:(r*numFrames)+1);
        edges = points(1):Args.Overlap:points(end);
        spikeCounts = sum(histcie(Spiketrain,edges,'DropLast'),1);
        spikeCounts = conv(spikeCounts,shape);
        spikeCounts(1:m-1)=[];spikeCounts(end-m:end)=[];
        if r>1
            FR = concatenate(FR,spikeCounts);
        else
            FR = spikeCounts;
        end
    end
    
elseif Args.Count % If a specified SpikeRate is desired, the output will be a single vector
    
    mean_sc = length(Spiketrain)/(FramePoints(end)/1000);
    Args.Binsize = (Args.Count/mean_sc)*1000; % Get the binsize for a desired rate in spikes/second
    for r = 1:numRepetitions
        points = FramePoints(r*numFrames-numFrames+1:(r*numFrames)+1);
        edges = points(1):Args.Binsize:points(end);
        spikeCounts = histcie(Spiketrain,edges,'DropLast')';
        if r>1
            FR = concatenate(FR,spikeCounts);
        else
            FR = spikeCounts;
        end
    end
    
  
else % only if Binsize is an input Argument or for the Repetitions Argument
    
    for r = 1:numRepetitions
        points = FramePoints(r*numFrames-numFrames+1:(r*numFrames)+1);
        edges = points(1):Args.Binsize:points(end); % Due to a variable binsize data may be lost from the last bin
        spikeCounts = (histcie(Spiketrain,edges,'DropLast'))'; % take out the last bin 
        if r>1
            FR = concatenate(FR,spikeCounts);
            timebins = concatenate(timebins,diff(edges));
        else
            FR = spikeCounts;
            timebins = diff(edges);
        end
    end         
    
    if Args.Rate
        FR = FR./(round(mean(timebins(:)))/1000);
    end
    
end