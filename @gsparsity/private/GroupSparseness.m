function [S,Values] = GroupSparseness(Args,spikes,FramePoints,stimInfo)
% The function will calculate the sparseness of the firing based on the
% psth at the frame resolution or a given binsize, using the vinje 
% and gallant sparseness equation.

if Args.Synchrony % Find the Synchronous spiking across cells
    if Args.Repetitions
        numFrames = stimInfo.data.catInfo.num_frames;
        numRepetitions = stimInfo.data.catInfo.num_repetitions;
        for ii = 1:numRepetitions
            RepFrames = FramePoints(ii*numFrames-numFrames+1:(ii*numFrames)+1);
            RepCounts = histcie(spikes,RepFrames(1):Args.Binsize:RepFrames(end));
            RepCounts(find(RepCounts>0)) = 1;
            RepCounts = prod(RepCounts,2);
            if ii > 1
                FrameCounts = concatenate(FrameCounts,RepCounts,'ColumnWise');
            else
                FrameCounts = RepCounts;
            end
        end
    else
        FrameCounts = histcie(spikes,FramePoints(1):Args.Binsize:FramePoints(end),'DropLast');
        FrameCounts(find(FrameCounts>0))=1;
        FrameCounts = prod(FrameCounts,2);
    end       
elseif Args.Sum % Add the Spike Counts so that the maximum value per bin is equal to the number of cells.
    FrameCounts = histcie(spikes,FramePoints,'DropLast');
    if Args.SpikeCounts
        if Args.Normalize % Normalize each Spike Count by the mean spike count then multiply each cell.
            FrameCounts = FrameCounts./repmat(max(FrameCounts),size(FrameCounts,1),1);
        end       
        FrameCounts = sum(FrameCounts,2);        
    else
        FrameCounts(find(FrameCounts>0)) = 1;
        FrameCounts = sum(FrameCounts,2);
    end
end

%%% Create Sparseness Value %%%%
if Args.Repetitions
    numFrames = stimInfo.data.catInfo.num_frames;
    if Args.Synchrony
        spike_matrix = FrameCounts';
    else
        spike_matrix = reshape(FrameCounts,numFrames,[])';
    end
    mscframes = mean(spike_matrix);
    Values = mscframes';
    S = (1-((sum(mscframes)/length(mscframes))^2)/sum((mscframes.*mscframes)/length(mscframes)))/(1-(1/length(mscframes)))*100;
else
    Values = FrameCounts';
    S = (1-((sum(FrameCounts)/length(FrameCounts))^2)/sum((FrameCounts.*FrameCounts)/length(FrameCounts)))/(1-(1/length(FrameCounts)))*100;
end