function [counts,meancounts,spike_matrix,timebins,binsize] = CalcFiringRate(obj,stimInfo,Args)

% The function will calculate the firing rates based on the specfied binsize
% or the frame duration. The Firing rate can be calculated with an Overlap
% value or an adapted binsize depending on the desired average spike rate.

binsize = Args.Binsize;
numFrames = stimInfo.data.catInfo.num_frames;
numRepetitions = stimInfo.data.catInfo.num_repetitions;
FrameDuration = stimInfo.data.catInfo.frame_duration;

% Load the SpikeTrain and FrameVector, adjusted or non-adjusted
if Args.Repetitions
    Spiketrain = obj.data.adjSpiketrain;
    FramePoints = obj.data.adjFramePoints;
else % for the ispikes object
    Spiketrain = obj.data.trial.cluster.spikes;
    FramePoints = stimInfo.data.framePoints/(stimInfo.data.catInfo.samplingrate/1000);
end

if Args.Count % Binsize based on the specified average spike number, typically 1, 2 or 10 spikes
    mean_sc = length(Spiketrain)/(FramePoints(end)/1000);
    binsize = (Args.Count/mean_sc)*1000; % Get the binsize for a desired rate in spikes/second
end

if ~isnumeric(binsize)  %% If a binsize is not avaliable use the frame duration %%%%

    binsize = FrameDuration;
    
    if Args.Overlap %% If a specified overlap is desired, the output will be a single vector

        m = round(binsize/Args.Overlap);
        shape = ones(1,m);
        edges = FramePoints(1):Args.Overlap:FramePoints(end);
        frameCounts = sum(histc(Spiketrain,edges,2),1);
        spike_matrix = conv(frameCounts,shape)';
        % get length of spike_matrix
        smlength = length(spike_matrix);
        % get length of edge effects
        trimlength = m - 1;
        % remove incomplete bins at the beginning of spike_matrix
        spike_matrix(1:trimlength)=[];
        % remove incomplete bins at the end of spike_matrix
        spike_matrix((smlength-trimlength+1):smlength) = [];
        counts = spike_matrix;%/(binsize/1000);
        timebins = [];
        meancounts = [];
      
    elseif Args.Repetitions
        
        frameCounts = histcie(Spiketrain,FramePoints,'DropLast');
        spike_matrix = reshape(frameCounts,numFrames,[])';
        % need this to convert spike counts to rate
        counts = (mean(spike_matrix))';
        % subtract first frame point to make sure we start from 0
        timebins = (FramePoints(1:numFrames) - FramePoints(1))';
        if(Args.NoiseCorr)
            % subtract the mean count from each bin and transpose so that
            % when we convert the matrix to a vector, it is in the right
            % order, i.e. frames of rep 1 followed by frames of rep 2
            tmpcounts = (spike_matrix - repmat(counts',numRepetitions,1))';
            meancounts = tmpcounts(:);
        else
            meancounts = counts;
        end
        if Args.Rate
            rateconv = binsize/1000;
            counts = counts/rateconv;
            meancounts = meancounts/rateconv;
            spike_matrix = spike_matrix/rateconv;
        end
        
    else

        frameCounts = histcie(Spiketrain,FramePoints,'DropLast');
        spike_matrix = reshape(frameCounts,numFrames,[])';
        % need this to convert spike counts to rate
        counts = (mean(spike_matrix))';
        % subtract first frame point to make sure we start from 0
        timebins = (FramePoints(1:numFrames) - FramePoints(1))';
        if(Args.NoiseCorr)
            % subtract the mean count from each bin and transpose so that
            % when we convert the matrix to a vector, it is in the right
            % order, i.e. frames of rep 1 followed by frames of rep 2
            tmpcounts = (spike_matrix - repmat(counts',numRepetitions,1))';
            meancounts = tmpcounts(:);
        else
            meancounts = counts;
        end
        if Args.Rate
            rateconv = binsize/1000;
            counts = counts/rateconv;
            meancounts = meancounts/rateconv;
            spike_matrix = spike_matrix/rateconv;
        end
    end

else

    if(Args.Overlap>0) % If a specified overlap is desired, the output will be a single vector

        m = round(binsize/Args.Overlap);
        shape = ones(m,1);
        % edges = FramePoints(1):Args.Overlap:FramePoints(end);
        % FR = sum(histcie(Spiketrain,edges,'DropLast'),1);
        % FR = conv(FR,shape)';
        % FR(1:m-1)=[];
        % FR(end-m:end)=[];
        % counts = FR;
        % timebins = [];
        % meancounts = [];
        overlap = Args.Overlap;
        %%% Create the spike count matrix for specified bin size %%%%
        % get number of bin sizes in each repetition by dividing the
        % duration of the first repetition by bin size and then taking the
        % floor so that we don't include spikes from the next frame in the last
        % bin
        rep1 = numRepetitions - 1;
        % create vector of time of first frame of each repetition
        repframe = FramePoints(1:numFrames:(numFrames*rep1+1));
        % find the smallest interval to compute number of bins otherwise
        % we could end up creating non-monotonically increasing bins
        numbins = floor( min(diff(repframe))/overlap );
        % add 1 to numbins since we start from 0
        nbins1 = numbins + 1;
        % create second matrix in matrix multiplication which contains a column of
        % ones, followed by a column increasing from 0 to duration in steps of
        % binsize.
        m1 = [ones(nbins1,1) (0:overlap:(numbins*overlap))'];
        % create first matrix in matrix multiplication which is the start time
        % of each repetition in the first row and ones in the second row
        m2 = [repframe; ones(1,numRepetitions)];
        % create histcie-edges using matrix multiplication
        binEdges = m1 * m2;
        % get timebins from the binEdges for the 1st repetition
        % subtract the first value since it might not be 0
        timebins = binEdges(:,1) - binEdges(1,1);
        % do histcie using all the binEdges. This will add 1 bin in each
        % repetition, which is the remainder from dividing the repetition
        % duration with binsize. Even if duration divides evenly into binsize,
        % there will be an extra bin since the next repetition always starts
        % from the end of the previous repetition. The histcie function is able
        % to handle bin edges that are equal so this works out well. We will
        % keep the last point since that will allow us to drop the last row
        % after we reshape. Otherwise, we will be missing the interval due to
        % dividing the last repetition by binsize
        frameCounts = histcie(Spiketrain,binEdges(:));
        % reshape into numbins+1 x num_repetitions since the last bin will be
        % remainder from dividing repetition duration by bin size.
        sm = reshape(frameCounts,nbins1,[]);
        % drop the last row since it is an incomplete window
        spike_matrix = sm(1:numbins,:);
        counts1 = mean(spike_matrix,2);
        % convolve to get correct window size
        counts2 = conv(counts1,shape);
        % drop the edges to get the correct number of bins
        counts = counts2(m:end);
        if(Args.NoiseCorr)
            % subtract the mean count from each bin and transpose so that
            % when we convert the matrix to a vector, it is in the right
            % order, i.e. frames of rep 1 followed by frames of rep 2
            tc = spike_matrix - repmat(counts1,1,numRepetitions);
            % pad tmpcounts with rows of 0's so we can safely convert
            % tmpcounts to a single vector and take the convolution without
            % mixing spike counts from different repetitions
            tc2 = [tc; zeros(m-1,numRepetitions)];
            tc3 = conv(tc2(:),shape);
            % reshape back to repetitions
            tc4 = reshape(tc3(m:end),[],100);
            % get rid of the bins resulting from the zero padding as well
            % as the bins in front of the first full window
            tc5 = tc4(1:numbins,:);
            meancounts = tc5(:);
            % take the transpose of spike_matrix so we can compute the
            % correlation coefficients of each column
            spike_matrix = spike_matrix';
        else
            % add zero to the end so that when we use stairs to plot the firing
            % rate, the last data point will stretch to the end of the last bin
            meancounts = [counts; 0];
        end
        if(Args.Rate)
            rateconv = binsize / 1000;
            counts = counts / rateconv;
            meancounts = meancounts / rateconv;
            spike_matrix = spike_matrix / rateconv;
        end
        
    elseif Args.Count % If a specified SpikeRate is desired, the output will be a single vector

        counts=[];
        timebins=[];
        for r = 1:numRepetitions
            points = FramePoints(r*numFrames-numFrames+1:(r*numFrames)+1);
            edges = points(1):Args.Binsize:points(end)';
            spikeCounts = histcie(Spiketrain,edges,'DropLast')';
            counts = [counts spikeCounts];
            timebins = [timebins diff(edges)];
        end
        if Args.Rate
            counts = counts./(timebins/1000);
        end
        meancounts = counts;

    elseif Args.Repetitions

        for r = 1:numRepetitions
            points = FramePoints(r*numFrames-numFrames+1:(r*numFrames)+1);
            edges = points(1):Args.Binsize:points(end); % Due to a variable binsize data may be lost from the last bin
            spikeCounts = (histcie(Spiketrain,edges))'; % leave the last bin so counts will match up with edges
            if r>1
                counts = concatenate(counts,spikeCounts);
                % timebins = concatenate(timebins,diff(edges));
            else
                counts = spikeCounts;
                timebins = edges;
            end
        end

        if Args.Rate
            counts = counts./(round(mean(timebins(:)))/1000);
        end
        counts = (mean(counts))'; % Repetitions argument converts the spikecounts across the Repetitions into a PSTH
        meancounts = counts;
        % timebins = (mean(timebins))';

    else % Firingrate calculation across the entire spiketrain, but first calculated on the repetitions

        if ~exist('numRepetitions')
            counts=[];
            timebins=[];
            for r = 1:numRepetitions
                points = FramePoints(r*numFrames-numFrames+1:(r*numFrames)+1);
                edges = points(1):Args.Binsize:points(end)';
                spikeCounts = histcie(Spiketrain,edges,'DropLast')';
                counts = [counts spikeCounts];
                timebins = [timebins diff(edges)];
            end
            if Args.Rate
                counts = (counts./(round(mean(timebins(:)))/1000))';
            end
            meancounts = counts;
            timebins = timebins';
        else
            timebins = (FramePoints(1):Args.Binsize:FramePoints(end))';
            counts = histcie(Spiketrain,timebins,'DropLast');
            meancounts = counts;
            if Args.Rate
                counts = counts./(round(mean(diff(timebins)))/1000);
            end
        end
        spike_matrix=counts;
    end
end
