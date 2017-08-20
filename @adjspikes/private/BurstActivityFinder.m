function [BurstRaster] = BurstActivityFinder(raster,Args)
% This private adjspikes function will calculate the bursting activity within the
% spiketrain based on the intra and inter burst intervals.

for r = 1:size(raster,1);
    spiketrain = raster(r,:);
    spiketrain = spiketrain(find(spiketrain>0));
    isis = diff(spiketrain);
    BurstSpikes = (find(isis>Args.InterBI))+1;
    if ~isempty(BurstSpikes)
        if BurstSpikes(1)~=1
            BurstSpikes = [1 BurstSpikes];
        end
        if BurstSpikes(end) == length(spiketrain)
            BurstSpikes(end) = [];
        end
        new_spikes = [];
        for ii = 1:length(BurstSpikes)
            n_isi = diff(spiketrain(BurstSpikes(ii):BurstSpikes(ii)+1));
            if n_isi > Args.IntraBI
                BurstSpikes(ii) = NaN;
            end
            count = 1;
            while n_isi < Args.IntraBI
                new_spikes = [new_spikes BurstSpikes(ii)+count];
                count = count + 1;
                if (BurstSpikes(ii)+count) > length(spiketrain)
                    n_isi = Args.IntraBI+1;
                else
                    n_isi = diff(spiketrain(BurstSpikes(ii)+(count-1):BurstSpikes(ii)+count)); 
                end
            end
        end    
        BurstSpikes = sort([BurstSpikes new_spikes]); % Add the new spikes
        BurstSpikes = BurstSpikes(find(BurstSpikes>0)); % Take out the NaNs
        BurstSpikes = unique(BurstSpikes); % Possibly doubled the counts.
        if r == 1
            BurstRaster = spiketrain(BurstSpikes);
        else
            BurstRaster = concatenate(BurstRaster,spiketrain(BurstSpikes));
        end
    else
        if r == 1
            BurstRaster(r) = NaN;
        else
            BurstRaster = concatenate(BurstRaster,NaN);
        end    
    end
end