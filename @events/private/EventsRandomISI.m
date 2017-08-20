function [randomized_rasters] = EventsRandomISI(spike_raster,rep_duration,num_random)

% The function will read in the rasters and create a designated number of
% randomized rasters with the same ISI distribution.
random_dist=[];
spike_dist=[];
if length(spike_raster) < 10    
    num_extra_rasters = 10 - length(spike_raster);    
    for rr = length(spike_raster)+1:10        
        rp = randperm(length(spike_raster));        
        spike_raster{rr} = spike_raster{rp(1)};        
    end    
end
for r = 1:num_random    
    for i_r = 1:length(spike_raster)        
        if size(spike_raster{1},1) == 1
            isispikes = spike_raster{i_r};
        else
            isispikes = spike_raster{i_r}';
        end        
        if isempty(isispikes)            
            diff_spikes = [];            
        elseif length(isispikes) == 1            
            diff_spikes = [isispikes(1) (rep_duration-isispikes(end))];
            num_isi = 2;            
        else            
            diff_spikes = [isispikes(1) diff(isispikes) (rep_duration-isispikes(end))];
            num_isi = length(diff_spikes);            
        end        
        if isempty(diff_spikes)
            num_isi = 1;
            random_isi{i_r} = spike_raster{i_r};            
        else            
            rp = randperm(num_isi);            
            rand_isi = diff_spikes(rp);            
            if sort(rand_isi)==sort(diff_spikes)
                diff_spikes = [];
                rp = [];
            else
                fprintf('Error in randomization')
                return
            end            
            cum_rand = cumsum(rand_isi);
            random_isi{i_r} = cum_rand(1:end-1);              
        end
    end        
    for i = 1:length(random_isi) 
        if i>1
            raster = concatenate(raster,random_isi{i});
        else
            raster = random_isi{i};
        end
    end
    randomized_rasters{r}.raster = raster;        
end %for random

