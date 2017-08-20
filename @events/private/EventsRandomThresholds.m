function [thresholds] = EventsRandomThresholds(events,randomized_rasters)
% The function requires the eventsrandomisi function, the randomized
% spiketrains in random, the bin to be used for creating the psth's
% required for the threshold to be used in the event analysis.
for r = 1:length(randomized_rasters)   
    [data] = slidingHist(randomized_rasters{r}.raster,events.SlideDuration,events.Binsize,events.RepetitionDuration);    
    % Since the edges can contain many spikes due to short intervals I take out the edges.
    max_psth(r) = max(data(100:length(data)-100));
    mean_psth(r) = mean(data(100:length(data)-100));
end
thresholds.meanupperthreshold = mean(max_psth);
thresholds.maxupperthreshold = max(max_psth);
thresholds.minupperthreshold = min(max_psth);
thresholds.meanlowerthreshold = mean(mean_psth);
thresholds.maxlowerthreshold = max(mean_psth); 
