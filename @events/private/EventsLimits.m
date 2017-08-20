function [events] = EventsLimits(raster,events,convraster)
% This is the first pass function which takes the raster as an input, then generates
% thresholds and a sliding psth, and will then find the initial event bounderies.
%%%%%%%% Obtain the sliding PSTH %%%%%%%%%%
[PSTH,PSTHBins] = slidingHist(convraster,events.SlideDuration,events.Binsize,events.RepetitionDuration);
% Change to spikes per second
events.PSTH = ((PSTH/events.NumRepetitions)/(events.Binsize/1000))';
events.PSTHBins = PSTHBins';

%%%%%%%% Thresholding the psth based on randomized spiketrains or the Percentile of the PSTH %%%%%
if strcmp(events.ThresholdType,'randomize')   
    num_random_spiketrains = 100;
    [randomized_rasters] = EventsRandomISI(raster,events.RepetitionDuration,num_random_spiketrains);    
    [Thresholds] = EventsRandomThresholds(events,randomized_rasters);    
    upperthreshold = Thresholds.meanupperthreshold;
    lowerthreshold = Thresholds.meanlowerthreshold;    
    Thresholds.NumRandomSpiketrains = num_random_spiketrains;
    Thresholds.ThresholdUsed = 'Mean Upper Threshold via Randomization';    
else %% Threshold based on the Percentile in the nonzero firing rate distribution %%
    prcpsth = events.PSTH(find(events.PSTH>0)); % Take out the 0 values 
    Thresholds.meanupperthreshold = prctile(prcpsth,events.Upper);  
    if isnumeric(events.Lower)
        Thresholds.meanlowerthreshold = prctile(prcpsth,events.Lower);
    else
        Thresholds.meanlowerthreshold = mean(prcpsth);
    end
    Thresholds.ThresholdUsed = 'Percentile Threshold';    
    upperthreshold = Thresholds.meanupperthreshold;
    lowerthreshold = Thresholds.meanlowerthreshold; 
    if Thresholds.meanupperthreshold<Thresholds.meanlowerthreshold 
        upperthreshold = prctile(events.PSTH,99);
        lowerthreshold = Thresholds.meanlowerthreshold;
        Thresholds.meanupperthreshold = prctile(events.PSTH,99);
        fprintf('Upper Threshold is Lower than the Lower Threshold')
    elseif Thresholds.meanupperthreshold==0 & Thresholds.meanlowerthreshold==0
        upperthreshold = prctile(events.PSTH,99);
        lowerthreshold = prctile(events.PSTH,99);
        Thresholds.meanupperthreshold = upperthreshold;
        Thresholds.meanlowerthreshold = lowerthreshold;
        fprintf('Upper Threshold and Lower threshold were zero so using the 99th percentile for both instead.\n')
    end               
end

% EVENT FINDING PART based on the thresholds assigned from the Randomized PSTH's.
limits = [events.PSTH' >= upperthreshold];
diff_limits = diff([0 limits]);
start_event = (find(diff_limits == 1))-events.SlideDuration; % Gives me the start of an event. Shift by slide duration due to diff function
end_event = (find(diff_limits == -1))-events.SlideDuration; % Gives me the end of an event.
%%%% Check the starting and ending event points.
if isempty(start_event);    
    start_event = [];
    end_event = [];
elseif length(start_event) > length(end_event)
    if start_event(end) > end_event(end)
        end_event = [end_event events.PSTHBins(end)];
    end    
end 
%%%%% if there is an error write it out.
if length(start_event) ~= length(end_event)
    error('Error in Positive or Negative Event Limits in EventLimits Function')
end
%%%%% Create the events structure
events.StartEvent = start_event';
events.EndEvent = end_event';
events.Thresholds = Thresholds;