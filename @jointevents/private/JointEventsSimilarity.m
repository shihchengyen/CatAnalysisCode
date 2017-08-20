function [expected_dist,measured_dist,H,P] = JointEventsSimilarity(rasters,JointEventsInfo,repetition_duration)
%%% This private function will take the two rasters and calculate an
%%% individual repetition event calculation to create two distributions of
%%% the expected overlap by chance and of the measured overlap. These
%%% distributions are compared through a ttest to look for significance.

NumRepetitions = length(rasters(1).raster);
for r = 1:NumRepetitions    
    %%%%%% Get the Event Points for the First Cell %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% PSTH based on the average spike counts per bin in time %%%%%%%
    [y,t] = slidingHist(rasters(1).raster{r},1,JointEventsInfo(1).Binsize,repetition_duration);
    DataPSTH_1 = y/(JointEventsInfo(1).Binsize/1000);    
    bins = 0:1.5:repetition_duration;
    random_spiketrain = histcie(rasters(1).raster{r},bins)';
    randmat = rand(100,length(random_spiketrain));  
    [srand,irand] = sort(randmat');
    randmat = irand';
    random_spiketrain = random_spiketrain(randmat);
    slide_mat = convmtx(ones(JointEventsInfo(1).Binsize),(length(bins)-JointEventsInfo(1).Binsize)+1);
    random_mat = random_spiketrain'.*slide_mat;
    
    JointEventsInfo(1).Thresholds.meanupperthreshold = mean(UpperThresholds);
    JointEventsInfo(1).Thresholds.meanlowerthreshold = mean(LowerThresholds);
    % EVENT FINDING PART based on the thresholds assigned from the Randomized PSTH's.
    limits = [DataPSTH_1 > JointEventsInfo(1).Thresholds.meanupperthreshold];
    diff_limits = diff([0 limits]);
    start_event = find(diff_limits == 1); % Gives me the start of an event.
    end_event = find(diff_limits == -1); % Gives me the end of an event.
    %%%% Check the starting and ending event points.
    if isempty(start_event);    
        start_event = [];
        end_event = [];
    elseif length(start_event) > length(end_event)
        if start_event(end) > end_event(end)
            %fprintf('Initial error in the last end event')
            end_event = [end_event length(DataPSTH_1)];
        end    
    end 
    %%%%% if there is an error write it out.
    if length(start_event) ~= length(end_event)
        fprintf('Error in Positive or Negative Event Limits in EventLimits Function', '/n')
        return
    end
    if ~isempty(start_event)
        events.DataPSTH = DataPSTH_1;
        events.StartEvent = start_event;
        events.EndEvent = end_event;
        events.Thresholds = JointEventsInfo(1).Thresholds;    
        [events] = EventsLimitsConfirm(events);
        [events] = EventsWiden(events);
        [events1] = EventsDurations(events); events=[];
    else
        events1.DataPSTH = DataPSTH_1;
        events1.StartEvent = [];
        events1.EndEvent = [];
        events1.Thresholds = JointEventsInfo(1).Thresholds;
        events1.EventProbability = 0;
        events1.EventDurations = [];
    end
    
    %%%%%% Get the Event Points for the Second Cell %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% PSTH based on the average spike counts per bin in time %%%%%%%
    [y,t] = slidingHist(rasters(2).raster{r},1,JointEventsInfo(2).Binsize,repetition_duration);    
    DataPSTH_2 = y/(JointEventsInfo(2).Binsize/1000); 
    bins = 0:1:repetition_duration;
    random_spiketrain = histcie(rasters(2).raster{r},bins)';
    randmat = rand(100,length(random_spiketrain));  
    [srand,irand] = sort(randmat');
    randmat = irand';
    random_spiketrain = random_spiketrain(randmat);
    for ii = 1:100
        [y,t] = slidingHist(find(random_spiketrain(ii,:)==1),1,JointEventsInfo(2).Binsize,repetition_duration);
        Random_DataPSTH = y/(JointEventsInfo(2).Binsize/1000);
        UpperThresholds(ii) = max(Random_DataPSTH);
        LowerThresholds(ii) = mean(Random_DataPSTH);
    end   
    JointEventsInfo(2).Thresholds.meanupperthreshold = mean(UpperThresholds);
    JointEventsInfo(2).Thresholds.meanlowerthreshold = mean(LowerThresholds);
    % EVENT FINDING PART based on the thresholds assigned from the Randomized PSTH's.
    limits = [DataPSTH_2 > JointEventsInfo(2).Thresholds.meanupperthreshold];
    diff_limits = diff([0 limits]);
    start_event = find(diff_limits == 1); % Gives me the start of an event.
    end_event = find(diff_limits == -1); % Gives me the end of an event.
    %%%% Check the starting and ending event points.
    if isempty(start_event);    
        start_event = [];
        end_event = [];
    elseif length(start_event) > length(end_event)
        if start_event(end) > end_event(end)
            %fprintf('Initial error in the last end event')
            end_event = [end_event length(DataPSTH_1)];
        end    
    end 
    %%%%% if there is an error write it out.
    if length(start_event) ~= length(end_event)
        fprintf('Error in Positive or Negative Event Limits in EventLimits Function', '/n')
        return
    end
    if ~isempty(start_event)
        events.DataPSTH = DataPSTH_1;
        events.StartEvent = start_event;
        events.EndEvent = end_event;
        events.Thresholds = JointEventsInfo(1).Thresholds;    
        [events] = EventsLimitsConfirm(events);
        [events] = EventsWiden(events);
        [events2] = EventsDurations(events); events=[];
    else
        events2.DataPSTH = DataPSTH_1;
        events2.StartEvent = [];
        events2.EndEvent = [];
        events2.Thresholds = JointEventsInfo(1).Thresholds;
        events2.EventProbability = 0;
        events2.EventDurations = [];
    end
    %%%%%%% Expected overlap in the to events based on chance %%%%%%%%%
    expected_dist(r) = events1.EventProbability*events2.EventProbability;
    
    eventline1=[];
    if ~isempty(events1.StartEvent)        
        for e = 1:length(events1.StartEvent)
            eventline1 = [eventline1 events1.StartEvent(e):events1.EndEvent(e)];
        end
    end
    eventline2=[];
    if ~isempty(events2.StartEvent)        
        for e = 1:length(events2.StartEvent)
            eventline2 = [eventline2 events2.StartEvent(e):events2.EndEvent(e)];
        end
    end
    
    if ~isempty(eventline1)
        if ~isempty(eventline2)
            overlap = intersect(eventline1,eventline2);    
            %%% Joint overlap Probability %%%%
            measured_dist(r) = length(overlap)/repetition_duration;
        else
            measured_dist(r) = 0;
        end
    else
        measured_dist(r) = 0;
    end    
    r
end

[H,P] = TTEST2(expected_dist,measured_dist);
