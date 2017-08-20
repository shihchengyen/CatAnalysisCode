function [events] = EventsStartingPoints(raster,events);
%%%% This function is called during the events analysis routine to find the
%%%% spikes located within the first bin that has a value of 75 percent of
%%%% the maximum psth value.
slide_duration = events.SlideDuration;
nreps = events.NumRepetitions;
for e = 1:length(events.StartEvent)    
    psth_vector = events.StartEvent(e):events.EndEvent(e);
    if psth_vector(1) == 0
        psth_vector(1) = [];
    end
    psth_values = events.PSTH(psth_vector);
    %%% Add the max psth and mean psth values in an event to look at the
    %%% heterogeneity of the values across events.
    events.MaxPSTHValue(e) = max(psth_values); 
    events.MeanPSTHValue(e) = mean(psth_values);
    max_psth_value = max(psth_values);    
    threshold_value = max_psth_value*.75;    
    bin_index = find(psth_values>=threshold_value);    
    event_bin = psth_vector(bin_index(1));    
    event_spike_times=[];
    for r = 1:length(raster)        
        s_ind = find(raster{r}>=event_bin-slide_duration & raster{r}<event_bin+slide_duration);        
        e_train = raster{r};        
        if ~isempty(s_ind)
            if size(e_train,2) == 1
                e_train = e_train';
            end
            event_spike_times = [event_spike_times e_train(s_ind)];
        end        
    end    
    if isempty(event_spike_times)
        events.EventStartingPoints(e) = event_bin;
    else
        events.EventStartingPoints(e) = mean(event_spike_times);
    end    
end

events.MaxPSTHValue=events.MaxPSTHValue';
events.MeanPSTHValue=events.MeanPSTHValue';
events.EventStartingPoints=events.EventStartingPoints';

%%%%%%%%% Find the event intervals and substimulus event intervals
if length(events.EventStartingPoints) > 1
    events.EventIntervalDurations = diff(events.EventStartingPoints);
    sub_ind = find(events.EventIntervalDurations < events.FrameDuration);
else
    events.EventIntervalDurations = [];
    sub_ind = [];
end
if ~isempty(sub_ind)
    events.EventSubStim = events.EventIntervalDurations(sub_ind);
    events.EventSubStimIndex = sub_ind;
else
    events.EventSubStim = [];
    events.EventSubStimIndex = [];
end

% find the percentage of repetitions that have spikes in both events
% making up a sub-stimulus interval
% get number of sub-stimulus intervals
nsubstim = length(sub_ind);
substimreps = zeros(nsubstim,1);
for i = 1:nsubstim
	% get index of sub-stimulus interval
	% this will be the interval between the subi event and the subi+1
	% event
	subi = sub_ind(i);
	subi2 = subi + 1;
	% get start and end points for the first event in a sub-stimulus
	% interval
	estart = events.StartEvent(subi) - slide_duration;
	eend = events.EndEvent(subi) + slide_duration;
	% get start and end points for the second event
	estart2 = events.StartEvent(subi2) - slide_duration;
	eend2 = events.EndEvent(subi2) + slide_duration;
	% get reps that have spikes in subi event
	for r = 1:nreps
		% get the raster for repetition r
		e_train = raster{r};
		% find spikes between start and end points of 1st event in
		% repetition r
		s_ind = find((e_train>=estart) & (e_train<eend));
		if(~isempty(s_ind))
			% check the 2nd event since there were spikes in the
			% first event
			s_ind = find((e_train>=estart2) & (e_train<eend2));
			if(~isempty(s_ind))
				% this rep had spikes in both events so increment
				% substimreps
				substimreps(i) = substimreps(i) + 1;
			end
		end
	end
end
% divide by nreps to get percentage
events.SubStimReps = substimreps / nreps;



