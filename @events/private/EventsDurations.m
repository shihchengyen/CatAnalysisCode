function [events] = EventsDurations(events);

% This function is called during the eventsanalysis functions and will
% calculate the event durations and the amount of time a cell fires within
% an event, the probability of observing a cell firing significantly.

stimulus_duration = ceil(events.RepetitionDuration);
event_durations=[];
for e = 1:length(events.StartEvent)    
    event_durations(e) = length((events.StartEvent(e)):events.EndEvent(e));    
end
events.EventProbability = sum(event_durations)/stimulus_duration;
events.EventDurations = event_durations';
