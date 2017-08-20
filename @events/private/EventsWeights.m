function [events] = EventsWeights(events);

% This function will load the end and start points of an event and assign
% a weight value depending on the integral under the psth for that
% particular event. Then it will weight the event based on its duration and
% area.
% This function runs under the eventsanalysis function.

warning off MATLAB:divideByZero
for e = 1:length(events.start_event)    
    p_values = events.data_psth(events.start_event(e):events.end_event(e));    
    event_durations(e) = length(events.start_event(e):events.end_event(e));    
    area_values(e) = trapz(p_values);    
    event_weights(e) = area_values(e)/event_durations(e);    
end
events.event_weights = round(event_weights/min(event_weights));
