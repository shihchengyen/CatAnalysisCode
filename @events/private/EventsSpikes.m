function [events] = EventsSpikes(raster,events)
% Take the event line and get the average spike times and std of the spike times for each event.
% The threshold for percentage of reps with a spike.
perc_station = events.NumRepetitions*(events.PercentageReps/100); 
event_index=[];
for e = 1:length(events.StartEvent)    
    event_spike_times=[]; % A vector of spike times for all spikes not including repetitions without any spikes.
    event_spike_counts=[]; % A vector of spike counts for each repetition, including repetitions without any spikes.
    for eve = 1:events.NumRepetitions                
        s_ind = find(raster{eve}>=events.StartEvent(e) & raster{eve}<events.EndEvent(e));        
        e_train = raster{eve};    
        if size(e_train,2) == 1;
            e_train = e_train';
        end
        if isempty(s_ind)
            event_spike_counts = [event_spike_counts 0];
        else
            event_spike_times = [event_spike_times e_train(s_ind)];
            event_spike_counts = [event_spike_counts length(e_train(s_ind))];
        end        
    end
    % This value will determine if an event is further considered.
    stationarity_index = length(find(event_spike_counts > 0));    
    if stationarity_index < perc_station        
        %Spike count calculations
        events.EventMeanSpikeCounts(e)=0;
        events.EventFanoFactor(e)=0;
        %Spike Time calculations
        events.EventMeanSpikeTimes(e)=0;
        event_index = [event_index e];
    else
        %Spike Count calculations
        events.EventMeanSpikeCounts(e) = mean(event_spike_counts);
        events.EventFanoFactor(e) = var(event_spike_counts)/mean(event_spike_counts);           
        %Spike Time calculations
        events.EventMeanSpikeTimes(e) = mean(event_spike_times);       
    end
end
%Remove events without the appropriate number of spikes for a given number
%of repetitions
events.StartEvent(event_index)=[];
events.EndEvent(event_index)=[];
events.EventMeanSpikeCounts(event_index)=[];
events.EventFanoFactor(event_index)=[];
events.EventMeanSpikeTimes(event_index)=[];
events.EventMeanSpikeCounts=events.EventMeanSpikeCounts';
events.EventFanoFactor=events.EventFanoFactor';
events.EventMeanSpikeTimes=events.EventMeanSpikeTimes';

if length(events.StartEvent) ~= length(events.EndEvent)
    ERROR('Error in Positive or Negative Event Limits in Eventfinder Program', '/n')
    return
end
