function [events] = EventsLimitsConfirm(events);
% This function will look at the low points between start and end points of
% adjacent events to see if they drop below the lower threshold. If they do
% not, then the two events will be merged together.
data_psth = events.PSTH';
start_event = events.StartEvent';
end_event = events.EndEvent';
thresholds = events.Thresholds;
if length(start_event) == 1    
    % Do nothing since there is only one event      
else
    for v = 1:length(start_event)-1       
        valley = min(data_psth(end_event(v):start_event(v+1)));        
        if valley > thresholds.meanlowerthreshold            
            end_event(v) = -1;
            start_event(v+1) = -1;
        end              
    end     
end %if
start_event(find(start_event<0)) = [];
end_event(find(end_event<0)) = [];
if length(start_event) ~= length(end_event)
    error('Error in Positive or Negative Event Limits in EventsLimitsConfirm Function')
end

%%%%%%%%%%% Now Widen the start and end points until the lower threshold is
%%%%%%%%%%% hit
% This function will widen the boundries and make the start and end event points at the
% lowerthreshold level.
for e = 1:length(start_event)
    if start_event(e) > 1
        valley = data_psth(start_event(e)-1); % look at preceding psth bin 
        neg_move = 0;
        while valley > thresholds.meanlowerthreshold
            neg_move = neg_move+1;
            if (start_event(e)-neg_move) == events.PSTHBins(1) %% Check to see if you hit the stimulus Beginning
                valley = thresholds.meanlowerthreshold;
            else
                valley = data_psth(start_event(e)-neg_move);
            end
            if data_psth(start_event(e)) <= thresholds.meanlowerthreshold
                neg_move=0;
                valley = thresholds.meanlowerthreshold;
            end
        end         
        start_event(e) = start_event(e)-(neg_move);    
    end    
    if end_event(e) < events.PSTHBins(end)
        valley = data_psth(end_event(e)+1); % look at next psth bin 
        pos_move = 0;
        while valley > thresholds.meanlowerthreshold
            pos_move = pos_move+1;
            if (end_event(e)+pos_move) == events.PSTHBins(end) %% Check to see if you hit the stimulus Beginning
                valley = thresholds.meanlowerthreshold;
            else
                valley = data_psth(end_event(e)+pos_move);
            end
            if data_psth(end_event(e)) <= thresholds.meanlowerthreshold
                pos_move=0;
                valley = thresholds.meanlowerthreshold;
            end
        end 
        end_event(e) = end_event(e)+(pos_move); 
    end     
end       

if length(start_event) ~= length(end_event)
    error('Error in Positive or Negative Event Limits in EventsLimitsConfirm Function')
end
%%%%%% Add the info back into the structure
events.StartEvent = start_event'; %% Due to the fact of histogram and binning the data
events.EndEvent = (end_event-1)';
