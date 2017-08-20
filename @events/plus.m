function r = plus(p,q,varargin)
%@events/plus Overloaded plus function for EVENTS objects.
%   R = plus(P,Q) combines sparsity objects P and Q and returns the
%   EVENTS object R.
% check for empty object
% combine objects
% get name of class
classname = mfilename('class');

% check if first input is the right kind of object
if(~isa(p,classname))
    % check if second input is the right kind of object
    if(~isa(q,classname))
        % both inputs are not the right kind of object so create empty
        % object and return it
        r = feval(classname);
    else
        % second input is the right kind of object so return that
        r = q;
    end
else
    if(~isa(q,classname))
        % p is the right kind of object but q is not so just return p
        r = p;
    elseif(isempty(p))
        % p is right object but is empty so return q, which should be
        % right object
        r = q;
    elseif(isempty(q))
        % p are q are both right objects but q is empty while p is not
        % so return p
        r = p;
    else
        % both p and q are the right kind of objects so add them
        % together
        % assign p to r so that we can be sure we are returning the right
        % object
        r = p;
        r.data.numSets = p.data.numSets + q.data.numSets;
        r.data.setNames = {p.data.setNames{:} q.data.setNames{:}};
        r.data.dlist = [p.data.dlist; q.data.dlist];
        r.data.events.PSTH = concatenate(r.data.events.PSTH,q.data.events.PSTH,'Columnwise');
        r.data.events.PSTHBins = concatenate(r.data.events.PSTHBins,q.data.events.PSTHBins,'Columnwise');
        r.data.events.StartEvent = concatenate(r.data.events.StartEvent,q.data.events.StartEvent,'Columnwise');
        r.data.events.EndEvent = concatenate(r.data.events.EndEvent,q.data.events.EndEvent,'Columnwise');
        r.data.events.EventProbability = concatenate(r.data.events.EventProbability,q.data.events.EventProbability,'Columnwise');
        r.data.events.EventDurations = concatenate(r.data.events.EventDurations,q.data.events.EventDurations,'Columnwise');
        r.data.events.Thresholds = concatenate(r.data.events.Thresholds,q.data.events.Thresholds,'Columnwise');
        r.data.events.EventMeanSpikeCounts = concatenate(r.data.events.EventMeanSpikeCounts,q.data.events.EventMeanSpikeCounts,'Columnwise');
        r.data.events.EventFanoFactor = concatenate(r.data.events.EventFanoFactor,q.data.events.EventFanoFactor,'Columnwise');
        r.data.events.EventMeanSpikeTimes = concatenate(r.data.events.EventMeanSpikeTimes,q.data.events.EventMeanSpikeTimes,'Columnwise');
        r.data.events.MaxPSTHValue = concatenate(r.data.events.MaxPSTHValue,q.data.events.MaxPSTHValue,'Columnwise');
        r.data.events.MeanPSTHValue = concatenate(r.data.events.MeanPSTHValue,q.data.events.MeanPSTHValue,'Columnwise');
        r.data.events.EventStartingPoints = concatenate(r.data.events.EventStartingPoints,q.data.events.EventStartingPoints,'Columnwise');
        r.data.events.FrameDuration = concatenate(r.data.events.FrameDuration,q.data.events.FrameDuration,'Columnwise');
        r.nptdata = plus(r.nptdata,q.nptdata);
    end
end