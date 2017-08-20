function r = plus(p,q,varargin)

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
		% check Args to make sure they are the same
		r.data.binSize = [p.data.binSize; q.data.binSize];
		r.data.WinStartTime = concatenate(p.data.WinStartTime,q.data.WinStartTime,'Columnwise');
		r.data.DataSyncSpikes = concatenate(p.data.DataSyncSpikes,q.data.DataSyncSpikes,'Columnwise');
		r.data.SurrSyncHist = concatenate(p.data.SurrSyncHist,q.data.SurrSyncHist,'Columnwise');
		r.data.DataZScore = concatenate(p.data.DataZScore,q.data.DataZScore,'Columnwise');
		r.data.SurrSyncMean = concatenate(p.data.SurrSyncMean,q.data.SurrSyncMean,'Columnwise');
		r.data.SurrSyncStd = concatenate(p.data.SurrSyncStd,q.data.SurrSyncStd,'Columnwise');
		r.data.WindowIndex = [p.data.WindowIndex; (p.data.WindowIndex(end) ...
			+ q.data.WindowIndex(2:end))];
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
