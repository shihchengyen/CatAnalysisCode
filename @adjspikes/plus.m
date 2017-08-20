function r = plus(p,q,varargin)
%@adjspikes/plus Overloaded plus function for adjspikes objects.
%   R = plus(P,Q) combines adjspikes objects P and Q and returns the
%   adjspikes object R.

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
		% r.data.dlist = [p.data.dlist; q.data.dlist];
        r.data.adjSpiketrain = concatenate(p.data.adjSpiketrain , q.data.adjSpiketrain);
        r.data.adjFramePoints = concatenate(p.data.adjFramePoints , q.data.adjFramePoints);
        pts = isfield(p.data,'TossedSpikes');
        qts = isfield(q.data,'TossedSpikes');
        if(pts && qts)
            r.data.TossedSpikes = concatenate(p.data.TossedSpikes , q.data.TossedSpikes);
        elseif(pts && ~qts)
            r.data.TossedSpikes = concatenate(p.data.TossedSpikes, []);
        elseif(~pts && qts)
            r.data.TossedSpikes = concatenate([],q.data.TossedSpikes);
        else
            r.data.TossedSpikes = concatenate([],[]);
        end
        r.data.raster = concatenate(p.data.raster,q.data.raster,'Columnwise');
        r.data.BurstRaster = concatenate(p.data.BurstRaster,q.data.BurstRaster,'ColumnWise');
        r.data.rasterIndex = [p.data.rasterIndex; (p.data.rasterIndex(end) ...
        	+ q.data.rasterIndex(2:end))];
        r.nptdata = plus(p.nptdata,q.nptdata);
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
