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
		r.data.SubBinSize = [p.data.SubBinSize; q.data.SubBinSize];
		r.data.ActualSubBinSize = [p.data.ActualSubBinSize; q.data.ActualSubBinSize];
		r.data.BinSums = {p.data.BinSums{:} q.data.BinSums{:}};
		r.data.refreshRate = [p.data.refreshRate; q.data.refreshRate];
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
