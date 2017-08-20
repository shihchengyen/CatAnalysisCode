function r = plus(p,q,varargin)
%@firingrate/plus Overloaded plus function for firingrate objects.
%   R = plus(P,Q) combines sparsity objects P and Q and returns the
%   firingrate object R.

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
		r.data.firingRate = concatenate(p.data.firingRate,q.data.firingRate,'Columnwise');
		r.data.meancounts = concatenate(p.data.meancounts,q.data.meancounts,'Columnwise');
		r.data.timebins = concatenate(p.data.timebins,q.data.timebins,'Columnwise');
		r.data.binsize = [p.data.binsize; q.data.binsize];
        r.data.Fit = [p.data.Fit,q.data.Fit];
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
    end
end