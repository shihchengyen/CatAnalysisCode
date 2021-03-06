function r = plus(p,q,varargin)
%@responsiveness/plus Overloaded plus function for RESPONSIVENESS objects.
%   R = plus(P,Q) combines sparsity objects P and Q and returns the
%   RESPONSIVENESS object R.

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
        r.data.responsiveness = [p.data.responsiveness; q.data.responsiveness];
		r.data.corrcoef = [p.data.corrcoef; q.data.corrcoef];
		r.data.corrp = [p.data.corrp; q.data.corrp];
		r.data.ctrlcoeff = [p.data.ctrlcoeff; q.data.ctrlcoeff];
		r.data.ctrlp = [p.data.ctrlp; q.data.ctrlp];
		r.data.distdiffr = [p.data.distdiffr; q.data.distdiffr];
		r.data.distdiffp = [p.data.distdiffp; q.data.distdiffp];
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
    end
end

