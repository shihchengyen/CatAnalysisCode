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
		% useful fields for most objects
		r.data.numSets = p.data.numSets + q.data.numSets;
		r.data.setNames = {p.data.setNames{:} q.data.setNames{:}};
		
		% object specific fields
        r.data.ncells = [p.data.ncells; q.data.ncells];        
		r.data.mediancorr = [p.data.mediancorr; q.data.mediancorr];
		r.data.lqcorr = [p.data.lqcorr; q.data.lqcorr];
		r.data.uqcorr = [p.data.uqcorr; q.data.uqcorr];
		r.data.medianvang = [p.data.medianvang; q.data.medianvang];
		r.data.lqvang = [p.data.lqvang; q.data.lqvang];
		r.data.uqvang = [p.data.uqvang; q.data.uqvang];
        r.data.medianjep = [p.data.medianjep; q.data.medianjep];
		r.data.lqjep = [p.data.lqjep; q.data.lqjep];
		r.data.uqjep = [p.data.uqjep; q.data.uqjep];
        
		r.data.pszscore = concatenate(p.data.pszscore,q.data.pszscore, ...
            'Columnwise');
		r.data.kwp = [p.data.kwp; q.data.kwp];
        r.data.mfr = concatenate(p.data.mfr,q.data.mfr,'Columnwise');

		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
