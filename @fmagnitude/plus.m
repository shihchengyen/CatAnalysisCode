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
		% both p and q are the right kind of objects and they are not 
        % empty so add them together
		% assign p to r so that we can be sure we are returning the right
		% object
		r = p;
        if(sum(p.data.f - q.data.f)~=0)
            fprintf('Warning %s: Frequency vectors not the same length!',classname);
        end
		r.data.nCells = p.data.nCells + q.data.nCells;
		r.data.cellname = {p.data.cellname{:} q.data.cellname{:}};
		r.data.mean = [p.data.mean q.data.mean];
		r.data.stdev = [p.data.stdev q.data.stdev];
		% concatenate fftMag
		r.data.fftMag = [p.data.fftMag q.data.fftMag];
		r.data.colIndex = [p.data.colIndex; p.data.colIndex(end) ...
			+ q.data.colIndex(2:end)];
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
