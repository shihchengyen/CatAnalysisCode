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
		if(p.data.binSize~=q.data.binSize)
			fprintf('Warning %s: Bin sizes not the same!\n',classname);
		end
        if(sum(p.data.lags - q.data.lags)~=0)
            fprintf('Warning %s: Lag vectors not the same!\n',classname);
        end
		r.data.xc = concatenate(p.data.xc,q.data.xc,'Columnwise');
		r.data.sxc = concatenate(p.data.sxc,q.data.sxc,'Columnwise');
		r.data.ecf = [p.data.ecf; q.data.ecf];
		r.data.kspvalue = [p.data.kspvalue; q.data.kspvalue];
		r.data.ksstat = [p.data.ksstat; q.data.ksstat];
		r.data.legendstr = {p.data.legendstr{:} q.data.legendstr{:}};
		r.data.fitgoodness = [p.data.fitgoodness; q.data.fitgoodness];
		r.data.fresult = {p.data.fresult{:} q.data.fresult{:}};
		r.data.xcwidth = [p.data.xcwidth; q.data.xcwidth];
		r.data.pciCounts = concatenate(p.data.pciCounts,q.data.pciCounts, ...
			'Columnwise');
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
