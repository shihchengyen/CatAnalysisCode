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
		r.data.frdiff = concatenate(p.data.frdiff,q.data.frdiff);
		r.data.setIndex = [p.data.setIndex; (p.data.setIndex(end) ...
			+ q.data.setIndex(2:end))];
		r.data.frdmean = [p.data.frdmean; q.data.frdmean];
		r.data.vmag = [p.data.vmag; q.data.vmag];
		r.data.vang = [p.data.vang; q.data.vang];
		r.data.vangp = [p.data.vangp; q.data.vangp];
		r.data.frdmax = [p.data.frdmax; q.data.frdmax];
        r.data.frdiffctrl = concatenate(p.data.frdiffctrl,q.data.frdiffctrl);
		r.data.frdctrlmean = [p.data.frdctrlmean; q.data.frdctrlmean];
		r.data.frdctrlmax = [p.data.frdctrlmax; q.data.frdctrlmax];
		r.data.corrcoef = [p.data.corrcoef; q.data.corrcoef];
		r.data.corrpvalue = [p.data.corrpvalue; q.data.corrpvalue];
		r.data.spcorrcoef = [p.data.spcorrcoef; q.data.spcorrcoef];
		r.data.spcorrpvalue = [p.data.spcorrpvalue; q.data.spcorrpvalue];
		r.data.noisecorr = [p.data.noisecorr; q.data.noisecorr];
		r.data.noisecorrpvalue = [p.data.noisecorrpvalue; q.data.noisecorrpvalue];
		r.data.winnoisecorr = concatenate(p.data.winnoisecorr,q.data.winnoisecorr, ...
                                'Columnwise');
		r.data.winnoisecorrp = concatenate(p.data.winnoisecorrp,q.data.winnoisecorrp, ...
                                'Columnwise');
		r.data.timebins = concatenate(p.data.timebins,q.data.timebins, ...
                                'Columnwise');
        r.data.frmat = concatenate(p.data.frmat,q.data.frmat,'Columnwise');
        r.data.frSetIndex = [p.data.frSetIndex; (p.data.frSetIndex(end) ...
			+ q.data.frSetIndex(2:end))];
		r.data.winstd = concatenate(p.data.winstd,q.data.winstd,'Columnwise');

		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
