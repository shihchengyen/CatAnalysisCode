function r = plus(p,q,varargin)
%@variability/plus Overloaded plus function for VARIABILITY objects.
%   R = plus(P,Q) combines VARIABILITY objects P and Q and returns the
%   VARIABILITY object R.

% check for empty object
if(isempty(q.data.cellname))
	r = p;
elseif(isempty(p.data.cellname))
	r = q;
else
	% combine objects
	% assign p to r so that we can be sure we are returning an object
	r = p;
	r.data.cellname = {r.data.cellname{:} q.data.cellname{:}};
	r.data.scmean = [r.data.scmean; q.data.scmean];
	r.data.scstd = [r.data.scstd; q.data.scstd];
	r.data.fano = [r.data.fano; q.data.fano];
	% get the last cell id in p
	plastid = p.data.cellid(end);
	% add plastid to q.data.cellid before adding 
	r.data.cellid = [r.data.cellid; q.data.cellid+plastid];
	r.data.fanoSurrPercent = [r.data.fanoSurrPercent; ...
		q.data.fanoSurrPercent];
	r.data.fanoSurrZScores = [r.data.fanoSurrZScores; ...
		q.data.fanoSurrZScores];
	r.data.framebins = [r.data.framebins; q.data.framebins];
	r.data.entropy = [r.data.entropy; q.data.entropy];
	r.data.entropySurrPercent = [r.data.entropySurrPercent; ...
		q.data.entropySurrPercent];
	r.data.entropySurrZScores = [r.data.entropySurrZScores; ...
		q.data.entropySurrZScores];
	r.nptdata = plus(p.nptdata,q.nptdata);
end
