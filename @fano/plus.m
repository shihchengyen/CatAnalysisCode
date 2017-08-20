function r = plus(p,q,varargin)
%@fano/plus Overloaded plus function for FANO objects.
%   R = plus(P,Q) combines FANO objects P and Q and returns the
%   FANO object R.

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
	r.data.surrPercent = [r.data.surrPercent; q.data.surrPercent];
	r.data.surrZScores = [r.data.surrZScores; q.data.surrZScores];
    r = set(r,'Number',r.data.cellid(end));
end
