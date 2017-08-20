function r = plus(p,q,varargin)
%@trentropy/plus Overloaded plus function for TRENTROPY objects.
%   R = plus(P,Q) combines TRENTROPY objects P and Q and returns the
%   TRENTROPY object R.

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
	r.data.framebins = [r.data.framebins; q.data.framebins];
	r.data.entropy = [r.data.entropy; q.data.entropy];
	r.data.surrPercent = [r.data.surrPercent; q.data.surrPercent];
	r.data.surrZScores = [r.data.surrZScores; q.data.surrZScores];
	% get the last cell id in p
	plastid = p.data.cellid(end);
	% add plastid to q.data.cellid before adding 
	r.data.cellid = [r.data.cellid; q.data.cellid+plastid];
    r = set(r,'Number',r.data.cellid(end));
end
