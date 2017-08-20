function [r,varargout] = get(obj,varargin)
%adjspikesxcorr/get Get function for adjspikesxcorr objects

Args = struct('AnalysisLevel',0,'GroupPlotProperties',0);
Args.flags = {'AnalysisLevel'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.AnalysisLevel)
	r = 'Pairs';
elseif(Args.GroupPlotProperties>0)
	r.separate = 'Horizontal';
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end
