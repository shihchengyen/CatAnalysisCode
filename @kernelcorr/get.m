function [r,varargout] = get(obj,varargin)
%kernelcorr/get Get function for frdiff objects

Args = struct('AnalysisLevel',0);
Args.flags = {'AnalysisLevel'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.AnalysisLevel)
	r = 'Pairs';
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end
