function [r,varargout] = get(obj,varargin)
%shufflesync/get Get function for shufflesync objects

Args = struct('AnalysisLevel',0,'flags',{{'AnalysisLevel'}},'Number',0, ...
	'Windowed',0);
Args.flags = {'Number','Windowed'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.AnalysisLevel)
	r = 'Pairs';
elseif(Args.Number)
	if(Args.Windowed)
		r = size(obj.data.WinStartTime,1);
	end
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end