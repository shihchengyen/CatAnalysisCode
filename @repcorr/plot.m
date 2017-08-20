function obj = plot(obj,varargin)

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b');
Args.flags = {'LabelsOff'};
[Args,varargin2] = getOptArgs(varargin,Args);

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
else
	% plot all data
end

if(~Args.LabelsOff)
	xlabel('X Axis')
	ylabel('Y Axis')
end
