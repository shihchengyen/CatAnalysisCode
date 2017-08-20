function obj = plot(obj,varargin)

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b');
Args.flags = {'LabelsOff'};
[Args,varargin2] = getOptArgs(varargin,Args);

if(~isempty(Args.NumericArguments))
	n = Args.NumericArguments{1};
	% plot one data set at a time
	hist(obj.data.psparse(:,n))

	if(~Args.LabelsOff)
		xlabel('Population Sparseness')
		ylabel('Number')
	end
else
	% plot all data
	gind = groupDirs(obj);
	gdata = nanindex(obj.data.jtprob,gind);
	bh = bar(gdata','Grouped');
	
	if(~Args.LabelsOff)
        ylabel('Joint Event Probability')
		xlabel('Recording Site')
    end
end

