function obj = plot(obj,varargin)

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
    'BoxPlotColors','bgry','BoxPlotSymbol','r+','BoxPlotMarkerSize',8);
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
    subplot('Position',[0.13 0.21 0.775 0.715])
	h = boxplot(obj.data.psparse,'notch','on','colors',Args.BoxPlotColors, ...
        'symbol',Args.BoxPlotSymbol);
    set(h(end,~isnan(h(end,:))),'MarkerSize',Args.BoxPlotMarkerSize)
    xl = xlim;
    set(gca,'XTickLabel',[])
    xlabel('')
	if(~Args.LabelsOff)
        ylabel('Population Sparseness')
    end
    subplot('Position',[0.13 0.11 0.775 0.1])
    bar(obj.data.ncells)
    yl = ylim;
    axis([xl 1 yl(2)])

	if(~Args.LabelsOff)
		xlabel('Sites')
        ylabel('Cells')
	end
end

