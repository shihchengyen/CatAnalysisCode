function obj = plot(obj,varargin)

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b','VAng',0, ...
    'ShowCont',0,'iSpikes',0);
Args.flags = {'LabelsOff','VAng','ShowCont','iSpikes'};
[Args,varargin2] = getOptArgs(varargin,Args);

if(~isempty(Args.NumericArguments))
	n = Args.NumericArguments{1};
	% plot one data set at a time
    if(Args.VAng)
        % plot vector angle
		boxplot(obj.data.vang(:,n),1)
        if(Args.ShowCont)
            hold on
            % plot continuous correlation coefficient
            plot(obj.data.contvang(n),'m.')
            hold off
        end
		if(~Args.LabelsOff)
			ylabel('Vector Angle (deg)')
		end
    elseif(Args.iSpikes)
        % get cluster directories
        % get current directory
        cwd = pwd;
        cd(obj.data.setNames{n})
        clusterDirs = getDataDirs('GetClusterDirs');
       	cd(clusterDirs{1})
        % get ispikes
        isp1 = ispikes('auto');
        cd(clusterDirs{2})
        isp2 = ispikes('auto');
        cd(cwd)
        subplot(2,1,1)
        stem(isp1.data.trial.cluster.spikes)
        subplot(2,1,2)
        stem(isp2.data.trial.cluster.spikes)        
    else
        % plot correlation coefficient
		boxplot(obj.data.ascorr(:,n),1)
        if(Args.ShowCont)
            hold on
            % plot continuous correlation coefficient
            plot(obj.data.contcoefs(n),'m.')
            hold off
        end
		if(~Args.LabelsOff)
			ylabel('Correlation Coefficients')
		end
    end
else
	% plot all data
    if(Args.VAng)
        % plot vector angles
        boxplot(obj.data.vang,1);
        GroupedIntraPlot(obj,varargin2{:});
        if(Args.ShowCont)
            hold on
            plot(obj.data.contvang,'m.')
            hold off
        end
		if(~Args.LabelsOff)
			xlabel('Pairs')
            ylabel('Vector Angle (deg)')
		end
    else
        % plot correlation coefficients
		boxplot(obj.data.ascorr,1);
        GroupedIntraPlot(obj,varargin2{:});
        if(Args.ShowCont)
            hold on
            plot(obj.data.contcoefs,'m.')
            hold off
        end
		if(~Args.LabelsOff)
			xlabel('Pairs')
            ylabel('Correlation Coefficients')
		end
    end
end
