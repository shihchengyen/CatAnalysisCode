function obj = plot(obj,varargin)

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
    'VAng',0,'Events',0,'XShift',[],'PSparse',0,'MeanStd',0);
Args.flags = {'LabelsOff','VAng','Events','PSparse','MeanStd'};
[Args,varargin2] = getOptArgs(varargin,Args);

if(~isempty(Args.NumericArguments))
	n = Args.NumericArguments{1};
	% plot one data set at a time
    if(Args.MeanStd)
        plot(obj.data.mfr(:,(n-1)*2+1),obj.data.mfr(:,n*2),'.')
    else
		hist(obj.data.psparse(:,n))
	
		if(~Args.LabelsOff)
			xlabel('Population Sparseness')
			ylabel('Number')
		end
    end
else
	% plot all data
    if(Args.VAng)
        % plot vector angle
        md = obj.data.medianvang;
        lq = md-obj.data.lqvang;
        uq = obj.data.uqvang-md;
        ylabelstr = 'Vector Angle';
    elseif(Args.Events)
        md = obj.data.medianjep;
        lq = md-obj.data.lqjep;
        uq = obj.data.uqjep-md;
        ylabelstr = 'Joint Event Probablity';
    elseif(Args.PSparse)
        md = obj.data.medianjep;
        lq = md-obj.data.lqjep;
        uq = obj.data.uqjep-md;
        ylabelstr = 'Population Sparseness';
    else
        % plot correlation
        md = obj.data.mediancorr;
        lq = md-obj.data.lqcorr;
        uq = obj.data.uqcorr-md;
        ylabelstr = 'PSTH Correlation';
    end
    xvals = 1:obj.data.numSets;
    if(~isempty(Args.XShift))
        xvals = xvals + Args.XShift;
    end
	errorbar(xvals,md,lq,uq,'rx');
	
	if(~Args.LabelsOff)
        ylabel(ylabelstr);
		xlabel('Recording Site');
    end
end

