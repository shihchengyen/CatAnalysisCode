function obj = plot(obj,varargin)
%@adjspikesxcorr/plot Plot function for adjspikesxcorr objects.

Args = struct('XLim',[],'NumBins',10,'Fit',0,'Sliding',0,'SurrHist',0, ...
	'ZScoreHist',0,'ZScoreThresh',3,'HistStep',0.5,'LabelsOff',0, ...
    'ZScore',0,'SurrMeanStd',0,'NoSig',0);
Args.flags = {'Fit','Sliding','SurrHist','ZScoreHist','LabelsOff','ZScore', ...
	'SurrMeanStd','NoSig'};
Args = getOptArgs(varargin,Args);

if(~Args.Sliding)
    lagslength = length(obj.data.lags);
	if(~isempty(Args.XLim))
		% make sure XLim has 2 entries
		if(length(Args.XLim)~=2)
			error('Error: XLim should have 2 numbers!');
		else
			% find first value in obj.data.f that exceeds XLim(1)
			fimin = find(obj.data.lags<Args.XLim(1));
			% find first value in obj.data.f that exceed XLim(2)
			fimax = find(obj.data.lags>Args.XLim(2));
			if(isempty(fimax) || isempty(fimin))
				% no points found so just use all points
				plotpts = 1:lagslength;
			else
				plotpts = (fimin(1)+1):(fimax(1)-1);
			end
		end
	else
		% plot using lags
		plotpts = 1:lagslength;
	end
	x = obj.data.lags(plotpts);

	if(isempty(Args.NumericArguments))
		% plots population ecf
		hist(obj.data.ecf,Args.NumBins);	
	else
		% plot data from individual cell
		n = Args.NumericArguments{1};
		% get bin size
		binSize = obj.data.binSize;
		x2 = [x(1)-binSize x x(end)+binSize];
		if(Args.Fit)
			% get xc
			xc = obj.data.xc(plotpts,n);
			sxc = obj.data.sxc(plotpts,n);
			dxc = xc - sxc;
			% plot difference between xcorr and shift-predictor
			bar(x,dxc,'histc');
			hold on
			% find bins where sxc was lower than pciCounts
			% sxclower = find(sxc < obj.data.pciCounts(:,n));
			% if(~isempty(sxclower))
				% plot significant bars with different color
				% bar(x(sxclower),dxc(sxclower),'histc','r');
				% end
			cfobj = obj.data.fresult{n};
			% plot the Gaussian fit if there is one
			if(~isempty(cfobj))
				plot(cfobj);
				% turn off legend put up by cfit/plot
				legend off
			end
			xlim(Args.XLim)
			hold off
		else
			% plot xcorr
			bar(x,obj.data.xc(plotpts,n),'histc');
			hold on
			% plot shift-predictor
			% add points at the start and the end to complete the plot
			stairs(x2,[0; obj.data.sxc(plotpts,n); 0],'r');
            if(~Args.NoSig)
				% plot significance level
				stairs(x2,[0; obj.data.pciCounts(:,n); 0],'m');
            end
			hold off
		end
		% get legend string
		lstring{1} = obj.data.legendstr{n};
		% if(~isempty(lstring))
		%     legend(lstring);
		% end
		lstring{2} = ['ECF: ' num2str(obj.data.ecf(n))];
		lstring{3} = ['KSSTAT: ' num2str(obj.data.ksstat(n))];
		ax1 = axis;
		text(ax1(2)*0.75,ax1(4)*0.9,lstring,'HorizontalAlignment','Center');
		title(getDataDirs('ShortName','DirString',obj.data.setNames{n}));
		xlabel('Time (ms)')
		ylabel('Counts')
	end
else
	if(isempty(Args.NumericArguments))
		% use the first data set
		n = 1;
	else
		n = Args.NumericArguments{1};
	end
	if(Args.SurrHist)
		% get indices for sxc corresponding to pair n
        if(isempty(obj.data.pciCounts))
            % use all of sxc
            sxc = obj.data.sxc;
            startind = 1;
            endind = size(sxc,2);
        else
            startind = obj.data.pciCounts(n) + 1;
            endind = obj.data.pciCounts(n+1);
			sxc = obj.data.sxc(:,startind:endind);
        end
        % get the indices corresponding to non-nans, which is the actual
        % range of synchoronous spikes for this pair
        sxci = ~isnan(sxc(:,1));
        % sxc = sxc(~isnan(sxc));
        sxc2 = sxc(sxci,:);
		lags = obj.data.lags(:,n);
        % find the non-nan values in lags which corresponds to
        % the number of windows for this pair (endind-startind+1)
        lagi = 1:(endind-startind+1);
        lags2 = lags(lagi);
		xc = obj.data.xc(:,n);
        % find the non-nan values in xc which are the same as lagi
        xc2 = xc(lagi);
		xtime = lags2+(mean(diff(lags2))/2);
        % xtime = xtime(~isnan(xtime));
		imagesc(xtime,0:size(sxc2,1),sxc(sxci,:));
		hold on
		% find points with z-score larger than the threshold
		zi = find(abs(obj.data.ecf(:,n))>Args.ZScoreThresh);
		plot(xtime,xc2,'k*')
		plot(xtime(zi),xc2(zi),'r*')
		hold off
        if(~Args.LabelsOff)
			xlabel('Time (ms)')
			ylabel('Number of synchronous spikes')
        end
	elseif(Args.ZScoreHist)
		data = obj.data.ecf(:,n);
        minbin = floor(min(data)/Args.HistStep) * Args.HistStep;
        maxbin = ceil(max(data)/Args.HistStep) * Args.HistStep;
		bins = minbin:Args.HistStep:maxbin;
		counts = histcie(data,bins);
		bar(bins,counts,'histc')
        hold on
        % find the indices that exceed the zscore threshold
        zi = find(bins>=Args.ZScoreThresh);
        if(~isempty(zi))
            h = bar(bins(zi),counts(zi),'histc');
            set(h,'FaceColor','r');
        end
        hold off
        xlim([minbin maxbin])
        if(~Args.LabelsOff)
			xlabel('Z-Score')
			ylabel('Number of occurrences')
        end
    elseif(Args.ZScore)
		lags = obj.data.lags(:,n);
        % append value to lags to indicate duration of last bin
		xtime = [lags; lags(end)+mean(diff(lags))];
        % plot zscores using stairs
        data = [obj.data.ecf(:,n); 0];
        stairs(xtime,data);
        hold on
        % get axis limits
        alim = xlim';
        plot(alim,[Args.ZScoreThresh; Args.ZScoreThresh],'r');
        plot(alim,[-Args.ZScoreThresh; -Args.ZScoreThresh],'r');
        hold off
        if(~Args.LabelsOff)
            xlabel('Time (ms)')
            ylabel('Z-Score')
        end
    elseif(Args.SurrMeanStd)
    	lags = obj.data.lags(:,n);
        % append value to lags to indicate duration of last bin
		xtime = [lags; lags(end)+mean(diff(lags))];
		% get the mean of the surrogates
		data1 = [obj.data.kspvalue(:,n); 0];
		% get the std of the surrogates
		data2 = [obj.data.ksstat(:,n); 0];
		[a,h1,h2] = plotyy(xtime,data1,xtime,data2,'stairs');
		if(~Args.LabelsOff)
			xlabel('Time (ms)')
			set(get(a(1),'Ylabel'),'String','Mean Surrogate Spikes');
			set(get(a(2),'Ylabel'),'String','Std Surrogate Spikes');
		end
	else
		% pad xc with a 0 so that the last data point will be drawn for the 
		% duration
		stairs(obj.data.lags(:,n),[obj.data.xc(:,n); 0])
	end
	title(getDataDirs('ShortName','DirString',obj.data.setNames{n}))
end

zoom xon
