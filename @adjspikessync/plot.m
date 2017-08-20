function obj = plot(obj,varargin)
%@adjspikessync/plot Plot function for adjspikessync objects.

Args = struct('SurrHist',0,'ZScoreHist',0,'ZScoreThresh',3,'HistStep',0.5, ...
	'LabelsOff',0,'ZScore',0,'SurrMeanStd',0);
Args.flags = {'SurrHist','ZScoreHist','LabelsOff','ZScore','SurrMeanStd'};
Args = getOptArgs(varargin,Args);

if(isempty(Args.NumericArguments))
	% use the first data set
	n = 1;
else
	n = Args.NumericArguments{1};
end
if(Args.SurrHist)
	% get indices for sxc corresponding to pair n
	if(isempty(obj.data.WindowIndex))
		% use all of sxc
		sxc = obj.data.SurrSyncHist;
		startind = 1;
		endind = size(sxc,2);
	else
		startind = obj.data.WindowIndex(n) + 1;
		endind = obj.data.WindowIndex(n+1);
		sxc = obj.data.SurrSyncHist(:,startind:endind);
	end
	% get the indices corresponding to non-nans, which is the actual
	% range of synchoronous spikes for this pair
	sxci = ~isnan(sxc(:,1));
	% sxc = sxc(~isnan(sxc));
	sxc2 = sxc(sxci,:);
	lags = obj.data.WinStartTime(:,n);
	% find the non-nan values in lags which corresponds to
	% the number of windows for this pair (endind-startind+1)
	lagi = 1:(endind-startind+1);
	lags2 = lags(lagi);
	xc = obj.data.DataSyncSpikes(:,n);
	% find the non-nan values in xc which are the same as lagi
	xc2 = xc(lagi);
	xtime = lags2+(mean(diff(lags2))/2);
	% xtime = xtime(~isnan(xtime));
	imagesc(xtime,0:size(sxc2,1),sxc(sxci,:));
	hold on
	% find points with z-score larger than the threshold
	zi = find(abs(obj.data.DataZScore(:,n))>Args.ZScoreThresh);
	plot(xtime,xc2,'k*')
	plot(xtime(zi),xc2(zi),'r*')
	hold off
	if(~Args.LabelsOff)
		xlabel('Time (ms)')
		ylabel('Number of synchronous spikes')
	end
elseif(Args.ZScoreHist)
	data = obj.data.DataZScore(:,n);
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
	lags = obj.data.WinStartTime(:,n);
	% append value to lags to indicate duration of last bin
	xtime = [lags; lags(end)+mean(diff(lags))];
	% plot zscores using stairs
	data = [obj.data.DataZScore(:,n); 0];
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
	lags = obj.data.WinStartTime(:,n);
	% append value to lags to indicate duration of last bin
	xtime = [lags; lags(end)+mean(diff(lags))];
	% get the mean of the surrogates
	data1 = [obj.data.SurrSyncMean(:,n); 0];
	% get the std of the surrogates
	data2 = [obj.data.SurrSyncStd(:,n); 0];
	[a,h1,h2] = plotyy(xtime,data1,xtime,data2,'stairs');
	if(~Args.LabelsOff)
		xlabel('Time (ms)')
		set(get(a(1),'Ylabel'),'String','Mean Surrogate Spikes');
		set(get(a(2),'Ylabel'),'String','Std Surrogate Spikes');
	end
else
	% pad xc with a 0 so that the last data point will be drawn for the 
	% duration
	stairs(obj.data.WinStartTime(:,n),[obj.data.DataSyncSpikes(:,n); 0])
end
title(getDataDirs('ShortName','DirString',obj.data.setNames{n}))

zoom xon
