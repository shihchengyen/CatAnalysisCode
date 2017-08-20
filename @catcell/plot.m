function obj = plot(obj,varargin)
%@catcell/plot Plot function for CATCELL object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.
%
%   OBJ = plot(OBJ,'VarMean') plots the variance versus the mean
%   of the spike counts for each frame. If the file specified in
%   the optional input argument 'SurrogateFF' (default: 'framesgFF.mat')
%   is present, the data from the surrogates are loaded and plotted in 
%   green first.
%
%   OBJ = plot(OBJ,'Mean') plots the mean spike counts for each frame,
%   as well as the surrogates if they are present. 
%
%   OBJ = plot(OBJ,'Variance') plots the spike count variance for each 
%   frame, as well as the surrogates if they are present. 
%
%   OBJ = plot(OBJ,'Fano') plots the spike count Fano Factor for each 
%   frame, as well as the surrogates if they are present. 
%
%   OBJ = plot(OBJ,N,'Surrogates') plots the spike trains of the Nth
%   set of surrogate data. The surrogates are assumed to be contained
%   in files with basenames specified by the 'SurrogateBasename'
%   optional input argument (default: 'framesg'). The number of sets of 
%   surrogates in each file is specified by 'SurrogateSetsPerFile' 
%   (default: 100).
%
%   OBJ = plot(OBJ,N,'Latency') plots all spikes relative to the onset
%   of a stimulus frame. The bin size (in ms) can be specified using the 
%   'SubBinSize' optional input argument (default: 1).
%
%   OBJ = plot(OBJ,N,'EventIntervals') plots all event intervals. The 
%   name of the file containing the event intervals of the surrogates 
%   can be specified using the 'EventSurrogate' optional input argument
%   (default: 'framesgEvents3.mat'). If the 'Hist' optional input
%   argument is present, a subplot with the histogram of the data 
%   plotted above the histogram of the surrogates is plotted. The 
%   x-limits of this subplot can be specified using the 'EventXVals'
%   optional input argument.
%
%   OBJ = plot(OBJ,N,'MTpsd') plots amplitude of frequency spectrum
%   computed using multi-tapers. These are the optional input arguments:
%      MTBinSize - specifies the bin size used to compute the histogram
%                  of the spike times.
%      MTNW - specifies the NW paramter for the multi-taper analysis
%             (default: 4).
%      RepDuration - specifies the duration of the time series in ms
%                    (default: 30000).
%      MTWindow - specifies the length of each segment of data in data
%                 points (default: 4000).
%      MTWinStep - specifies the number of data points to take for the
%                  next window (default: 2000).
%      MTFreqLimits - specifies the range of frequencies to display
%                     (default: [0 200]).
%
%   OBJ = plot(OBJ,N,'MTspec') plots the time-frequency power spectrum
%   computed using multi-tapers. The optional input arguments are the
%   same as for 'MTpsd'.
%
%   Example usage:
%   nd = ProcessDays(nptdata,'SessionCells','nptCellCmd', ...
%      'robj = nptdata(0,0,''Eval'',''ispresent(''''catcell.mat'''', ...
%      ''''file'''')'');');
%   InspectGUI(nd,'Objects',{'catcell',{}})
%   InspectGUI(nd,'Objects',{'catcell',{'VarMean'}})
%   InspectGUI(nd,'Objects',{'catcell',{'Mean'};'catcell', ...
%      {'Variance'};'catcell',{'Fano'}},'SP',[3 1])
%   cd 'a4/a402/a402g5c1';load catcell; InspectGUI(cc,'Surrogates')
%
%   obj = plot(obj);
%   obj = plot(obj,'VarMean','Mean','Variance','Fano',...
%       'SurrogateFF','framesgFF.mat');
%   obj = plot(obj,n,'Surrogates','SurrogateBasename','framesg', ...
%             'SurrogateSetsPerFile',100);

Args = struct('VarMean',0,'SurrogateFF','framesgFF.mat', ...
	'Surrogates',0,'Mean',0,'Variance',0,'Fano',0, ...
	'SurrogateBasename','framesg','SurrogateSetsPerFile',100, ...
	'Latency',0,'SubBinSize',1,'EventIntervals',0,'EventSurrogate',...
	'framesgEvents3.mat','EventXVals',0:5:50,'MTpsd',0,'MTBinSize',1,...
	'MTNW',4,'RepDuration',30000,'MTWindow',4000, ...
	'MTWinStep',2000,'MTFreqLimits',[0 200],'Hist',0,'MTspec',0);

Args = getOptArgs(varargin,Args,'flags', ...
		{'VarMean','Surrogates','Mean','Variance','Fano','Latency', ...
		'EventIntervals','MTpsd','Hist','MTspec'});

if(Args.VarMean)
	% check for presence of surrogate data spike counts and std
	if(ispresent(Args.SurrogateFF,'file'))
		load(Args.SurrogateFF)
		plot(sgSC.scmean,(sgSC.scstd).^2,'g.')
		hold on
	end
	% plot the variance versus the mean
	plot(obj.data.scmean,obj.data.scstd.^2,'Color',[0.5 0.5 0.5], ...
		'Marker','.','LineStyle','none');
	hold on
	% get axis limits
	ax1 = axis;
	% find the smaller of xmax or ymax
	lmax = min([ax1(2) ax1(4)]);
	% draw y = x line
	line([0 lmax],[0 lmax],'Color','k');
	% draw min variance scallops
	x = 0:0.1:1;
	minvar = x .* (1 - x);
	for i = 0:(ax1(2)-1)
		plot(i+x,minvar,'k')
	end
	hold off
	xlabel('Mean spike count')
	ylabel('Variance')
elseif(Args.Mean)
	% check for presence of surrogate data spike counts and std
	if(ispresent(Args.SurrogateFF,'file'))
		load(Args.SurrogateFF)
		plot(sgSC.scmean,'g.')
		hold on
	end
	% plot the mean
	plot(obj.data.scmean,'r.');
	hold off
	ylabel('Mean spike count')
	xlabel('Frame number')
elseif(Args.Variance)
	% check for presence of surrogate data spike counts and std
	if(ispresent(Args.SurrogateFF,'file'))
		load(Args.SurrogateFF)
		plot((sgSC.scstd).^2,'g.')
		hold on
	end
	% plot the mean
	plot((obj.data.scstd).^2,'r.');
	hold off
	ylabel('Spike count variance')
	xlabel('Frame number')
elseif(Args.Fano)
	warning off MATLAB:divideByZero
	% check for presence of surrogate data spike counts and std
	if(ispresent(Args.SurrogateFF,'file'))
		load(Args.SurrogateFF)
		plot(((sgSC.scstd).^2)./sgSC.scmean,'g.')
		hold on
	end
	% plot the mean
	plot(((obj.data.scstd).^2)./obj.data.scmean,'r.');
	hold off
	ylabel('Spike count Fano Factor')
	xlabel('Frame number')
	warning on MATLAB:divideByZero
elseif(Args.Surrogates)
	% get number of plot
	n = Args.NumericArguments{1} - 1;
	% figure which file to load
	filen = floor(n/Args.SurrogateSetsPerFile) + 1;
	% figure which set in the file to plot
	setn = rem(n,Args.SurrogateSetsPerFile) + 1;
	% load surrogates
	sptrain = readSurrogateBin([Args.SurrogateBasename num2str(filen) '.bin']);
    cla
	plotRasters(sptrain{setn})
elseif(Args.Latency)
	% check if the data is already computed
	if(obj.data.latencySubBinSize~=Args.SubBinSize)
		% requested bin size is different from that already computed so
		% recompute
		% divide frame vector into 1 ms bins
		[subbins,sbinsize,nsbins] = divideBins( ...
			obj.data.stimulus_info.adjusted_frames_vector, ...
			'SubBinSize',1);
		% take histogram using 1 ms bins
		binsc = histcie(obj.data.cell_info.adjusted_spiketrain, ...
			subbins,'DropLast');
		% reshape into nsbins x (frames*reps)
		binsc1 = reshape(binsc,nsbins,[]);
		% take the sum of each bin and append 0 at the end so we can use
		% 'histc' option with bar
		binsums = [sum(binsc1,2); 0];
	else
		% data already computed so just extract from data structure
		binsums = obj.data.latencyBinSums;
		sbinsize = obj.data.latencyActualSubBinSize;
		nsbins = length(binsums);
	end
	% get xvals
    finalxval = (nsbins-1) * sbinsize;
	xvals = 0:sbinsize:finalxval;
	% plot distribution
	bar(xvals,binsums,'histc')
    hold on
    % replicate distribution to illustrate response latencies beyond the
    % frame limit
	bar(finalxval+xvals,binsums,'histc','c');
	% load data
	load(obj.data.cellname)
	% plot frequencies related to frequencies
	refresh = data.stimulus_info.video_refresh;
	refreshDuration = 1000/refresh;
	% get axis limits
	ax1 = axis;
	if(refresh>130)
		% 150 Hz so there are 6 refreshes
		line(repmat(refreshDuration * [1:6],2,1),repmat([ax1(3);ax1(4)],1,6), ...
			'Color','r')
	elseif(refresh>90)
		% 120 Hz so there are 4 refreshs
		line(repmat(refreshDuration * [1:4],2,1),repmat([ax1(3);ax1(4)],1,4), ...
			'Color','r')
	else
		% 85 Hz so there are 3 refreshs
		line(repmat(refreshDuration * [1:3],2,1),repmat([ax1(3);ax1(4)],1,3), ...
			'Color','r')
	end
    hold off
elseif(Args.EventIntervals)
	% load event intervals for surrogates
	if(ispresent(Args.EventSurrogate,'file'))
		m = load(Args.EventSurrogate);
		if(size(m.event.eventIntervals,2)>1)
			% the event intervals were computed for some surrogates with
			% a 1000x1000 cell array instead of 1000x1 so we have to 
			% extract the first 1000 elements
			eventIntervals = cell(1000,1);
			for i = 1:1000
				eventIntervals{i} = m.event.eventIntervals{i};
			end
			eIntervals = cell2array(eventIntervals);
		else
			eIntervals = cell2array(m.event.eventIntervals);
		end
		% find the smallest interval for the 95th percentile of the 
		% surrogates
		% sort intervals for each surrogate
		eIs = sort(eIntervals);
		% first row is the smallest interval for each surrogate so grab
		% the row and sort it to get the distribution of smallest
		% intervals for all surrogates
		mineIs = sort(eIs(1,:)');
		% grab the 95th percentile, which would be 0.05 * 1000 = 50th
		% point
		surr95 = mineIs(50);
		if(Args.Hist)
			subplot(2,1,2)
			if(~isempty(eIntervals))
				n2 = histcie(eIntervals(:),Args.EventXVals);
				bar(Args.EventXVals,n2,'histc');
				ax1 = axis;
				line([surr95 surr95],[ax1(3) ax1(4)],'Color','r')
				title(sprintf('Smallest interval for the 95th percentile of the surrogates: %f\n',surr95))
			end
			subplot(2,1,1)
		end
	end
	% get event intervals for data
	load(obj.data.cellname)
	cla
    if(~isempty(data.cell_info.events.event_interval_durations))
    	if(Args.Hist)
			n = histcie(data.cell_info.events.event_interval_durations, ...
				Args.EventXVals);
			bar(Args.EventXVals,n,'histc');
		else
			eventIntervals = data.cell_info.events.event_interval_durations;
			line(repmat(eventIntervals,2,1),repmat([0;1],1, ...
			    length(eventIntervals)),'Color','b');
            if(ispresent(Args.EventSurrogate,'file'))			
                % plot 95th percentile of surrogate intervals
			    line([surr95 surr95],[0 1],'Color','r')
            end
			% plot frame interval
			frameInt = 1000/data.stimulus_info.video_refresh ...
                * data.stimulus_info.refreshes_per_frame;
			line([frameInt frameInt],[0 1],'Color','m','LineStyle',':')
			xlim([0 50])
		end
    end
elseif(Args.MTpsd)
	% load data
	load(obj.data.cellname)
	% compute and plot multi-taper power spectrum density
	[pxx,freq] = mtpsd(histcie(data.cell_info.original_spiketrain, ...
		0:Args.MTBinSize:Args.RepDuration,'DropLast'),Args.MTNW, ...
		(2*Args.MTNW)-1,2.^(ceil(log2(Args.MTWindow))), ...
		1000/Args.MTBinSize,Args.MTWindow,Args.MTWinStep);
	psdplot(pxx,freq,'Hz','db');
	% plot frequencies related to frequencies
	refresh = data.stimulus_info.video_refresh;
	% get axis limits
	ax1 = axis;
	if(refresh>130)
		% 150 Hz so there are 6 refreshes
		line(repmat(refresh/6 * [1:6],2,1),repmat([ax1(3);ax1(4)],1,6), ...
			'Color','r')
	elseif(refresh>90)
		% 120 Hz so there are 4 refreshs
		line(repmat(refresh/4 * [1:4],2,1),repmat([ax1(3);ax1(4)],1,4), ...
			'Color','r')
	else
		% 85 Hz so there are 3 refreshs
		line(repmat(refresh/3 * [1:3],2,1),repmat([ax1(3);ax1(4)],1,3), ...
			'Color','r')
	end
	xlim(Args.MTFreqLimits)
elseif(Args.MTspec)
	% load data
	load(obj.data.cellname)
	% compute and plot multi-taper power spectrum density
	[pxx,freq,tvec] = mtspec(histcie(data.cell_info.original_spiketrain, ...
		0:Args.MTBinSize:Args.RepDuration,'DropLast'),Args.MTNW, ...
		(2*Args.MTNW)-1,2.^(ceil(log2(Args.MTWindow))), ...
		Args.MTWindow,Args.MTWinStep,1000/Args.MTBinSize);
	imagesc(tvec, freq, pxx)
	xlabel('Time (seconds)');
	ylabel('Frequency (Hz)');
	% plot frequencies related to frequencies
	refresh = data.stimulus_info.video_refresh;
	% get axis limits
	ax1 = axis;
	if(refresh>130)
		% 150 Hz so there are 6 refreshes
		line(repmat([ax1(1);ax1(2)],1,6),repmat(refresh/6 * [1:6],2,1), ...
			'Color','r')
	elseif(refresh>90)
		% 120 Hz so there are 4 refreshs
		line(repmat([ax1(1);ax1(2)],1,4),repmat(refresh/4 * [1:4],2,1), ...
			'Color','r')
	else
		% 85 Hz so there are 3 refreshs
		line(repmat([ax1(1);ax1(2)],1,3),repmat(refresh/3 * [1:3],2,1), ...
			'Color','r')
	end
	ylim(Args.MTFreqLimits)
else
	% plot the rasters
	cla
	plotRasters(obj.data.cell_info.raster)
end

title(obj.data.cellname)
