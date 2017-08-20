function [obj, varargout] = plot(obj,varargin)
%@variability/plot Plot function for VARIABILITY object.
%   OBJ = plot(OBJ,VARARGIN) plots the variance versus the mean of the 
%   spike counts for all frames. The frames that have Fano Factors 
%   lower than 95% of the surrogates are plotted in red. These are the 
%   optional input arguments:
%      Numeric - plots the data for cell N.
%      'Indices' - plots the data for selected frames.
%      'TitleString' - specifies the title of the plot, usually used
%                      with the 'Indices' option.
%      'SurrPercentThreshold' - specifies a surrogate threshold 
%                               (default: 0.95).
%      'SpikeCountThreshold' - specifies the spike count threshold to 
%                              use in selecting windows for plotting 
%                              (default: 1).
%
%   OBJ = plot(OBJ,N,VARARGIN) plots the variance versus the mean for
%   all frames in the Nth cell.
%
%   OBJ = plot(OBJ,'Fano',VARARGIN) plots the Fano Factors for all
%   frames in order. 
%
%   OBJ = plot(OBJ,'Entropy',VARARGIN) plots the Time-Reliability
%   entropy for all frames in order.
%
%   These are the optional input arguments (in addition to those 
%   listed above): 
%      'SurrPercent' plots the percentage of surrogates with higher 
%                    Fano Factors for each frame in cell N.
%
%      'SurrZScores' plots the z-scores of the Fano Factor in each
%                    frame compared to the surrogate data.
%
%      'Hist'        plots a histogram of the Fano Factors instead.
%
%   OBJ = plot(OBJ,N,'Fano','ShowSurrogateData') Displays Fano 
%   Factors from all frames of cell N, along with the data from the
%   surrogates.
%
%   OBJ = plot(OBJ,N,'Entropy','ShowSurrogateData') Displays Time- 
%   Reliability entropy from all frames of cell N, along with the 
%   data from the surrogates.
%
%   OBJ = plot(OBJ,N,...,'XTime') plots the data for each frame as
%   a function of time.
%
%   OBJ = plot(OBJ,'Hist',VARARGIN) plots the distribution of Fano
%   Factors for windows with a mean spike count of at least 1, along
%   with the distribution of Fano Factors for windows with TR-Entropy
%   percentiles above 95%. Be sure to change the step size of the
%   histogram bins with the 'HistCStep' argument. 
%
%   OBJ = plot(OBJ,'FanoEntZScores',VARARGIN) plots the Fano Factors
%   versus the entropy Z-Scores for each window.
%
%   OBJ = plot(OBJ,'MeanEnt',VARARGIN) plots the mean spike count versus
%   the entropy for each window.
%
%   OBJ = plot(OBJ,'MeanFano',VARARGIN) plots the mean spike count versus
%   the Fano Factor for each window.
%
%   In addition to the optional input arguments specified above, the
%   following arguments are applicable to the last four types of plots:
%      'SurrZScoreThreshold' - specifies the z-score threshold (default:
%                              3).
%      'Hist' - flag indicating that a histgram of the distribution of
%               data should be plotted.
%      'HistCStep' - specifies the step size for the histogram (default:
%                    0.05).

Args = struct('SurrPercentThreshold',0.95, ...
				'SurrZScoreThreshold',3, ...
				'SpikeCountThreshold',1, ...
				'UsePercentThreshold',0, ...
				'UseZScoreThreshold',0, ...
				'Indices',[], ...
				'TitleString','', ...
				'Fano',0, ...
				'Entropy',0, ...
				'SurrPercent',0, ...
				'SurrZScores',0, ...
				'FanoEntZScores',0, ...
				'MeanEnt',0, ...
				'MeanFano',0, ...
				'XTime',0, ...
				'Hist',0,'HistXStep',0.05,'HistYStep',0.05, ...
				'HistXMax',[],'HistYMax',[], ...
				'HistXMin',[],'HistYMin',[], ...
				'ShowSurrogateData',0, ...
				'SurrogateFF','framesgFF.mat', ...
				'SurrogateTRE','framesgTRE10.mat', ...
                'NoCellSep',0, ...
				'LineStyle','none', ...
				'DataBMPlotColor',[0.5 0.5 0.5], ...
				'DataBMSPlotColor','b', ...
				'DataPlotColor','k', ...
				'DataSigPlotColor','r', ...
				'DataPlotSymbol','.', ...
				'DataPlotSize',6, ...
				'SurrogatePlotColor','g', ...
				'SurrogatePlotSymbol','.', ...
				'SurrogatePlotSize',6, ...
				'FontSize',10, ...
                'ReturnVars', {''}, ...
                'GroupPlots',1,'GroupPlotIndex',1, ...
                'ArgsOnly',0,'ShowSurrogateRasters',{''});
Args.flags = {'Fano','Entropy','SurrPercent','SurrZScores', ...
	'FanoEntZScores', 'MeanEnt','MeanFano','XTime','Hist', ...
	'ShowSurrogateData','UsePercentThreshold','UseZScoreThreshold','NoCellSep', ...
    'ArgsOnly'};
Args = getOptArgs(varargin,Args, ...
	'aliases',{'HistCStep' {'HistXStep','HistYStep'}; ...
				'HistMax' {'HistXMax','HistYMax'}; ...
				'HistMin' {'HistXMin','HistYMin'}});
            
% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    varargout{1} = {'Args',Args};
    return;
end

% define constants
NONE = 0;
HIST = 1;
DATA1 = 2;
DATA2 = 3;

% initialize flags
plotAllData = 0;
plotOneData = 0;
plotType = NONE;

plotCellSep = 0;

if(~isempty(Args.Indices))
	indices = Args.Indices;	
	titlestr = Args.TitleString;
elseif(~isempty(Args.NumericArguments))
	n = Args.NumericArguments{1};
	% get the indices that go with cell n
	indices = find(obj.data.cellid==n);
	% pass varargin to getDataDirs in case there is a prefix we want to 
	% remove
	titlestr = getDataDirs('ShortName','DirString',obj.data.cellname{n});
	plotOneData = 1;
else
	% default is to plot all data
	indices = 1:length(obj.data.scmean);
	titlestr = ['Data from ' num2str(length(obj.data.cellname)) ' cells'];
	plotAllData = 1;
end

% find indices that have mean spike counts above and below
% the spike count threshold
mscindices = obj.data.scmean(indices) >= Args.SpikeCountThreshold;
amindices = find(mscindices);
bmindices = find(mscindices==0);

% find significant indices according to specified significance criteria
if(Args.Entropy)
	if( Args.UseZScoreThreshold ...
		|| (Args.SurrZScores && ~Args.UsePercentThreshold) )
		sigindices = find(-obj.data.entropySurrZScores(indices) ...
			>= Args.SurrZScoreThreshold);
	else
		sigindices = find(obj.data.entropySurrPercent(indices) ...
			>= Args.SurrPercentThreshold);
	end
else
	if( Args.UseZScoreThreshold ...
		|| (Args.SurrZScores && ~Args.UsePercentThreshold) )
		sigindices = find(-obj.data.fanoSurrZScores(indices) ...
			>= Args.SurrZScoreThreshold);
	else
		sigindices = find(obj.data.fanoSurrPercent(indices) ...
			>= Args.SurrPercentThreshold);
	end
end

% separate indices into significant and non-significant and above and
% below spike count threshold
% get the non-significant indices
nsigindices = setdiff((1:length(indices))',sigindices);
% get object indices for significant points above spike count
% threshold
amsigindices = intersect(sigindices,amindices);
% get object indices for significant points below spike count
% threshold
bmsigindices = intersect(sigindices,bmindices);
% get object indices for non-significant points above spike 
% count threshold
amnsigindices = intersect(nsigindices,amindices);
% get object indices for non-significant points below spike 
% count threshold
bmnsigindices = intersect(nsigindices,bmindices);
% concatenate indices into matrix for histcie
indmat = concatenate(concatenate(concatenate( ...
								amnsigindices,bmnsigindices,'Columnwise'), ...
					amsigindices,'Columnwise'), ...
		bmsigindices,'Columnwise');
% set column variables so it will be easy to change
AMNSIG = 1;
BMNSIG = 2;
AMSIG = 3;
BMSIG = 4;

% initialize SData so it will be empty if we don't find any surrogate
% data
SData = [];

% check if we need to go to cluster directory, which is only
% relevant when we are plotting individual data set
if( (Args.XTime || Args.ShowSurrogateData || ~isempty(Args.ShowSurrogateRasters) ) && plotOneData )
	% get current directory
	cwd = pwd;
	% change to relevant directory
	cd(obj.nptdata.SessionDirs{n});
	if(Args.XTime)
		% get frame times in milliseconds from adjspikes object
		% instantiate adjspikes object
		adjs = adjspikes('auto',varargin{:});
		% get frame times for first repetition from adjspikes
		% use length of indices since that should correspond to 
		% number of frames
		xvals = adjs.data.adjFramePoints(1:length(indices));
		% adjust for the fact that framepoints does not always
		% start from 0
		if(xvals(1)~=0)
			xvals = xvals - xvals(1);
		end
		% shift xvals so that data point will be plotted in
		% center of frame window when used with plot
		xvals = xvals + (xvals(2)-xvals(1))/2;
		xlabelstr = 'Time (ms)';
	else
		% get frame numbers
		xvals = 1:length(indices);
		xlabelstr = 'Frame Number';
	end
	if(Args.ShowSurrogateData)
		% we might need fano surrogate data for the default option
		% which is to plot mean vs variance so just check to make sure
		% we are not plotting Entropy data
		if(~Args.Entropy)
			% check for presence of surrogate data
			if(ispresent(Args.SurrogateFF,'file'))
				fprintf('Loading surrogate data...\n');
				mat = load(Args.SurrogateFF);
				fn = fieldnames(mat);
				SData = getfield(mat,fn{1});
			end
        else
			% check for presence of surrogate data
			if(ispresent(Args.SurrogateTRE,'file'))
				fprintf('Loading surrogate data...\n');
				mat = load(Args.SurrogateTRE);
				fn = fieldnames(mat);
				SData = getfield(mat,fn{1});
			end
		end
    elseif(~isempty(Args.ShowSurrogateRasters))
        % open the raster file
        sb = readSurrogateBin(Args.ShowSurrogateRasters{1});
        SData = cell2array(sb{Args.ShowSurrogateRasters{2}});
	end
	% return to original directory
	cd(cwd);
else
	xvals = 1:length(indices);
	xlabelstr = 'Frame Number';
end

cla
if(Args.Fano)
	if(Args.SurrPercent)
		% plot percent of surrogates with lower Fano
		data1 = obj.data.fanoSurrPercent(indices);
		if(Args.Hist)
			plotType = HIST;
			histcedges = 0:Args.HistXStep:1;
			xlabelstr = 'Percent of surrogates with higher Fano Factor';
			ylabelstr = 'Number of Frames';
		else
			plotType = DATA1;					
			ylabelstr = 'Percent of surrogates with higher Fano Factor';
		end
	elseif(Args.SurrZScores)
		% plot Fano z-scores compared to surrogates
		data1 = obj.data.fanoSurrZScores(indices);
		if(Args.Hist)
			plotType = HIST;		
			histcedges = (floor(min(data1)/Args.HistXStep) ...
				* Args.HistXStep):Args.HistXStep:(ceil(max(data1) ...
				/Args.HistXStep)*Args.HistXStep);
			xlabelstr = 'Fano Z-Score';
			ylabelstr = 'Number of Frames';
		else
			plotType = DATA1;
			ylabelstr = 'Fano Z-Score';
		end
	else
		% plot Fano Factors
		plotType = DATA1;
		if(Args.ShowSurrogateData)
			if(~isempty(SData))
				% get Fano Factors from mean and std
				warning off MATLAB:divideByZero
				surrFanos = ((SData.scstd).^2)./SData.scmean;
				warning on MATLAB:divideByZero
				% get mean Fano
				sFmean = nanmean(surrFanos,2);
				% get Fano std
				sFstd = nanstd(surrFanos,0,2);
				% plot mean and std
				% use line instead of errorbar to remove ticks at
				% the end of the errorbars which overlap quite a bit
				cla
				line(repmat(xvals,2,1),[sFmean-sFstd sFmean+sFstd]', ...
					'Color',Args.SurrogatePlotColor);
				hold on
				plot(xvals,sFmean,'Color',Args.SurrogatePlotColor, ...
					'Marker',Args.SurrogatePlotSymbol, ...
					'MarkerSize',Args.SurrogatePlotSize, ...
					'LineStyle',Args.LineStyle);
			end
		end
		data1 = obj.data.fano(indices);
		ylabelstr = 'Fano Factors';
	end
elseif(Args.Entropy)
	if(Args.SurrPercent)
		% plot percent of surrogates that have higher entropy
		data1 = obj.data.entropySurrPercent(indices);
		if(Args.Hist)
			plotType = HIST;
			histcedges = 0:Args.HistXStep:1;
			xlabelstr = 'Percent of surrogates with higher Entropies';
			ylabelstr = 'Number of Frames';
		else
			plotType = DATA1;
			ylabelstr = 'Percent of surrogates with higher Entropies';
		end
	elseif(Args.SurrZScores)
		% plot z-scores of entropy compared to surrogate data
		data1 = obj.data.entropySurrZScores(indices);
		if(Args.Hist)
			plotType = HIST;
			histcedges = (floor(min(data1)/Args.HistXStep) ...
				* Args.HistXStep):Args.HistXStep:(ceil(max(data1) ...
				/ Args.HistXStep)*Args.HistXStep);
			xlabelstr = 'Entropy Z-Score';
			ylabelstr = 'Number of Frames';
		else
			plotType = DATA1;
			ylabelstr = 'Entropy Z-Score';
		end
	else
		% plot entropy
        plotType = DATA1;
		if(Args.ShowSurrogateData)
			if(~isempty(SData))
				% get mean entropy
				sTmean = nanmean(SData);
				% get entropy std
				sTstd = nanstd(SData);
				% plot mean and std
				% use line instead of errorbar to remove ticks at
				% the end of the errorbars which overlap quite a bit
				cla
				line(repmat(xvals,2,1),[sTmean-sTstd; sTmean+sTstd], ...
					'Color',Args.SurrogatePlotColor);
				hold on
				plot(xvals,sTmean,'Color',Args.SurrogatePlotColor, ...
					'Marker',Args.SurrogatePlotSymbol, ...
					'MarkerSize',Args.SurrogatePlotSize, ...
					'LineStyle',Args.LineStyle);
			end
		end
		data1 = obj.data.entropy(indices);
		ylabelstr = 'Entropy';
	end
elseif(Args.FanoEntZScores)
% 	% get windows with mean spike count >= 1
% 	msc1 = find(obj.data.scmean(indices)>=Args.SpikeCountThreshold);
% 	% get FF for those windows
% 	msc1ff = obj.data.fano(indices(msc1));
% 	% get entropy Z-scores for those windows
% 	msc1entz = obj.data.entropySurrZScores(indices(msc1));
% 	% construct data
% 	data = [msc1entz msc1ff];
    % get windows with mean spike count >=1 and fano factor <1
    mscff = find(obj.data.scmean(indices)>=Args.SpikeCountThreshold & obj.data.fano(indices)<1);
    % get fano Z-scores for those windows
    FFZScore = obj.data.fanoSurrZScores(indices(mscff));
    % get entropy for those windows
    EnZScore = obj.data.entropySurrZScores(indices(mscff));
    % construct data
    data = [EnZScore FFZScore];
	if(Args.Hist)
		% get 2D histogram
		[n,bins] = histn(data);
		% plot histogram
		imagesc(bins(:,1),bins(:,2),n')
		% flip the y axis
		set(gca,'YDir','normal')
	else
		plot(EnZScore,FFZScore,'Color',Args.DataPlotColor, ...
			'Marker',Args.DataPlotSymbol, ...
			'MarkerSize',Args.DataPlotSize, ...
			'LineStyle',Args.LineStyle);
		% get Pearson correlation coefficient
		[r1,p1] = corrcoef(data);
		% get Spearman correlation coefficient
		[r2,p2] = spcorr(data);
        if(length(r1)>1)
			% get axis limits so we can position text
			ax1 = axis;
			% put up coefficients
			text(ax1(1),0.95*ax1(4), ...
				['Pearson r = ' num2str(r1(1,2)) ' p = ' num2str(p1(1,2))])
			text(ax1(1),0.9*ax1(4), ...
				['Spearman r = ' num2str(r2(1,2)) ' p = ' num2str(p2(1,2))])
        end
	end
% 	ylabelstr = 'Fano Factors';
%   xlabelstr = 'TR-Entropy z-scores';
    ylabelstr = 'Fano Z-Score';
	xlabelstr = 'Entropy Z-Score';
elseif(Args.MeanEnt)
	plot(obj.data.scmean(indices),obj.data.entropy(indices), ...
		'Color',Args.DataPlotColor, ...
		'Marker',Args.DataPlotSymbol, ...
		'MarkerSize',Args.DataPlotSize, ...
		'LineStyle',Args.LineStyle)
	xlabelstr = 'Mean spike count';
	ylabelstr = 'TR-Entropy';
elseif(Args.MeanFano)
	if(Args.Hist)
		% get data
		scmean = obj.data.scmean(indices);
		fano = obj.data.fano(indices);
		
		% get hist steps for scmean
		scmstep = Args.HistXStep;
		if(isempty(Args.HistXMax))
			% no max specified so make sure the bins cover the last 
			% point in scmean
			scmmax = ceil(max(scmean)/scmstep) * scmstep;
		else
			% use specified max but make sure bins include specified max
			scmmax = ceil(Args.HistXMax/scmstep) * scmstep;
		end
		if(isempty(Args.HistXMin))
			% no min specified so just start from 0
			scmmin = 0;
		else
			% use specified min but make sure bins include specified min
			scmmin = floor(Args.HistXMin/scmstep) * scmstep;
		end
		% set up bins for histogram starting from 0
		scmbins = scmmin:scmstep:scmmax;
		
		% get hist steps for fano
		fanostep = Args.HistYStep;
		if(isempty(Args.HistYMax))
			% no max specified so make sure the bins cover the last 
			% point in fano
			fanomax = ceil(max(fano)/fanostep) * fanostep;
		else
			% use specified max but make sure bins include specified max
			fanomax = ceil(Args.HistYMax/fanostep) * fanostep;
		end
		if(isempty(Args.HistYMin))
			% no min specified so just start from 0 
			fanomin = 0;
		else
			% use specified min but make sure bins include specified min
			fanomin = floor(Args.HistYMin/fanostep) * fanostep;
		end
		% set up bins for histogram starting from 0
		fanobins = fanomin:fanostep:fanomax;
		% histn puts the 1st column in y so put [fano scmean]
		n1 = histn([fano scmean],concatenate(fanobins,scmbins)');
		imagesc(scmbins,fanobins,n1)
		% flip the y-axis
		set(gca,'YDir','normal');
		xlabelstr = 'Mean spike count';
		ylabelstr = 'Fano Factors';
	else
		plotType = DATA2;
		data1 = obj.data.scmean(indices);
		data2 = obj.data.fano(indices);
		xlabelstr = 'Mean spike count';
		ylabelstr = 'Fano Factors';
    end
elseif(~isempty(Args.ShowSurrogateRasters))
        sptsize = size(SData);
        yi = repmat(1:sptsize(2),sptsize(1),1);
        plot(SData,yi,'k.')
        xlabelstr = 'Time';
        ylabelstr = 'Repetitions';
else
	% plot mean spike count versus spike count variance
	if(Args.Hist)
		% get windows with mean spike count >= 1
		msc1 = find(obj.data.scmean(indices)>=Args.SpikeCountThreshold);
        if(~isempty(msc1))
			% get FF for those windows
			msc1ff = obj.data.fano(indices(msc1));
			% get max FF
			maxmsc1ff = ceil(max(msc1ff)/Args.HistXStep)*Args.HistXStep;
			% get bins for histogram
			ffedges = 0:Args.HistXStep:maxmsc1ff;
			% calculate histogram
			n = histcie(msc1ff,ffedges);
			bar(ffedges,n,'histc')
			hold on
			% find windows that are above the threshold for 
			% EntSurrThresh
			msc1tre95 = find(obj.data.entropySurrPercent(indices(msc1)) ...
				>=Args.SurrPercentThreshold);
            if(~isempty(msc1tre95))
				% get FF for those windows
				msc1tre95ff = obj.data.fano(indices(msc1(msc1tre95)));
     			% calculate the histogram
				n2 = histcie(msc1tre95ff,ffedges);
				% bar(ffedges,[n2 (n-n2)],'stacked','histc');
				h = bar(ffedges,n2,'histc');
                set(h,'FaceColor','r');
            end
        end
		xlabelstr = 'Fano Factors';
		ylabelstr = 'Number of Frames';
		hold off
	else
		if(Args.ShowSurrogateData)
			if(~isempty(SData))
				% convert data to single vector so the number of 
				% children in the figure is reduced 
				plot(SData.scmean(:),(SData.scstd(:)).^2, ...
					'Color',Args.SurrogatePlotColor, ...
					'Marker',Args.SurrogatePlotSymbol, ...
					'MarkerSize',Args.SurrogatePlotSize, ...
					'LineStyle',Args.LineStyle)
				hold on
			end
		end
		% plot the variance versus the mean
		% plot non-significant plots above spike count threshold
		plot(obj.data.scmean(indices(amnsigindices)),obj.data.scstd(indices(amnsigindices)).^2, ...
			'Color',Args.DataPlotColor, ...
			'Marker',Args.DataPlotSymbol, ...
			'MarkerSize',Args.DataPlotSize,'LineStyle',Args.LineStyle);
		hold on
		% plot non-significant plots below spike count threshold
		plot(obj.data.scmean(indices(bmnsigindices)),obj.data.scstd(indices(bmnsigindices)).^2, ...
			'Color',Args.DataBMPlotColor, ...
			'Marker',Args.DataPlotSymbol, ...
			'MarkerSize',Args.DataPlotSize,'LineStyle',Args.LineStyle);
        % plot significant plots above spike count threshold
		plot(obj.data.scmean(indices(amsigindices)),obj.data.scstd(indices(amsigindices)).^2, ...
			'Color',Args.DataSigPlotColor, ...
			'Marker',Args.DataPlotSymbol, ...
			'MarkerSize',Args.DataPlotSize, ...
			'LineStyle',Args.LineStyle)
		% plot significant plots below spike count threshold
		plot(obj.data.scmean(indices(bmsigindices)),obj.data.scstd(indices(bmsigindices)).^2, ...
			'Color',Args.DataBMSPlotColor, ...
			'Marker',Args.DataPlotSymbol, ...
			'MarkerSize',Args.DataPlotSize,'LineStyle',Args.LineStyle);
        
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
		xlabelstr = 'Mean spike count';
		ylabelstr = 'Variance';
	end
end

if(plotType==HIST)
	n = histcie(nanindex(data1,indmat),histcedges,'DataCols','DropLast');
	% plot all conditions
	% we really want to use the histc option but we can't use
	% both 'stacked' and 'histc' so we will use stacked so we
	% can plot all 4 data types with one command and then just
	% shift the lables to make it look like the 'histc' option
	% was used
	h = logbar(histcedges(1:(end-1)),n,'stacked');
	% change colors
	set(h(AMSIG),'FaceColor',Args.DataSigPlotColor);
	set(h(BMSIG),'FaceColor',Args.DataBMSPlotColor);
	set(h(AMNSIG),'FaceColor',Args.DataPlotColor);
	set(h(BMNSIG),'FaceColor',Args.DataBMPlotColor);
	% change tick marks so that it corresponds to histc
	axis tight
	xticks = [get(gca,'XTick') histcedges(end)];
	xticks2 = xticks - (histcedges(2)-histcedges(1))/2;
	set(gca,'XTick',xticks2,'XTickLabel',xticks);
	% plot points at labels 0 and 1 so that the final axis tight
	% ocmmand will not change the axis limits
	hold on
	% get min y value
	y = ylim;
	plot([xticks2(1) xticks2(end)],[y(1) y(1)],'Marker','.', ...
		'Color','k','LineStyle','none','MarkerSize',1);
	hold off
elseif(plotType==DATA1)
	plot(xvals(amnsigindices),data1(amnsigindices), ...
		'Color',Args.DataPlotColor, ...
		'Marker',Args.DataPlotSymbol, ...
		'MarkerSize',Args.DataPlotSize, ...
		'LineStyle',Args.LineStyle)
	hold on
	% plot non-significant plots below spike count threshold
	plot(xvals(bmnsigindices),data1(bmnsigindices), ...
		'Color',Args.DataBMPlotColor, ...
		'Marker',Args.DataPlotSymbol, ...
		'MarkerSize',Args.DataPlotSize, ...
		'LineStyle',Args.LineStyle)
	% plot significant plots above spike count threshold
	plot(xvals(amsigindices),data1(amsigindices), ...
		'Color',Args.DataSigPlotColor, ...
		'Marker',Args.DataPlotSymbol, ...
		'MarkerSize',Args.DataPlotSize, ...
		'LineStyle',Args.LineStyle)
	% plot significant plots below spike count threshold
	plot(xvals(bmsigindices),data1(bmsigindices), ...
		'Color',Args.DataBMSPlotColor, ...
		'Marker',Args.DataPlotSymbol, ...
		'MarkerSize',Args.DataPlotSize, ...
		'LineStyle',Args.LineStyle)
	if(plotAllData && ~Args.NoCellSep)
		% draw lines between cells
		% get number of frames for each cell and then add 0.5 so
		% we draw the lines between points instead of on top of points
		cellLimits = find(diff(obj.data.cellid)) + 0.5;
        if(~isempty(cellLimits))
    		plot(repmat(cellLimits',2,1),ylim,'k');
        end
	end
	hold off
elseif(plotType==DATA2)
	% plot non-significant plots above spike count threshold
	plot(data1(amnsigindices),data2(amnsigindices), ...
		'Color',Args.DataPlotColor, ...
		'Marker',Args.DataPlotSymbol, ...
		'MarkerSize',Args.DataPlotSize, ...
		'LineStyle',Args.LineStyle);
	hold on
	% plot non-significant plots below spike count threshold
	plot(data1(bmnsigindices),data2(bmnsigindices), ...
		'Color',Args.DataBMPlotColor, ...
		'Marker',Args.DataPlotSymbol, ...
		'MarkerSize',Args.DataPlotSize, ...
		'LineStyle',Args.LineStyle);
	% plot significant plots above spike count threshold
	plot(data1(amsigindices),data2(amsigindices), ...
		'Color',Args.DataSigPlotColor, ...
		'Marker',Args.DataPlotSymbol, ...
		'MarkerSize',Args.DataPlotSize, ...
		'LineStyle',Args.LineStyle)
	% plot significant plots below spike count threshold
	plot(data1(bmsigindices),data2(bmsigindices), ...
		'Color',Args.DataBMSPlotColor, ...
		'Marker',Args.DataPlotSymbol, ...
		'MarkerSize',Args.DataPlotSize, ...
		'LineStyle',Args.LineStyle)
	hold off
end

if(plotCellSep)
	hold on
	% draw lines between cells
	% get number of frames for each cell and then add 0.5 so
	% we draw the lines between points instead of on top of points
	cellLimits = find(diff(obj.data.cellid)) + 0.5;
	plot(repmat(cellLimits',2,1),ylim,'k');
	hold off
end

axis tight
set(gca,'FontSize',Args.FontSize);
xlabel(xlabelstr,'FontSize',Args.FontSize);
ylabel(ylabelstr,'FontSize',Args.FontSize);
title(titlestr,'FontSize',Args.FontSize);

rvarl = length(Args.ReturnVars);
if(rvarl>0)
     % assign requested variables to varargout
     for rvi = 1:rvarl
     	 rvidx = rvi * 2;
         varargout{1}{rvidx-1} = Args.ReturnVars{rvi};
         varargout{1}{rvidx} = eval(Args.ReturnVars{rvi});
     end
end
