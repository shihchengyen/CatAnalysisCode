function [obj,varargout] = plot(obj,varargin)
%@frdiff/plot Plot function for frdiff objects

Args = struct('HideTitle',0,'Color','b','GroupPlots',1,'GroupPlotIndex',1, ...
    'ArgsOnly',0,'ReturnVars',{''}, ...
    'Hist',0,'Max',0,'VMag',0,'VAng',0,'Corr',0,'HistBins',[], ...
	'YYMeanVMag',0,'YYVAngMax',0,'YScale','linear','XScale','linear', ...
	'NumHistBins',10,'HistBinSize',[],'YYMeanMax',0,'YYVAngVMag',0, ...
    'GroupedIntra',0,'CorrPThresh',0.05,'DisplaceBarSig',0.02, ...
    'VAngCorr',0,'NoiseCorr',0,'WinNoiseCorr',0,'FRMean',0, ...
    'HistWinFRNoiseCorr',0,'HistBinsY',[],'RatePercentile',75, ...
    'WinFRNoiseCorr',0,'WinFRNoiseScatter',0, ...
    'CLims',[],'PercentileWithZeros',0,'WinCVNoiseScatter',0,'CV',0, ...
    'ScatterWithInsigP',0,'ScatterDotSize',50,'Std',0,'WinStdNoiseScatter',0, ...
    'XLim',[],'YLim',[],'FontName','Helvetica','FontSize',10,'NoTitle',0, ...
    'TickDir','in','DiffThreshBins',0,'SpCorr',0,'NoDiff',0,'Norm',0, ...
    'Median',0);
Args.flags = {'HideTitle','ArgsOnly', ...
    'Hist','Max','VMag','VAng','Corr','YYMeanVMag','YYVAngMax', ...
        'YYMeanMax','YYVAngVMag','GroupedIntra','VAngCorr','NoiseCorr', ...
        'WinNoiseCorr','FRMean','HistWinFRNoiseCorr','WinFRNoiseCorr', ...
        'WinFRNoiseScatter','PercentileWithZeros','WinCVNoiseScatter','CV', ...
        'ScatterWithInsigP','Std','WinStdNoiseScatter','NoTitle', ...
        'DiffThreshBins','SpCorr','NoDiff','Norm','Median'};
[Args,varargin2] = getOptArgs(varargin,Args);

frdmeanlabel = 'Mean FR Difference (s/s)';
frdmaxlabel = 'Max FR Difference (s/s)';
vmaglabel = 'Mean Unit Vector Magnitude';
vanglabel = 'Vector Angle (degrees)';
corrlabel = 'PSTH Correlation';
spcorrlabel = 'Spearman PSTH Correlation';
noisecorrlabel = 'Noise Correlation';
winnoisecorrlabel = 'Windowed Noise Correlation';
meanratedifflabel = 'Mean Rate Differences (spikes/second)';
ratedifflabel = 'Rate Differences (spikes/second)';
frdmedianlabel = 'Median FR Difference (s/s)';

if(~isempty(Args.NumericArguments))
	% plot one group at a time
	n = Args.NumericArguments{1};
	% get the indices that go with set n
	indices = (obj.data.setIndex(n)+1):obj.data.setIndex(n+1);
	% get length of indices
	nindices = length(indices);
	% get the relevant indices for data like frmat and winstd, which has at
	% least 2 cols for each set
	colidx = (obj.data.frSetIndex(n)+1):obj.data.frSetIndex(n+1);
		% get the firing rates of the cells that go with set n
	% get firing rates for both cells
	frmat = obj.data.frmat(:,colidx);
	% save a copy of frmat with 0's so that FRMean and CV options work
	% properly
	frmat1 = frmat;
	if(~Args.PercentileWithZeros)
		% replace all 0's with nan's
		frmat(frmat==0) = nan;
	end
	% get size of frmat
	[frrows,frcols] = size(frmat);
	% compute prctile
	rthresh = prctile(frmat,Args.RatePercentile);
	% find bins that exceed threshold
	threshbins = sum(frmat > repmat(rthresh,frrows,1),2);
	% find bins that exceed threshold for all cells
	tbinsidx = threshbins > (frcols-1);
	if(Args.Hist)
		if(Args.DiffThreshBins)
			% only compute bins that exceed percentile threshold
			if(nindices>1)
				data = nanmean(obj.data.frdiff(indices,tbinsidx));
				xlabel(meanratedifflabel, ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			else
				data = obj.data.frdiff(indices,tbinsidx);
				xlabel(ratedifflabel, ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)			
			end
			ctrldata = obj.data.frdiffctrl(n,tbinsidx);
		else			
			% take the mean of the differences if there is more than 1
			% this is different from the frdmean field which is the mean
			% difference per frame
			if(nindices>1)
				data = nanmean(obj.data.frdiff(indices,:));
				xlabel(meanratedifflabel, ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			else
				data = obj.data.frdiff(indices,:);
				xlabel(ratedifflabel, ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			end
			ctrldata = obj.data.frdiffctrl(n,:);
		end
		% plot histogram
		if(~isempty(Args.HistBins))
			bins = Args.HistBins;
		else
            alldata = [data(:); ctrldata(:)];
			mindata = min(alldata);
			maxdata = max(alldata);
			if(~isempty(Args.HistBinSize))
				bins = mindata:Args.HistBinSize:maxdata;
			else
				bins = mindata:(maxdata-mindata)/Args.NumHistBins:maxdata;
			end
		end
		counts = histcie(data,bins);
        counts2 = histcie(ctrldata,bins);
		stairs(bins,counts)
        hold on
        stairs(bins,counts2,'r');
        hold off
		ylabel('Number of occurrences', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		set(gca,'XScale',Args.XScale,'YScale',Args.YScale)
		if(~Args.NoTitle)
	    	title(getDataDirs('ShortName','DirString',obj.data.setNames{n}), ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
    elseif(Args.HistWinFRNoiseCorr || Args.WinNoiseCorr ...
    		|| Args.WinFRNoiseCorr || Args.WinFRNoiseScatter ...
    		|| Args.WinCVNoiseScatter || Args.WinStdNoiseScatter)
		% length of ydata
		nbins = sum(~isnan(frmat));
    	% get noise correlation for all windows
		ydata = obj.data.winnoisecorr(:,n);
		% get p values for all windows
		ydatap = obj.data.winnoisecorrp(:,n);
		% get size of frmat
		[frmrows,frmcols] = size(frmat);
        % sort frmat to get indices we can sort again to convert rates to
        % percentile for each cell
        [frms1,frmsi1] = sort(frmat);
        [frms2,frmsi2] = sort(frmsi1);
        % normalize by the number of bins so values go between 0 and 100
        frmnorm = frmsi2 ./ repmat(nbins,frmrows,1) * 100;        
		% find windows that have significant p values
		sigbins = ydatap < Args.CorrPThresh;
		% get data for windows above threshold and significant
		idx1 = tbinsidx & sigbins;
		if(Args.ScatterWithInsigP)
			plotidx = tbinsidx;
		else
			plotidx = idx1;
		end
		if(Args.WinNoiseCorr || Args.WinFRNoiseCorr ...
			|| Args.WinFRNoiseScatter || Args.WinCVNoiseScatter ...
			|| Args.WinStdNoiseScatter)
			% get data for windows above threshold but not significant
			idx2 = tbinsidx & (~sigbins);
			% get data for windows below threshold but significant
			idx3 = (~tbinsidx) & sigbins;
			% get data for windows below threshold but not significant
			% ~A & ~B is equivalent to ~(A | B) which is ~idx1
			idx4 = ~idx1;
			if(Args.WinNoiseCorr)
				% create vector with 0's instead of nan's since the bar command
				% puts an asterisk in places where there are nan's
				ynan = zeros(frmrows,1);
				y1 = ynan;
				y1(idx1) = ydata(idx1);
				y2 = ynan;
				y2(idx2) = ydata(idx2);
				y3 = ynan;
				y3(idx3) = ydata(idx3);
				y4 = ynan;
				y4(idx4) = ydata(idx4);
				% create bar plot of these four vectors using histc
				tb = obj.data.timebins(:,n);
				% shift x values so that it will look like we used
				% bar(...,'histc')
				xt = tb + (tb(2)-tb(1))/2;
				b4 = bar(xt,y4,1,'c');
				hold on
				b3 = bar(xt,y3,1,'m');
				b2 = bar(xt,y2,1);
				b1 = bar(xt,y1,1,'r');
				hold off
				set(b1,'EdgeColor','none');
				set(b2,'EdgeColor','none');
				set(b3,'EdgeColor','none');
				set(b4,'EdgeColor','none');
				xlabel('Time (ms)', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
				ylabel(winnoisecorrlabel, ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			elseif(Args.WinFRNoiseScatter)
				% create scatter plot. Use ydata(plotidx,:) so that we can 
				% either just specify the values and allow the scatter 
				% function to map the color, or specify the colors directly if
				% CLim is specified.
                scatter(frmnorm(plotidx,1),frmnorm(plotidx,2),Args.ScatterDotSize,ydata(plotidx),'filled')
                if(~isempty(Args.CLims))
                    caxis(Args.CLims)
                end
				h = colorbar;
				set(h,'FontName',Args.FontName,'FontSize',Args.FontSize);
				exlim = isempty(Args.XLim);
				eylim = isempty(Args.YLim);
				if(exlim && eylim)
					alims = [Args.RatePercentile 100];
					xlim(alims)
					ylim(alims)
				else
					if(~exlim)
						xlim(Args.XLim);
					end
					if(~eylim)
						ylim(Args.YLim);
					end
				end
				xlabel('Rate 1 (%)', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
				ylabel('Rate 2 (%)', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			elseif(Args.WinCVNoiseScatter)
				% get firing rates for both cells
				cvmat = obj.data.winstd(:,(obj.data.frSetIndex(n)+1) ...
					: obj.data.frSetIndex(n+1)) ./ frmat;
				scatter(cvmat(plotidx,1),cvmat(plotidx,2),Args.ScatterDotSize,ydata(plotidx),'filled')
                if(~isempty(Args.CLims))
                    caxis(Args.CLims)
                end
				h = colorbar;
				set(h,'FontName',Args.FontName,'FontSize',Args.FontSize);
				if(~isempty(Args.XLim))
					xlim(Args.XLim)
				end
				if(~isempty(Args.YLim))
					ylim(Args.YLim)
				end
				xlabel('CV 1', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
				ylabel('CV 2', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			elseif(Args.WinStdNoiseScatter)
				% get firing rates for both cells
				stdmat = obj.data.winstd(:,(obj.data.frSetIndex(n)+1) ...
					: obj.data.frSetIndex(n+1));
				scatter(stdmat(plotidx,1),stdmat(plotidx,2),Args.ScatterDotSize,ydata(plotidx),'filled')
                if(~isempty(Args.CLims))
                    caxis(Args.CLims)
                end
				h = colorbar;
				set(h,'FontName',Args.FontName,'FontSize',Args.FontSize);
				if(~isempty(Args.XLim))
					xlim(Args.XLim)
				end
				if(~isempty(Args.YLim))
					ylim(Args.YLim)
				end
				xlabel('Standard Deviation 1 (Hz)', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
				ylabel('Standard Deviation 2 (Hz)', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			else % if(Args.WinFRNoiseCorr)
				% get mean firing rates
				% frm = power(obj.data.frdiffctrl(n,:),2);
                % get mean percentiles
                frm = mean(frmnorm,2);
				% plot mean FR versus windowed noise correlation
				plot(frm(idx4),ydata(idx4),'c.')
				hold on
				plot(frm(idx3),ydata(idx3),'m.')
				plot(frm(idx2),ydata(idx2),'.')
				plot(frm(idx1),ydata(idx1),'r.')
				hold off
				xlabel('Mean Firing Rate (Hz)', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
				ylabel(winnoisecorrlabel, ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			end
			if(~Args.NoTitle)
				title(getDataDirs('ShortName','DirString',obj.data.setNames{n}), ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			end
        else % if(Args.HistWinFRNoiseCorr)
			% set up 2D histogram
			% xdata = power(obj.data.frdiffctrl(n,:),2)';
			xdata = mean(frmnorm,2);
			if(isempty(Args.HistBins))
				% get max and min of data
				xmax = max(xdata);
				xmin = min(xdata);
				if(isempty(Args.HistBinSize))
					binsize = (xmax-xmin)/Args.NumHistBins;
				else
					binsize = Args.HistBinSize;
				end
				% make sure final bin includes max
				% xlastbin = ceil((xmax-xmin)/binsize)*binsize;
				% set up bins
				xbins = xmin:binsize:xmax;
			else
				xbins = Args.HistBins;
			end
			% create vector with nan
			y1 = repmat(nan,frmrows,1);
            % y1(idx1) = ydata(idx1);
            y1(tbinsidx) = ydata(tbinsidx);
			if(isempty(Args.HistBinsY))
				% get max and min of data
				ymax = max(y1);
				ymin = min(y1);
				ybinsize = (ymax-ymin)/Args.NumHistBins;
				% make sure final bin includes max
				% xlastbin = ceil((xmax-xmin)/binsize)*binsize;
				% set up bins
				ybins = ymin:ybinsize:ymax;
			else
				ybins = Args.HistBinsY;
			end
			if(isempty(ybins) || sum(isnan(ybins)))
				cla
			else
				% histn puts the 1st column in y so put [ydata xdata]
				n1 = histn([y1 xdata],concatenate(ybins,xbins)');
				imagesc(xbins,ybins,n1)
				% flip the y-axis
				set(gca,'YDir','normal');
				h = colorbar;
				set(h,'FontName',Args.FontName,'FontSize',Args.FontSize);
			end
			xlabel('Mean Rate (%)', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			ylabel(winnoisecorrlabel, ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			if(~Args.NoTitle)
				title(getDataDirs('ShortName','DirString',obj.data.setNames{n}), ...
					'FontName',Args.FontName,'FontSize',Args.FontSize)
			end
		end % if(Args.WinNoiseCorr)
	elseif(Args.FRMean)
		% get the time bins
        xt = obj.data.timebins(:,n);
        odfrdiff = obj.data.frdiff(indices,:);
		% plot the rates
        if(Args.DiffThreshBins)
            ntbinsidx = ~tbinsidx;
            xt(ntbinsidx) = nan;
            frmat1(ntbinsidx,:) = nan;
            odfrdiff(:,ntbinsidx) = nan;
        end
        stairs(xt,frmat1)
		hold on
        if(~Args.NoDiff)
            % plot the difference
            stairs(xt,odfrdiff,'k')
        end
		% plot the mean
		% stairs(xt,power(obj.data.frdiffctrl(n,:),2),'k')
		% get the x-axis limits
		xl = xlim;
		% plot the percentile rates
		line(xl',repmat(rthresh,2,1),'LineStyle',':')
		hold off
		xlabel('Time (ms)', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		ylabel('Rate (Hz)', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		if(~Args.NoTitle)
			title(getDataDirs('ShortName','DirString',obj.data.setNames{n}), ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
	elseif(Args.CV)
		% get the data
		cvmat = obj.data.winstd(:,colidx) ./ frmat1;
		% get the time bins
        xt = obj.data.timebins(:,n);
		% plot the rates
		stairs(xt,cvmat)
		xlabel('Time (ms)', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		ylabel('CV', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		if(~Args.NoTitle)
			title(getDataDirs('ShortName','DirString',obj.data.setNames{n}), ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
	elseif(Args.Std)
		% get the data
		stdmat = obj.data.winstd(:,colidx);
		% get the time bins
        xt = obj.data.timebins(:,n);
		% plot the rates
		stairs(xt,stdmat)
		xlabel('Time (ms)', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		ylabel('Standard Deviation (Hz)', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		if(~Args.NoTitle)
			title(getDataDirs('ShortName','DirString',obj.data.setNames{n}), ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
	else
		data = obj.data.frdiff(indices,:);
		ctrldata = obj.data.frdiffctrl(n,:);
		xt = obj.data.timebins(:,n);
		if(Args.DiffThreshBins)
			% get indices that we want to set to nan
			ntbinsidx = ~tbinsidx;
			data(:,ntbinsidx) = nan;
			ctrldata(ntbinsidx) = nan;
	        xt(ntbinsidx) = nan;
		end
		dlabel = 'Rate Difference (Hz)';
		% plot differences for each bin
		% check if there is only 1 difference
		if(nindices>1)
			% plot all the individual differences
			% plot(data','.-')
            stairs(xt,data');
			hold on
			% plot the mean of the differences
			rsum = nanmean(data);
		else
			rsum = data;
		end
		% plot the running sum of the differences
		% plot(rsum,'k.-');
        stairs(xt,rsum,'k');
        if(~isempty(ctrldata))
			hold on
			% plot the differences for the control data
			% plot(ctrldata,'m.-');
			stairs(xt,ctrldata,'m--');
			hold off
		end
		% put the groupname and the cummulative differences in the title
		xlabel('Time (ms)', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		ylabel(dlabel, ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		if(~Args.NoTitle)
			title([getDataDirs('ShortName','DirString',obj.data.setNames{n}) ' D = ' num2str(obj.data.frdmean(n))], ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
	end
else
	% get number of sets
	nsets = obj.data.numSets;
    % set ctrldata to [] so we only need to change it for valid options
    ctrldata = [];
	% get firing rates for all cells
	frmat = obj.data.frmat;
	if(Args.PercentileWithZeros)
		% get number of bins for each pair. Use frdiff instead of 
		% winnoisecorr since the latter has nan values in bins where there
		% are no spikes in one of the cells
		nbins = sum(~isnan(obj.data.frdiff'));
		% duplicate nbins since there are 2 spike trains for each value in
		% nbins
		nbins2 = reshape(repmat(nbins,2,1),1,[]);
	else
		% replace all 0's with nan's
		frmat(frmat==0) = nan;
		% find number of non-nan values in each column
		nbins2 = sum(~isnan(frmat));			
	end
	% get size of frmat
	[frmrows,frmcols] = size(frmat);
	rthresh = prctile(frmat,Args.RatePercentile);
    frthresh = frmat > repmat(rthresh,frmrows,1);
	% only works for pairs
	% save the odd and even column indices so we can use it again to
	% calculate the average rate percentile
	oddcols = 1:2:frmcols;
	evencols = 2:2:frmcols;
	% find all the odd columns and reshape into 1 column
	frtc1 = reshape(frthresh(:,oddcols),[],1);
	frtc2 = reshape(frthresh(:,evencols),[],1);		
	threshbins = reshape(sum([frtc1 frtc2],2),frmrows,[]);
	% find the bins that are above threshold for both cells
	tbinsidx = threshbins > 1;
	% find bins that are below threshold so we can use groupDirs for Mean
	% and Max data
	ntbinsidx = ~tbinsidx;
	
	% check if Norm is specified and if so modify frdiff and
	% frdiffctrl 
	if(Args.Norm)
		% find the max of each cell
		cellmax = max(frmat);
		% normalize to max of each cell
		frmnorm = frmat ./ repmat(cellmax,frmrows,1);
		% extract odd and even cols so we can do the diff easily
		frmn1 = frmnorm(:,oddcols);
		frmn2 = frmnorm(:,evencols);
		frdiffs = abs(frmn1 - frmn2);
		% convert frmn1 and frmn2 to single column and then take the mean 
		% and the square root before reshaping back to frmrows and frmcols/2
		frdiffctrl = reshape(sqrt(mean([frmn1(:) frmn2(:)],2)),frmrows,[]);
        frdmean = nanmean(frdiffs)';
        frdmax = max(frdiffs)';
        frdmedian = prctile(frdiffs,50)';
        frdctrlmean = nanmean(frdiffctrl)';
        frdctrlmax = max(frdiffctrl)';
        frdctrlmedian = prctile(frdiffctrl,50)';
		frdmeanlabel = 'Mean Normalized Difference Per Frame';
		frdmaxlabel = 'Max Normalized Difference';
	else
		frdiffs = obj.data.frdiff';			
		frdiffctrl = obj.data.frdiffctrl';
        frdmean = obj.data.frdmean;
        frdmax = obj.data.frdmax;
        frdmedian = prctile(frdiffs,50)';
        frdctrlmedian = prctile(frdiffctrl,50)';
        frdctrlmean = obj.data.frdctrlmean;
        frdctrlmax = obj.data.frdctrlmax;
	end
	cdata = obj.data.corrcoef;
	cdatap = obj.data.corrpvalue;
	scdata = obj.data.spcorrcoef;
	scdatap = obj.data.spcorrpvalue;
	vadata = obj.data.vang;
	vadatap = obj.data.vangp;
	if(Args.DiffThreshBins)
		% frdiff is orgainzed by row so we have to take the transpose
		% of ntbinsidx to make sure we set the right bins to nan
		frdiffs(ntbinsidx) = nan;
		frdiffctrl(ntbinsidx) = nan;
		nsetsi = 1:nsets;
		nsetsi2 = nsetsi * 2;
		nsetsi3 = nsetsi2 - 1;
		if(Args.Corr)
			% compute corrcoef only for bins that are jointly above the
			% percentile threshold for both cells
			% d = frmat;
			% set bins below threshold to nan
			% ntb2 = reshape(repmat(ntbinsidx,2,1),frmrows,[]);
			% d(ntb2) = nan;
			% initialize data, which will store the corr coeffs
			cdata = repmat(nan,nsets,1);
			cdatap = cdata;
			for seti = nsetsi
				% compute the signal correlation coefficient between cells
				% [ce1,cep1] = corrcoef(d(:,nsetsi3(seti):nsetsi2(seti)));
				[ce1,cep1] = corrcoef(frmat(tbinsidx(:,seti),nsetsi3(seti):nsetsi2(seti)));
				% get the non-redundant values
				ce = ce1(2,1);
				cep = cep1(2,1);
				cdata(seti) = ce;
				cdatap(seti) = cep;
			end
		elseif(Args.SpCorr)
			scdata = repmat(nan,nsets,1);
			scdatap = scdata;
			for seti = nsetsi
				% compute the signal correlation coefficient between cells
				% [ce1,cep1] = corrcoef(d(:,nsetsi3(seti):nsetsi2(seti)));
				[ce1,cep1] = spcorr(frmat(tbinsidx(:,seti),nsetsi3(seti):nsetsi2(seti)));
				% get the non-redundant values
				ce = ce1(2,1);
				cep = cep1(2,1);
				scdata(seti) = ce;
				scdatap(seti) = cep;
			end
		elseif(Args.VAng)
			vadata = repmat(nan,nsets,1);
			vadatap = vadata;
			for seti = nsetsi
				% compute the signal correlation coefficient between cells
				[vmag,vang,vangp] = vecsimilarity(frmat(tbinsidx(:,seti),nsetsi3(seti):nsetsi2(seti))',varargin{:});
				% get the non-redundant values
				vadata(seti) = vang;
				vadatap(seti) = vangp;
			end
		end
	end

	if(Args.VAngCorr)
		% find corrcoef that have p-values greater than the threshold
		sigind = cdatap < Args.CorrPThresh;
		insigind = ~sigind;
		plot(vadata(insigind),cdata(insigind),'.')
		hold on
		plot(vadata(sigind),cdata(sigind),'r.')
		hold off
		xlabel('Vector Angle (degrees)', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		ylabel('PSTH Correlation', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		if(~Args.NoTitle)
			title(['Data from ' num2str(nsets) ' pairs of cells'], ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
		set(gca,'TickDir',Args.TickDir,'FontName',Args.FontName, ...
			'FontSize',Args.FontSize)		
		return
	elseif(Args.HistWinFRNoiseCorr || Args.WinFRNoiseScatter ...
			|| Args.WinCVNoiseScatter || Args.WinStdNoiseScatter)
    	% get noise correlation for all windows
		ydata = obj.data.winnoisecorr;
		% get p values for all windows
		ydatap = obj.data.winnoisecorrp;
        % sort frmat to get indices we can sort again to convert rates to
        % percentile for each cell
        [frms1,frmsi1] = sort(frmat);
        [frms2,frmsi2] = sort(frmsi1);
        % normalize by the number of bins so values go between 0 and 100
        frmnorm = frmsi2 ./ repmat(nbins2,frmrows,1) * 100;
        % bins with nan that were padded by the plus function will have
        % values greater than 100 so remove those
        frmnorm(frmnorm>100) = nan;
		% find windows that have significant p values
		sigbins = ydatap < Args.CorrPThresh;
		% get data for windows above threshold and significant
		idx1 = tbinsidx & sigbins;
		% get the rate percentile for the 1st cell in a pair		
		frp1 = reshape(frmnorm(:,oddcols),[],1);
		% get the rate percentile for the 2nd cell in a pair		
		frp2 = reshape(frmnorm(:,evencols),[],1);
		if(Args.ScatterWithInsigP)
			plotidx = tbinsidx;
		else
			plotidx = idx1;
		end
		if(Args.HistWinFRNoiseCorr)
			% get the average percentile rate
			frpavg = reshape(mean([frp1 frp2],2),frmrows,[]);
			y1 = ydata(idx1);
			xdata = frpavg(idx1);
			if(isempty(Args.HistBins))
				% get max and min of data
				xmax = max(xdata);
				xmin = min(xdata);
				if(isempty(Args.HistBinSize))
					binsize = (xmax-xmin)/Args.NumHistBins;
				else
					binsize = Args.HistBinSize;
				end
				% make sure final bin includes max
				% xlastbin = ceil((xmax-xmin)/binsize)*binsize;
				% set up bins
				xbins = xmin:binsize:xmax;
			else
				xbins = Args.HistBins;
			end
			if(isempty(Args.HistBinsY))
				% get max and min of data
				ymax = max(y1);
				ymin = min(y1);
				ybinsize = (ymax-ymin)/Args.NumHistBins;
				% make sure final bin includes max
				% xlastbin = ceil((xmax-xmin)/binsize)*binsize;
				% set up bins
				ybins = ymin:ybinsize:ymax;
			else
				ybins = Args.HistBinsY;
			end			
			n1 = histn([y1 xdata],concat(ybins,xbins)');
			if(isempty(Args.CLims))
				imagesc(xbins,ybins,n1)
			else
				imagesc(xbins,ybins,n1,Args.CLims)
			end
			% flip the y-axis
			set(gca,'YDir','normal');
			h = colorbar;
			set(h,'FontName',Args.FontName,'FontSize',Args.FontSize);
			xlabel('Mean Rate (%)', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			ylabel(winnoisecorrlabel, ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		elseif(Args.WinFRNoiseScatter)
            % create scatter plot. Use ydata(plotidx,:) so that we can 
            % either just specify the values and allow the scatter 
            % function to map the color, or specify the colors directly if
            % CLim is specified.
			scatter(frp1(plotidx),frp2(plotidx),Args.ScatterDotSize,ydata(plotidx),'filled')
			if(~isempty(Args.CLims))
				caxis(Args.CLims)
			end
			h = colorbar;
			set(h,'FontName',Args.FontName,'FontSize',Args.FontSize);
			xlabel('Rate 1 (%)', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			ylabel('Rate 2 (%)', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		elseif(Args.WinCVNoiseScatter)
			% get data for windows above threshold and significant
			% idx2 = tbinsidx & (~sigbins);
            % calculate coefficient of variation
            cv1 = obj.data.winstd(:,oddcols) ./ frmat(:,oddcols);
            cv2 = obj.data.winstd(:,evencols) ./ frmat(:,evencols);
			scatter(cv1(plotidx),cv2(plotidx),Args.ScatterDotSize,ydata(plotidx),'filled')
			if(~isempty(Args.CLims))
				caxis(Args.CLims)
			end
			h = colorbar;
			set(h,'FontName',Args.FontName,'FontSize',Args.FontSize);
			xlabel('CV 1', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			ylabel('CV 2', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		else % if(Args.WinStdNoiseScatter)
			% get data for windows above threshold and significant
			% idx2 = tbinsidx & (~sigbins);
            % calculate coefficient of variation
            std1 = obj.data.winstd(:,oddcols);
            std2 = obj.data.winstd(:,evencols);
            % create scatter plot. Use ydata(plotidx,:) so that we can 
            % either just specify the values and allow the scatter 
            % function to map the color, or specify the colors directly if
            % CLim is specified.
			scatter(std1(plotidx),std2(plotidx),Args.ScatterDotSize,ydata(plotidx),'filled')
			if(~isempty(Args.CLims))
				caxis(Args.CLims)
			end
			h = colorbar;
			set(h,'FontName',Args.FontName,'FontSize',Args.FontSize);
			xlabel('Standard Deviation 1 (Hz)', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			ylabel('Standard Deviation 2 (Hz)', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
		if(~Args.NoTitle)
			title(['Data from ' num2str(nsets) ' pairs of cells'], ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
        if(~isempty(Args.XLim))
            xlim(Args.XLim)
        end
        if(~isempty(Args.YLim))
            ylim(Args.YLim)
        end
		set(gca,'TickDir',Args.TickDir,'FontName',Args.FontName, ...
			'FontSize',Args.FontSize)		
		return
	elseif(Args.YYMeanVMag)
		data1 = frdmean;
		data2 = obj.data.vmag;
		ylabel1 = frdmeanlabel;
		ylabel2 = vmaglabel;
	elseif(Args.YYVAngMax)
		data1 = vadata;
		data2 = frdmax;
		ylabel1 = vanglabel;
		ylabel2 = frdmaxlabel;
	elseif(Args.YYMeanMax)
		data1 = frdmean;
		data2 = frdmax;
		ylabel1 = frdmeanlabel;
		ylabel2 = frdmaxlabel;
	elseif(Args.YYVAngVMag)
		data1 = vadata;
		data2 = obj.data.vmag;
		ylabel1 = vanglabel;
		ylabel2 = vmaglabel;
	elseif(Args.Max)
		if(Args.DiffThreshBins)
			data = max(frdiffs)';
			ctrldata = max(frdiffctrl)';
		else
			data = frdmax;
			ctrldata = frdctrlmax;
		end
		dlabel = frdmaxlabel;
	elseif(Args.Median)
		if(Args.DiffThreshBins)
			data = prctile(frdiffs,50)';
			ctrldata = prctile(frdiffctrl,50)';
		else
			data = frdmedian;
			ctrldata = frdctrlmedian;
		end
		dlabel = frdmedianlabel;
	elseif(Args.VMag)
		data = obj.data.vmag;
		dlabel = vmaglabel;
	elseif(Args.VAng)
        % the GroupedIntra option expects all the data so separate
        % significant and non-significant values only if we are plotting
        % histograms
        if(Args.Hist)
            % find p-values that are not significant
            pi = vadatap < Args.CorrPThresh;
            data = vadata(pi);
            ctrldata = vadata(~pi);
        else
    		data = vadata;
        end
		dlabel = vanglabel;
	elseif(Args.Corr)
        % the GroupedIntra option expects all the data so separate
        % significant and non-significant values only if we are plotting
        % histograms
        if(Args.Hist)
			pi = cdatap < Args.CorrPThresh;
			data = cdata(pi);
			ctrldata = cdata(~pi);
        else % if(~Args.Hist)
            data = cdata;
        end
		dlabel = corrlabel;
	elseif(Args.SpCorr)
        % the GroupedIntra option expects all the data so separate
        % significant and non-significant values only if we are plotting
        % histograms
        if(Args.Hist)
            % find p-values that are not significant
            pi = scdatap < Args.CorrPThresh;
			data = scdata(pi);
            ctrldata = scdata(~pi);
        else
            data = scdata;
        end
		dlabel = spcorrlabel;
    elseif(Args.NoiseCorr)
        % the GroupedIntra option expects all the data so separate
        % significant and non-significant values only if we are plotting
        % histograms
        if(Args.Hist)
            % find p-values that are not significant
            pi = obj.data.noisecorrpvalue < Args.CorrPThresh;
			data = obj.data.noisecorr(pi);
            ctrldata = obj.data.noisecorr(~pi);
        else
            data = obj.data.noisecorr;
        end
        dlabel = noisecorrlabel;
    elseif(Args.WinNoiseCorr)
        % the GroupedIntra option expects all the data so separate
        % significant and non-significant values only if we are plotting
        % histograms
    	% get noise correlation for all windows
		ydata = obj.data.winnoisecorr;
        % length of ydata
        nbins = size(ydata,1);
        if(Args.Hist || Args.GroupedIntra)
			% get p values for all windows
			ydatap = obj.data.winnoisecorrp;
			% find windows that have significant p values
			sigbins = ydatap < Args.CorrPThresh;
			% get data for windows above threshold and significant
			idx1 = tbinsidx & sigbins;
			% get data for windows above threshold but not significant
			idx2 = tbinsidx & (~sigbins);
			data = obj.data.winnoisecorr(idx1);
            ctrldata = obj.data.winnoisecorr(idx2);
        else
            data = obj.data.winnoisecorr(tbinsidx);
        end
        dlabel = winnoisecorrlabel;
	else % default is to get the mean of the differences
        if(Args.DiffThreshBins)
			data = nanmean(frdiffs)';
			ctrldata = nanmean(frdiffctrl)';
		else
			data = frdmean;
			ctrldata = frdctrlmean;		
        end
		dlabel = frdmeanlabel;
	end
	
	% plot all data
	if(Args.YYMeanVMag || Args.YYVAngMax || Args.YYMeanMax || Args.YYVAngVMag)
		xdata = 1:nsets;
		[ax,h1,h2] = plotyy(xdata,data1,xdata,data2);
		set(get(ax(1),'Ylabel'),'String',ylabel1)
		set(get(ax(2),'Ylabel'),'String',ylabel2)
        set(h1,'Marker','.')
        set(h2,'Marker','.')
		xlabel('Set Number', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		if(~Args.NoTitle)
	        title(['PSTH differences for ' num2str(nsets) ' sets of cells'], ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
	elseif(Args.Hist)
		if(~isempty(Args.HistBins))
			bins = Args.HistBins;
		else
			alldata = [data; ctrldata];
			mindata = min(alldata);
			maxdata = max(alldata);
			if(~isempty(Args.HistBinSize))
				bins = mindata:Args.HistBinSize:maxdata;
			else
				bins = mindata:(maxdata-mindata)/Args.NumHistBins:maxdata;
			end
		end
		counts = histcie(data,bins);
		if(~isempty(ctrldata))
			counts2 = histcie(ctrldata,bins);
            if(Args.Corr || Args.SpCorr || Args.NoiseCorr || Args.WinNoiseCorr || Args.VAng)
                % shift x values so that it will look like we used
                % bar(...,'histc')
                bins2 = bins + (bins(2)-bins(1))/2;
                % do stacked bar plot so it is clear that it is the same
                % distribution with some values being insignificant
                % arrange counts2 in the first column so it will be plotted
                % at the bottom
                hb = bar(bins2,[counts2 counts],1,'stack');
                % reverse the colors so that insignificant values are in
                % red instead of everything else being in red
                set(hb(1),'FaceColor','r')
                set(hb(2),'FaceColor','b')
                % make sure that the tick marks look like we used histc as
                % well
                set(gca,'XTick',bins);
                xlim([bins(1) bins(end)])
            else
                % plot the real distribution using bar
                bar(bins,counts,'histc');
				hold on	
                % use stairs so it is clear that the control data is a
                % different distribution
				stairs(bins,counts2,'r')
				hold off
            end
        else
            % use bar plot for plots without control data since we don't
            % have to worry about visibility of the two plots
			bar(bins,counts,'histc');
		end
		xlabel(dlabel, ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		ylabel('Number of Cell Pairs', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		if(~Args.NoTitle)
			title(['Histogram for ' num2str(nsets) ' sets of cells. Mean ' num2str(nanmean(data)) ', Median ' num2str(prctile(data,50))], ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
    elseif(Args.GroupedIntra)
        if(Args.WinNoiseCorr)
        	% create matrix of nan's
        	wncmat = repmat(nan,size(obj.data.winnoisecorr));
        	% insert values above threshold and significant
        	wncmat(tbinsidx) = obj.data.winnoisecorr(tbinsidx);
        	% creat boxplot
        	boxplot(wncmat,1);
            % add the psth noise correlation values
            hold on
            % get x values
            xvals = 1:nsets;
            % find p-values that are significant
            pi = obj.data.noisecorrpvalue < Args.CorrPThresh;
            npi = ~pi;
            plot(xvals(pi),obj.data.noisecorr(pi),'ok');
            plot(xvals(npi),obj.data.noisecorr(npi),'om');
            hold off
            GroupedIntraPlot(obj,varargin2{:})
        	% add a column of nan's so we can easily add nan to the data
        	% passed to boxplot
        	% wncmat2 = [wncmat repmat(nan,nbins,1)];
        	% convert gind to vector
        	% gindvec = vecc(concatenate(gind,repmat(nan,2,1)));
        	% find occurrences of nan in gindvec
        	% gin = isnan(gindvec);
        	% replace occurrences of nan with column number we added above
        	% gindvec(gin) = nsets + 1;
        	% create boxplot
        	% boxplot(wncmat2(:,gindvec),1)
            % givl = gindrows + 2;
            % set(gca,'XTick',givl/2:givl:gindcols*(gindrows+2),'XTickLabel',1:gindcols)
            % cxlim = xlim;
            % set(gca,'XTick',5:5:cxlim(2))
			xlabel('Intra-Group Pairs', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			ylabel(dlabel, ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			if(~Args.NoTitle)
				title('Heterogeneity of intra-group pairs', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize);
			end
        else
            % group pairs according to site
            gind = groupDirs(obj);
            % sort the correlation coefficients within each site so the
            % plot is more readable
            % save the indices so we can use it later to mark the
            % insignificant values
            [gsdata,gsdi] = sort(nanindex(data,gind),'descend');            
			gdata = flipud(gsdata);
            gsdi2 = flipud(gsdi);
            [gdrows,gdcols] = size(gdata);
            % add 2 rows of nan's so we can easily put a little bit of
            % space between neighboring sites
            gdata = [gdata; repmat(nan,2,gdcols)];
            gdrows2 = gdrows + 2;
            gdsize = [gdrows2,gdcols];
            % find where the values are in gdata
            gidx = ~isnan(gdata);
            % count the number of pairs at each site
            gdsum = sum(gidx);
            % add 2 to account for the space between sites
            gdnum2 = gdsum + 2;
            gdc = 1:gdcols;
            gdmax = nanmax(gdata(:));
            fakeval = 2*gdmax;
            % put a fake value in gdata so it will be extracted and thus
            % create the spacing between sites
            % use 1.1 since correlation coefficients cannot exceed 1
            gdata(sub2ind(gdsize,gdsum+1,gdc)) = fakeval;
            gdata(sub2ind(gdsize,gdnum2,gdc)) = fakeval;
            % now find all non-nan values, which will include the 2 fake
            % values added to each site
            gidx2 = ~isnan(gdata);
            % create a similar matrix to store the x values so we can plot
            % neighboring sites with slightly different colors
            gdx = gdata;
            % put in the x values that will be used in the bar plot
            gdx(gidx2) = 1:sum(gdnum2);
            % plot the odd and even columns separately so we can use
            % different colors
            ocols = 1:2:gdcols;
            ecols = 2:2:gdcols;
            % get the odd columns first
            % get the non-nan values
            gidx2o = gidx2(:,ocols);
            % get the data (e.g. correlation coefficients)
            gdata1o = gdata(:,ocols);
            % pull out the values into a column vector
            gdata1 = gdata1o(gidx2o);
            % convert the fake values back to nan
            gdata1(gdata1==fakeval) = nan;
            % get the x values for plotting
            gdxo = gdx(:,ocols);
            % pull out the x values into a column vector
            gdx1 = gdxo(gidx2o);
            % now create the bar plot
            bar(gdx1,gdata1,1,'w');
            hold on
            % now do the same for even columns
            gidx2e = gidx2(:,ecols);
            gdata1e = gdata(:,ecols);
            gdata1 = gdata1e(gidx2e);
            gdata1(gdata1==fakeval) = nan;
            gdxe = gdx(:,ecols);
            gdx1 = gdxe(gidx2e);
            bar(gdx1,gdata1,1,'FaceColor',repmat(0.5,1,3));
            
			if(Args.Corr || Args.SpCorr || Args.NoiseCorr || Args.VAng)
				% get the p-values
				if(Args.Corr)
					pdata = nanindex(cdatap,gind);
				elseif(Args.SpCorr)
					pdata = nanindex(scdatap,gind);
				elseif(Args.NoiseCorr)
					pdata = nanindex(obj.data.noisecorrpvalue,gind);
                elseif(Args.VAng)
					pdata = nanindex(vadatap,gind);
                end
                % rearrange the p values since we sorted the data earlier
                pdata2 = pdata(sub2ind([gdrows,gdcols],gsdi2,repmat(gdc,gdrows,1)));
				% find the values above threshold, i.e. not significant
                % pad with two rows of 0's so that pnsig will be the
                % correct size
				pnsig = [(pdata2 >= Args.CorrPThresh); repmat(logical(0),2,gdcols)];
    			% displace the significant mark so it will be visible
                gdata2 = gdata + (gdata ./ abs(gdata) * (Args.DisplaceBarSig * gdmax));
				% plot asterisks for values that are not significant
                plot(gdx(pnsig),gdata2(pnsig),'.','MarkerSize',6)
			end
            hold off
			xlabel('Intra-Group Pairs', ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			ylabel(dlabel, ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
			if(~Args.NoTitle)
				title('Heterogeneity of intra-group pairs', ...
					'FontName',Args.FontName,'FontSize',Args.FontSize);
			end
		end
    else
		plot(data,'.-')
        if(~isempty(ctrldata))
            hold on
            plot(ctrldata,'r.-')
            hold off
        end
		ylabel(dlabel, ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		xlabel('Set Number', ...
			'FontName',Args.FontName,'FontSize',Args.FontSize)
		if(~Args.NoTitle)
			title(['PSTH differences for ' num2str(nsets) ' sets of cells'], ...
				'FontName',Args.FontName,'FontSize',Args.FontSize)
		end
	end
end
set(gca,'TickDir',Args.TickDir,'FontName',Args.FontName, ...
	'FontSize',Args.FontSize)
rvarl = length(Args.ReturnVars);
if(rvarl>0)
    % assign requested variables to varargout
    for rvi = 1:rvarl
        idx = (rvi-1)*2+1;
        varargout{1}{idx} = Args.ReturnVars{rvi};
        varargout{1}{idx+1} = eval(Args.ReturnVars{rvi});
    end
end
