function obj = plot(obj,varargin)
%@shufflesync/plot Plot function for shufflesync objects.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1, ...
	'SyncSpikes',0,'ZScoreHist',0,'ZScoreThresh',3,'HistStep',0.5, ...
	'ZScore',0,'SurrMeanStd',0,'Rasters',0,'PSTH',0, ...
	'XCorr',0,'PsthPThresh',0.05,'NoXCorrData',0,'NoPsthData',0, ...
	'Try',0,'YLim',[],'NormTest','','NormTestP',0.05, ...
    'CLim',[],'ColorBar',0,'Windowed',0,'AllZ3',0, ...
    'HighlightRasters',0,'XLim',[],'PSTHWindow',11,'SmoothedPSTH',0, ...
    'Z3Norm',0,'Z3Psth',0,'ExpandPSTH',0,'MinXCCount',10,'XCorrDataOnly',0, ...
    'NoTitle',0,'Summary',0,'SummaryFile','shufflesummary');
Args.flags = {'LabelsOff','SurrHist','ZScoreHist','ZScore', ...
	'SurrMeanStd','Rasters','PSTH','XCorr','NoXCorrData','NoPsthData', ...
	'Try','ColorBar','SyncSpikes','Windowed','AllZ3', ...
	'HighlightRasters','SmoothedPSTH','Z3Norm','Z3Psth', ...
    'ExpandPSTH','XCorrDataOnly','NoTitle','Summary'};
[Args,varargin2] = getOptArgs(varargin,Args, ...
    'remove',{'SyncSpikes','ZScoreHist','ZScoreThresh','HistStep', ...
        'ZScore','SurrMeanStd','Rasters','PSTH','XCorr','NoXCorrData'}, ...
    'aliases',{'PValue',{'PsthPThresh','NormTestP'}});

if(isempty(Args.NumericArguments))
	% use the first data set
	n = 1;
else
	n = Args.NumericArguments{1};
end
if(Args.SyncSpikes)
	% pad xc with a 0 so that the last data point will be drawn for the 
	% duration
	stairs(obj.data.WinStartTime(:,n),obj.data.DataSyncSpikes(:,n))
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
elseif(Args.Rasters)
	if(Args.Windowed && n>1)
        nwindows = size(obj.data.WinStartTime,1);
        % h = findobj(gca,'Color','m');
		xlim(obj.data.WinStartTime([max((n-2),1),min((n+4),nwindows)],1));
        xs = obj.data.WinStartTime(n,1);
        xe = obj.data.WinStartTime(min(n+2,nwindows),1);
        % set(h(1),'XData',[xe xe])
        % set(h(2),'XData',[xs xs])
        % line(repmat(obj.data.WinStartTime([n (n+2)],1)',2,1),repmat(ylim',1,2),'Color','m');
		return
	end
    cwd = pwd;
    cd(obj.nptdata.SessionDirs{n});
    % create nptgroup object
    ng = nptgroup('auto','GetClusterDirs');
    % plot adjspikes using nptgroup object and the list of highlighted
    % spikes
    if(Args.HighlightRasters)
	    plot(ng,'Object',{'adjspikes',{varargin2{:},'Highlight', ...
                obj.data.RasterMarkMatrices{n}}});
    else
	    plot(ng,'Object',{'adjspikes',{varargin2{:}}});
	end    
    cd(cwd);
	if(Args.Windowed && n==1)
		xlim(obj.data.WinStartTime([1,5],1));
        line(repmat(obj.data.WinStartTime([1 3],1)',2,1),repmat(ylim',1,2),'Color','m');
		return
	end
elseif(Args.XCorr)
    cwd = pwd;
    cd(obj.nptdata.SessionDirs{n});
    % load data file
    load(obj.data.Args.DataFile);
    % load surrogate data file
    load(obj.data.Args.SurrData);
    % length of lags
    lbl = length(lagbins);
    % difference in lagbins
    dlb = lagbins(2) - lagbins(1);
    % load result file
    load(obj.data.Args.ResultFile);
	winstime = obj.data.WinStartTime(:,n);
	xtime = winstime + (mean(diff(winstime)) ...
		* Args.GroupPlotIndex/(Args.GroupPlots+1));
	if(~Args.NoXCorrData)
        % fig = figure;
        set(gcf,'Renderer','painters');
        h = pcolor([winstime; endbins(end)],[lagbins lagbins(lbl)+dlb]-dlb/2,[[dxcorr zeros(lbl,1)]; zeros(1,numwindows+1)]);
        set(h,'EdgeColor','none');
        if(~isempty(Args.CLim))
            caxis(Args.CLim)
        end
        if(Args.ColorBar)
            colorbar
        end
	    hold on
	end
    if(~Args.XCorrDataOnly)        
        xt2 = repmat(xtime',lbl,1);
        lb2 = repmat(lagbins',1,numwindows);
        % find zscores that exceed the threshold
        zmat = abs(dxcz) > Args.ZScoreThresh;
        % check if we need to use the results of the normal distribution test
        if(~isempty(Args.NormTest))
            if(strcmp(Args.NormTest,'jbtest'))
                % use jbmat
                % ntmat = ~jbmat;
                ntmat = jbmat > Args.NormTestP;
            elseif(strcmp(Args.NormTest,'lillietest'))
                % use limat
                % ntmat = ~limat;
                ntmat = limat > Args.NormTestP;
            else
            	% find bins with counts greater than MinXCCount
            	ntmat = dxcorr > Args.MinXCCount;
            end
        else % if(~isempty(Args.NormTest))
            ntmat = ones(size(zmat));
        end % if(~isempty(Args.NormTest))
        % find windows where the p values for the difference in firing rate 
        % distributions are smaller than the threshold for both cells
        if(Args.SmoothedPSTH)
        	% we want to reject windows where there are two or more bins in
        	% each cell that exceeds the threshold
        	% so first transpose dspsthz so we can reshape to break the 
        	% z-scores into the 2 cells
        	dsp1 = reshape(dspsthz',windowlength,2*numwindows);
        	% now see how many bins are above zscore threshold for each cell
        	% and for each window
        	dsp2 = sum(dsp1>Args.ZScoreThresh) < 2;
        	% now reshape so we can check the result of both cells for each
        	% window
        	pvals = sum(reshape(dsp2,2,[])) > 1;
			% pvals = ~(sum(abs(dspsthz)>Args.ZScoreThresh,2) > 0)';
        else % if(Args.SmoothedPSTH)
			pvals = nansum(hpval>Args.PsthPThresh) > 1;
		end % if(Args.SmoothedPSTH)
		% psth difference is computed for each window while zmat contains a 
		% z-score for each lag so we need to use repmat to create a matrix that
		% is the same size as zmat
		pvmat = repmat(pvals,lbl,1);
        % find bins where the zscore is above threshold and the psth were
        % considered similar
        zpmat = zmat & pvmat;
        % find bins where the zscore is above threshold and the distribution
        % of the surrogates was normal
		wxcmat = zpmat & ntmat;
		if(Args.Summary)
			% find the significant bins
			[lb2i,wini] = find(wxcmat);
			lwxci = length(lb2i);
			if(Args.GroupPlotIndex>1)
				% load the data file
				load(Args.SummaryFile);
                % append the data in reverse order since we are going to
                % use unique and unique returns the last row for redundant
                % values, i.e. significant bins at multiple shuffle values
				shsummary = [lb2(lb2i) wini ...
					repmat(shuffle,lwxci,1); shsummary];
			else
				% create variables
				shsummary = [lb2(lb2i) wini repmat(shuffle,lwxci,1)];
			end
			% save to data file
			save(Args.SummaryFile,'shsummary');
			if(Args.GroupPlotIndex==Args.GroupPlots)
				ybins = 1:2:9;
                ybins2 = ybins+1; 
                xbins = Args.XLim(1):Args.XLim(2);
                xbins2 = xbins + 0.5;
                if(~isempty(shsummary))
					% find the smallest shuffle for each window
					[b,i,j] = unique(shsummary(:,[1 2]),'rows');
					% plot the data file
                    nsummary = histn(shsummary(i,[3 1]),concat(ybins',xbins','Columnwise'));
                    imagesc(xbins2,ybins2,nsummary);
                else
                    imagesc(xbins2,ybins2,zeros(length(ybins),length(xbins)))
                end
                if(~isempty(Args.CLim))
                    caxis(Args.CLim)
                end
			end % if(Args.GroupPlotIndex==Args.GroupPlots)
		else % if(Args.Summary)
			if(Args.AllZ3 || Args.Z3Norm)
				% find bins where the zscore is above threshold but the distribution of
				% the surrogates was not normal
				zpnot = zpmat & ~ntmat;
				% turn hold on again in case we are doing NoXCorrData and the first
				% hold on command was not issued
				hold on
				plot(xt2(zpnot),lb2(zpnot),'Marker','.','Color',[1 1 1], ...
					'LineStyle','none','MarkerSize',16);
			end
			if(Args.AllZ3 || Args.Z3Psth)
				% find windows where the zscore exceeded threshold but the psth were
				% considered different
				wxcmat2 = zmat & ntmat & (~pvmat);
				hold on
				plot(xt2(wxcmat2),lb2(wxcmat2),'Marker','.','Color',[0.5 0.5 0.5], ...
					'LineStyle','none','MarkerSize',16);
			end
			plot(xt2(wxcmat),lb2(wxcmat),'k.','MarkerSize',16)
			hold off
		end % if(Args.Summary)
    end % if(~Args.XCorrDataOnly)        
    if(~isempty(Args.XLim))
        xlim(Args.XLim)
    end
    if(~isempty(Args.YLim))
        ylim(Args.YLim)
    end
    if( Args.LabelsOff )
        set(gca,'XTickLabel',[])
    end
	cd(cwd);
elseif(Args.Windowed)
    cwd = pwd;
    ddir = obj.nptdata.SessionDirs{1};
    cd(ddir);
	surrprefix = obj.data.Args.SurrFilePrefix;
	l = load([surrprefix num2str(n,'%04d')]);
	if(Args.PSTH)
        if(~isempty(l.sptrains))
			% set some variables 
			numsurr = obj.data.Args.NumSurrogates;
            ns1 = numsurr + 1;
			defaultspt = repmat(nan,1,ns1);
			btdiff = obj.data.binSize;
			% convert cell arrays to matrices
			if(~isempty(l.sptrains{1}))
				sp1 = full(cell2mat(l.sptrains{1}));
			else
				sp1 = defaultspt;
			end
			if(~isempty(l.sptrains{2}))
				sp2 = full(cell2mat(l.sptrains{2}));
			else
				sp2 = defaultspt;
			end
			% load data from two windows before and two windows after
			% since the windows are stepped in half a frame
			l0sp1 = [];
			l0sp2 = [];
			if(n>2)
				l0 = load([surrprefix num2str(n-2,'%04d')]);
                if(~isempty(l0.sptrains))
                	if(~isempty(l0.sptrains{1}))
						l0spt1 = full(cell2mat(l0.sptrains{1}));
						l0sp1 = l0spt1(:,1);
					end
					if(~isempty(l0.sptrains{2}))
						l0spt2 = full(cell2mat(l0.sptrains{2}));
						l0sp2 = l0spt2(:,1);
					end
                end
			end
			l1sp1 = [];
			l1sp2 = [];
			if(n<(size(obj.data.WinStartTime,1)-1))
				l1 = load([surrprefix num2str(n+2,'%04d')]);
                if(~isempty(l1.sptrains))
                    if(~isempty(l1.sptrains{1}))
						l1spt1 = full(cell2mat(l1.sptrains{1}));
                        l1sp1 = l1spt1(:,1);
                    end
                    if(~isempty(l1.sptrains{2}))
						l1spt2 = full(cell2mat(l1.sptrains{2}));
                        l1sp2 = l1spt2(:,1);		
                    end
                end
			end
			% concatenate data together in order to perform histogram just
			% once
			smat = concat(sp1,sp2,l0sp1,l1sp1, ...
				l0sp2,l1sp2,'Columnwise');
            % get window size
            winsize = Args.PSTHWindow;
            ws1 = winsize - 1;
            % get half window size
            halfwinsize = ws1 / 2;
			% get the start of the frame before and the end of the frame
			% after
			d = load(obj.data.Args.DataFile);
            % get length of central window
            lcentralwin = d.windowlength;
            % compute length of psthbins and save for use later
            lpsthbins = lcentralwin + ws1;
            % shift by 0.5 at the beginning and the end so we can use
            % histcie instead of hist as hist includes all values outside
            % the range in the first and last bins
            psthbins = (d.startbins(n)-halfwinsize):(d.endbins(n)+halfwinsize);
			psthhistcbins = [psthbins-0.5 psthbins(lpsthbins)+0.5];
            % get indices representing window
            cwinidx = 1:lcentralwin;
            cwinidx2 = cwinidx+lcentralwin;
			psth = histcie(smat,psthhistcbins,'DropLast','DataCols');
			% add data from frame before and frame after before computing
			% smoothed PSTH
			ps1ind = 1:ns1;
			ps1ind2 = (2 * ns1) + [1 2];
			psth1 = psth(:,ps1ind) + repmat(sum(psth(:,ps1ind2),2),1,ns1);
			psth2 = psth(:,ps1ind+ns1) + repmat(sum(psth(:,ps1ind2+2),2),1,ns1);
			% now compute the sliding window smoothed psth
			smoothmtx = convmtx(ones(winsize,1),lpsthbins);
            % compute smoothed psth and then convert to rate
			spsth = smoothmtx * [psth1 psth2];
            % compute the relevant indices
            didx = cwinidx + ws1;
            % convert to bins so they can be used to label the axis
            % use subbins to convert from bin number to time
            bintimes = d.subbins;
            pbtimes = bintimes(psthbins(cwinidx+halfwinsize))';
            % shift by half a bin so the line plots of the data will be
            % centered in each bin
            pbtimes2 = pbtimes + btdiff/2;
            % pcolor takes bin edges rather than bin centers like imagesc
            xpsthbins = [pbtimes pbtimes(lcentralwin)+btdiff];
            % compute the scaling factor to convert to rate
            ratescale = 1000 / obj.data.reps / (winsize * btdiff);
            % grab the data
            dspsth1 = spsth(didx,1);
            dspsth2 = spsth(didx,ns1+1);
            % grab the surrogates
            surridx1 = 2:ns1;
            sspsthmat = [spsth(didx,surridx1)' spsth(didx,surridx1+ns1)'];
            % set the histogram bins
            aylim = Args.YLim;
            if(~isempty(aylim))
                sbins = aylim(1):aylim(2);
            else
                sbins = 0:max(sspsthmat(:));
            end
            % get the length of histogram bins
            sbl = length(sbins);
			% compute histogram of the surrogate psth for each bin
            % keep the last bin since pcolor needs one more bin at the end
            % to plot the value of the last data point
            sshist = histcie(sspsthmat,sbins);
            % compute mean and std of the surrogates
            sspm = mean(sspsthmat)';
            sspstd = std(sspsthmat)';
            % compute z-score for the data
            dspsthz = ([dspsth1; dspsth2] - sspm) ./ sspstd;
            % check if any are above 3
            zti = abs(dspsthz) > Args.ZScoreThresh;
            zti1 = zti(cwinidx);
            zti2 = zti(cwinidx2);
            % convert sbins to form useable as axis labels
            ysbins = ([sbins sbins(sbl)+1]-0.5) * ratescale;
            % create pad vectors so we don't have to generate them twice
            pad1 = zeros(sbl,1);
            pad2 = zeros(1,lcentralwin+1);
            % divide allocated axis into 2
            haxis = separateAxis('Vertical',2);
            set(gca,'Visible','off')
            % select first axis
            ha1 = axes('position',haxis(2,:));
            % plot distribution of surrogates for 1st cell
			h = pcolor(xpsthbins,ysbins,[[sshist(:,cwinidx) pad1]; pad2]);
			set(h,'EdgeColor','none');
			hold on
			% plot the data from cell 1
			dsrate1 = dspsth1 * ratescale;
            plot(pbtimes2,dsrate1,'k.-','MarkerSize',16)
            plot(pbtimes2(zti1),dsrate1(zti1),'ws','MarkerSize',16)
            hold off
            title(getDataDirs('ShortName','DirString',ddir))
            zoom on
			nwindows = size(obj.data.WinStartTime,1);
            if(Args.ExpandPSTH)
                if(n>1)
		            xlim(obj.data.WinStartTime([max((n-2),1),min((n+4),nwindows)],1));
		        end
                set(ha1,'XTick',[],'Color',[0 0 0.5])
            end
            if(Args.ColorBar)
    			colorbar
            end
            % select first axis
            ha2 = axes('position',haxis(1,:));
            % plot distribution of surrogates for 2nd cell
			h = pcolor(xpsthbins,ysbins,[[sshist(:,cwinidx2) pad1]; pad2]);
			set(h,'EdgeColor','none');
			hold on
			% plot the data from cell 2
			dsrate2 = dspsth2 * ratescale;
            plot(pbtimes2,dsrate2,'k.-','MarkerSize',16)
            plot(pbtimes2(zti2),dsrate2(zti2),'ws','MarkerSize',16)
            hold off
            zoom on
            if(Args.ExpandPSTH)
                if(n>1)
		            xlim(obj.data.WinStartTime([max((n-2),1),min((n+4),nwindows)],1));
		        end
                set(ha2,'Color',[0 0 0.5])
            end
            if(Args.ColorBar)
    			colorbar
            end
		else % if(~isempty(l.sptrains))
            cla
        end % if(~isempty(l.sptrains))
    else % if(Args.PSTH)
		s = load(obj.data.Args.SurrData);
		lagbins = s.lagbins;
		% length of lags
		lbl = length(lagbins);
		% difference in lagbins
		dlb = lagbins(2) - lagbins(1);
        xlagbins = [lagbins lagbins(lbl)+dlb]-dlb/2;
        if(~isempty(l.sptrains))
			surrxc = full(l.xc(:,3:1002))';
			xcdata = full(l.xc(:,1));
            % get ylim
            aylim = Args.YLim;
            if(~isempty(aylim))
                sxcbins = aylim(1):aylim(2);
            else
    			sxcbins = 0:max([surrxc(:); xcdata]);
            end
			sxcl = length(sxcbins);
			r = load(obj.data.Args.ResultFile);
			% get lags
			surrn = hist(surrxc,sxcbins);
			h = pcolor(xlagbins,[sxcbins sxcbins(sxcl)+1]-0.5,[[surrn zeros(sxcl,1)]; zeros(1,lbl+1)]);
			set(h,'EdgeColor','none');
			hold on
			% get z-scores for this window
			zsc = r.dxcz(:,n);
			% get indices that exceed threshold
			zsci = abs(zsc) > Args.ZScoreThresh;
			nzsci = ~zsci;
			% plot(lagbins(zsci),xcdata(zsci),'.','Color',[1 1 1])
			% plot(lagbins(nzsci),xcdata(nzsci),'.','Color',[0 0 0])
			plot(lagbins,xcdata,'k.-','MarkerSize',16)
			plot(lagbins(zsci),xcdata(zsci),'ws','MarkerSize',16)
			% plot(lagbins,l.xc(:,2),'o','Color',[0 0 0])
			hold off
			if(~isempty(Args.XLim))
				xlim(Args.XLim)
			end
            if(Args.ColorBar)
    			colorbar
            end
		else
            % do pcolor plot of empty array
            h = pcolor(xlagbins,-0.5:10.5,zeros(12,lbl+1));
            set(h,'EdgeColor','none');
        end
        title(getDataDirs('ShortName','DirString',ddir))
        zoom on
	end % if(Args.PSTH)
    cd(cwd);
    % return
elseif(Args.PSTH)
    cwd = pwd;
    cd(obj.nptdata.SessionDirs{n});
    % load data file
    load(obj.data.Args.DataFile);
    % load surrogate data file
    load(obj.data.Args.SurrData);
	bintimediff = subbins(2) - subbins(1);
	% get size of subbins
	sbl = length(subbins);
    % get conversion factor to go from mean spike counts to Hz
    pconvert = 1000 / bintimediff;
    % compute the PSTH for each cell
    ps1 = mean(binsc(:,1:100),2) * pconvert;
	ps2 = mean(binsc(:,101:200),2) * pconvert;
    psbins = subbins(1:(sbl-1));
	plot(psbins,ps1,'.-')
	hold on
	plot(psbins,ps2,'r.-')    
	% load the result file
	load(obj.data.Args.ResultFile);
	% psthbins will range from 1-shuffle to endbins(end)+shuffle
	% so pad subbins to add points in front and points after
	xt = [(-shuffle:-1)'*bintimediff; subbins; ...
		subbins(sbl)+( 1:(endbins(numwindows)-sbl+shuffle) )'*bintimediff];
	% plot surrogate means
	plot(nanindex(xt,psthbins+shuffle),psths{1}./reps*pconvert,'c.-')
	plot(nanindex(xt,psthbins+shuffle),psths{2}./reps*pconvert,'m.-')
	xlim([xt(1) xt(end)])
	hold off
	cd(cwd);
elseif(Args.Try)
		cd /Users/syen/Documents/ShihCheng/Data/Neural/Cat/newcatdata/t2/site02/session07/combinations/g3c2sg3c3s
	    sidx = 290;
		l = load(['shsysurr' num2str(sidx,'%04d')]);
		d = load('shsydata')
		s = load('shsysdata')
		surrxc = full(l.xc(:,3:1002))';
		sxcbins = 0:max(surrxc(:));
		n = hist(surrxc,sxcbins);		
	    xjb = ones(77,1);
	    xli = xjb;
		for i = 1:77
			[dummy,xjb(i)] = jbtest(surrxc(:,i));
			[dummy,xli(i)] = lillietest(surrxc(:,i));
		end
        pval = 0.01;
		xjbi = (xjb>pval);
		xlii = (xli>pval);	 
		figure
		plot(sxcbins,n(:,xjbi),'.-')
        title('Normal')
		figure
		plot(sxcbins,n(:,~xjbi),'.-')	
        title('Not Normal')
        figure
        imagesc(s.lagbins,sxcbins,n);
        hold on
        plot(s.lagbins(xjbi),l.xc(xjbi,1),'k.')
        plot(s.lagbins(~xjbi),l.xc(~xjbi,1),'ko')
        hold off
        figure
        plot(s.lagbins,xjb,'.-')
        hold on
        plot(s.lagbins(xjbi),xjb(xjbi),'r.')
        hold off
        title('p values from jbtest')
        figure
        i = 1;
        ev = event(1,77);
		while 1
            if(xjbi(i))
                plot(sxcbins,n(:,i),'r.-')
                title(['Lag ' num2str(s.lagbins(i)) ' Normal p = ' num2str(xjb(i))])
            else
                plot(sxcbins,n(:,i),'.-')
                title(['Lag ' num2str(s.lagbins(i)) ' Not Normal p = ' num2str(xjb(i))])
            end
			% get keyboard input to see what to do next
			key = input('RETURN - Next; p - Previous; N - N; q - Quit: ','s');
			i = str2num(key);
			if strcmp(key,'p')
				[i,ev] = Decrement(ev);
			elseif strcmp(key,'q')
				break;
			elseif ~isempty(i)
				ev = SetEventNumber(ev,i);
			else
				[i,ev] = Increment(ev);
			end
        end

        shuffle2 = 2;
        shuffle4 = 4;
        shuffle6 = 6;
        shuffle8 = 8;
		d2 = load('shsydata');
		l2 = load(['shsysurr' num2str(sidx,'%04d')]);
		d4 = load('shsydataA');
		l4 = load(['shsysurrA' num2str(sidx,'%04d')]);
		d6 = load('shsydataB');
		l6 = load(['shsysurrB' num2str(sidx,'%04d')]);
		d8 = load('shsydataC');
		l8 = load(['shsysurrC' num2str(sidx,'%04d')]);
        l2.psthbins = (d2.startbins(sidx)-shuffle2):(d2.endbins(sidx)+shuffle2);
        l2.psth1 = hist(cell2mat(l2.sptrains{1}),l2.psthbins);
        l4.psthbins = (d4.startbins(sidx)-shuffle4):(d4.endbins(sidx)+shuffle4);
        l4.psth1 = hist(cell2mat(l4.sptrains{1}),l4.psthbins);
        l6.psthbins = (d6.startbins(sidx)-shuffle6):(d6.endbins(sidx)+shuffle6);
        l6.psth1 = hist(cell2mat(l6.sptrains{1}),l6.psthbins);
        l8.psthbins = (d8.startbins(sidx)-shuffle8):(d8.endbins(sidx)+shuffle8);
        l8.psth1 = hist(cell2mat(l8.sptrains{1}),l8.psthbins);
	    % plot the raw data
	    figure
		plot(l2.psthbins,l2.psth1(:,1),'.-')
		hold on
		% get the mean of the surrogates
		pmean2 = mean(l2.psth1(:,2:1001),2);
		pmean4 = mean(l4.psth1(:,2:1001),2);
		pmean6 = mean(l6.psth1(:,2:1001),2);
		pmean8 = mean(l8.psth1(:,2:1001),2);
		% plot the surrogate means
		plot(l2.psthbins,pmean2,'r.-')
		plot(l4.psthbins,pmean4,'m.-')
		plot(l6.psthbins,pmean6,'g.-')
		plot(l8.psthbins,pmean8,'c.-')		
        hold off
		% get the std of the surrogates
		pstd2 = std(l2.psth1(:,2:1001),0,2);
		pstd4 = std(l4.psth1(:,2:1001),0,2);
		pstd6 = std(l6.psth1(:,2:1001),0,2);
		pstd8 = std(l8.psth1(:,2:1001),0,2);
        figure
        subplot(4,1,4)
		plot(l2.psthbins,l2.psth1(:,1),'.-')
		hold on
        errorbar(l8.psthbins,pmean8,pstd8,'c.-')
        hold off
        xl = xlim;
        subplot(4,1,1)
		plot(l2.psthbins,l2.psth1(:,1),'.-')
		hold on
        errorbar(l2.psthbins,pmean2,pstd2,'r.-')
        hold off
        xlim(xl)
        subplot(4,1,2)
		plot(l2.psthbins,l2.psth1(:,1),'.-')
		hold on
        errorbar(l4.psthbins,pmean4,pstd4,'m.-')
        hold off
        xlim(xl)
        subplot(4,1,3)
		plot(l2.psthbins,l2.psth1(:,1),'.-')
		hold on
        errorbar(l6.psthbins,pmean6,pstd6,'g.-')
        hold off
        xlim(xl)
		
		% compute smoothed psth
		l2mtx = convmtx(ones(11,1),39);
		l2p1 = l2mtx*l2.psth1;
		l2pbins = [((-5:-1) + l2.psthbins(1)) l2.psthbins (l2.psthbins(end) + (1:5))];
		l4mtx = convmtx(ones(11,1),43);
		l4p1 = l4mtx*l4.psth1;
		l4pbins = [((-5:-1) + l4.psthbins(1)) l4.psthbins (l4.psthbins(end) + (1:5))];
		l6mtx = convmtx(ones(11,1),47);
		l6p1 = l6mtx*l6.psth1;
		l6pbins = [((-5:-1) + l6.psthbins(1)) l6.psthbins (l6.psthbins(end) + (1:5))];
		l8mtx = convmtx(ones(11,1),51);
		l8p1 = l8mtx*l8.psth1;
		l8pbins = [((-5:-1) + l8.psthbins(1)) l8.psthbins (l8.psthbins(end) + (1:5))];
		% get the mean of the surrogates
		plmean2 = mean(l2p1(:,2:1001),2);
		plmean4 = mean(l4p1(:,2:1001),2);
		plmean6 = mean(l6p1(:,2:1001),2);
		plmean8 = mean(l8p1(:,2:1001),2);
		% plot the smoothed psth of the data
		figure
		plot(l2pbins,l2p1(:,1),'.-')
		hold on
		% plot the surrogate means
		plot(l2pbins,plmean2,'r.-')
		plot(l4pbins,plmean4,'m.-')
		plot(l6pbins,plmean6,'g.-')
		plot(l8pbins,plmean8,'c.-')	    
        hold off
		% get the std of the surrogates
		plstd2 = std(l2p1(:,2:1001),0,2);
		plstd4 = std(l4p1(:,2:1001),0,2);
		plstd6 = std(l6p1(:,2:1001),0,2);
		plstd8 = std(l8p1(:,2:1001),0,2);
        figure
        subplot(4,1,4)
		plot(l2.psthbins,l2.psth1(:,1),'.-')
		hold on
        errorbar(l8.psthbins,pmean8,pstd8,'c.-')
        hold off
        xl = xlim;
        subplot(4,1,1)
		plot(l2.psthbins,l2.psth1(:,1),'.-')
		hold on
        errorbar(l2.psthbins,pmean2,pstd2,'r.-')
        hold off
        xlim(xl)
        subplot(4,1,2)
		plot(l2.psthbins,l2.psth1(:,1),'.-')
		hold on
        errorbar(l4.psthbins,pmean4,pstd4,'m.-')
        hold off
        xlim(xl)
        subplot(4,1,3)
		plot(l2.psthbins,l2.psth1(:,1),'.-')
		hold on
        errorbar(l6.psthbins,pmean6,pstd6,'g.-')
        hold off
        xlim(xl)
        
		l20 = load(['shsysurrH' num2str(sidx,'%04d')]);
		m1s2 = cell2mat(l2.sptrains{1});
		v1s2 = m1s2(:,1);
		v2as2 = m1s2(:,2:end);
		v2s2 = v2as2(:);
		[hs2,ps2] = kstest2(v1s2,v2s2);
		[n1s2,bins1s2] = hist(v1s2,min(v1s2):max(v1s2));
		[n2s2,bins2s2] = hist(v2s2,min(v2s2):max(v2s2));
		m1s4 = cell2mat(l4.sptrains{1});
		v2as4 = m1s4(:,2:end);
		v2s4 = v2as4(:);
		[hs4,ps4] = kstest2(v1s2,v2s4);
		[n2s4,bins2s4] = hist(v2s4,min(v2s4):max(v2s4));
		m1s6 = cell2mat(l6.sptrains{1});
		v2as6 = m1s6(:,2:end);
		v2s6 = v2as6(:);
		[hs6,ps6] = kstest2(v1s2,v2s6);
		[n2s6,bins2s6] = hist(v2s6,min(v2s6):max(v2s6));
		m1s8 = cell2mat(l8.sptrains{1});
		v2as8 = m1s8(:,2:end);
		v2s8 = v2as8(:);
		[hs8,ps8] = kstest2(v1s2,v2s8);
		[n2s8,bins2s8] = hist(v2s8,min(v2s8):max(v2s8));
		m1s20 = cell2mat(l20.sptrains{1});
		v2as20 = m1s20(:,2:end);
		v2s20 = v2as20(:);
		[hs20,ps20] = kstest2(v1s2,v2s20);
		[n2s20,bins2s20] = hist(v2s20,min(v2s20):max(v2s20));
		figure
		plot(bins1s2,n1s2,'.');
        xlim([bins2s8(1) bins2s8(end)])
		a1 = gca;
		set(a1,'XColor','b','YColor','b')
		a2 = axes('Position',get(a1,'Position'), ...
			'XColor','k','YColor','k');
		plot(bins2s2,n2s2,'r.-','Parent',a2);
		hold on
		plot(bins2s4,n2s4,'m.-','Parent',a2);
		plot(bins2s6,n2s6,'g.-','Parent',a2);
		plot(bins2s8,n2s8,'c.-','Parent',a2);
		plot(bins2s20,n2s20,'k.-','Parent',a2);
		legend(['hs2=' num2str(hs2) ',ps2=' num2str(ps2)], ...
				['hs4=' num2str(hs4) ',ps4=' num2str(ps4)], ...
				['hs6=' num2str(hs6) ',ps6=' num2str(ps6)], ...
				['hs8=' num2str(hs8) ',ps8=' num2str(ps8)], ...
				['hs20=' num2str(hs20) ',ps20=' num2str(ps20)]);
        set(a2,'Color','none','YAxisLocation','right')
        xlim([bins2s8(1) bins2s8(end)])
		% [a,h1,h2] = plotyy(bins1s2,n1s2,bins2s2,n2s2);
		
	sidx = 203;
	d2 = load('shsydata');
	l2 = load(['shsysurr' num2str(sidx,'%04d')]);
	d4 = load('shsydataA');
	l4 = load(['shsysurrA' num2str(sidx,'%04d')]);
	d6 = load('shsydataB');
	l6 = load(['shsysurrB' num2str(sidx,'%04d')]);
	d8 = load('shsydataC');
	l8 = load(['shsysurrC' num2str(sidx,'%04d')]);
	smean2 = mean(l2.xc(:,3:end),2);
	sstd2 = std(l2.xc(:,3:end),0,2);
	figure
	subplot(3,1,1)
	errorbar(d2.lagbins,smean2,sstd2)
	hold on
	plot(d2.lagbins,l2.xc(:,2),'m.-')	
	plot(d2.lagbins,l2.xc(:,1),'r.-')	
	zs2 = (l2.xc(:,1)-smean2)./sstd2;
	zs2i = abs(zs2) > 3;
	plot(d2.lagbins(zs2i),l2.xc(zs2i,1),'r*')
	subplot(3,1,2)
	plot(d2.lagbins,zs2,'.')
	hold on
	plot(d2.lagbins(zs2i),zs2(zs2i),'r.')
	subplot(3,1,3)
	surrxc2 = full(l2.xc(:,3:end));
	sxcbins2 = 0:max(surrxc2(:));
	n2 = hist(surrxc2',0:max(surrxc2(:)));
	imagesc(d2.lagbins,sxcbins2,n2)
	
	smean4 = mean(l4.xc(:,3:end),2);
	sstd4 = std(l4.xc(:,3:end),0,2);	
	figure
	subplot(3,1,1)
	errorbar(d4.lagbins,smean4,sstd4)
	hold on
	plot(d4.lagbins,l4.xc(:,2),'m.-')	
	plot(d4.lagbins,l4.xc(:,1),'r.-')	
	zs4 = (l4.xc(:,1)-smean4)./sstd4;
	zs4i = abs(zs4) > 3;
	plot(d4.lagbins(zs4i),l4.xc(zs4i,1),'r*')
	subplot(3,1,2)
	plot(d4.lagbins,zs4,'.')
	hold on
	plot(d4.lagbins(zs4i),zs4(zs4i),'r.')
	subplot(3,1,3)
	surrxc4 = full(l4.xc(:,3:end));
	sxcbins4 = 0:max(surrxc4(:));
	n4 = hist(surrxc4',0:max(surrxc4(:)));
	imagesc(d4.lagbins,sxcbins4,n4)
	
	smean6 = mean(l6.xc(:,3:end),2);
	sstd6 = std(l6.xc(:,3:end),0,2);	
	figure
	subplot(3,1,1)
	errorbar(d6.lagbins,smean6,sstd6)
	hold on
	plot(d6.lagbins,l6.xc(:,2),'m.-')	
	plot(d6.lagbins,l6.xc(:,1),'r.-')	
	zs6 = (l6.xc(:,1)-smean6)./sstd6;
	zs6i = abs(zs6) > 3;
	plot(d6.lagbins(zs6i),l6.xc(zs6i,1),'r*')
	subplot(3,1,2)
	plot(d6.lagbins,zs6,'.')
	hold on
	plot(d6.lagbins(zs6i),zs6(zs6i),'r.')
	subplot(3,1,3)
	surrxc6 = full(l6.xc(:,3:end));
	sxcbins6 = 0:max(surrxc6(:));
	n6 = hist(surrxc6',0:max(surrxc6(:)));
	imagesc(d6.lagbins,sxcbins6,n6)
	
	smean8 = mean(l8.xc(:,3:end),2);
	sstd8 = std(l8.xc(:,3:end),0,2);	
	figure
	subplot(3,1,1)
	errorbar(d8.lagbins,smean8,sstd8)
	hold on
	plot(d8.lagbins,l8.xc(:,2),'m.-')	
	plot(d8.lagbins,l8.xc(:,1),'r.-')	
	zs8 = (l8.xc(:,1)-smean8)./sstd8;
	zs8i = abs(zs8) > 3;
	plot(d8.lagbins(zs8i),l8.xc(zs8i,1),'r*')
	subplot(3,1,2)
	plot(d8.lagbins,zs8,'.')
	hold on
	plot(d8.lagbins(zs8i),zs8(zs8i),'r.')
	subplot(3,1,3)
	surrxc8 = full(l8.xc(:,3:end));
	sxcbins8 = 0:max(surrxc8(:));
	n8 = hist(surrxc8',0:max(surrxc8(:)));
	imagesc(d8.lagbins,sxcbins8,n8)
	
    % get number of surrogate files
    numwindows = obj.data.WindowIndex(n+1) - obj.data.WindowIndex(n);
    % load the data file
    d = load(obj.data.Args.DataFile);
    % get the amount of shuffle
    shuffle = obj.data.Args.Shuffle;
    % get vector that picks out just the data without the shuffle edges
    sadd = shuffle + 1;
    % get the time of the subbins
    bintimes1 = d.subbins;
    % get length of bintimes1
    bt1l = length(bintimes1);
    % get binsize so we can add the necessary bins for the shuffles that 
    % go off the edge
    bsize = bintimes1(2)-bintimes1(1);
    % get extra bins
    extrabins = (1:shuffle)'*bsize;
    bintimes = [-flipud(extrabins); bintimes1; bintimes1(bt1l) + extrabins];
    dvec = sadd:(d.endbins(1)-d.startbins(1)+sadd);
    setholdflag = 1;
    % load the surrogate files
    for sidx = 1:numwindows
	    l = load([obj.data.Args.SurrFilePrefix num2str(sidx,'%04d')]);
        if(~isempty(l.sptrains))
            l.psthbins = (d.startbins(sidx)-shuffle):(d.endbins(sidx)+shuffle);
            % convert sptrains to matrices
            spt1 = cell2mat(l.sptrains{1});
            spt2 = cell2mat(l.sptrains{2});
            % check to see if we need to pad sptrains so the histogram is
            % taken in the right dimension
            if(size(spt1,1)==1)
                spt1 = concatenate(spt1,nan);
            end
            if(size(spt2,1)==1)
                spt2 = concatenate(spt2,nan);
            end
            l.psth1 = hist(spt1,l.psthbins);
            l.psth2 = hist(spt2,l.psthbins);
			% get the mean of the surrogates
			pmean1 = mean(l.psth1(:,2:1001),2);
			pmean2 = mean(l.psth2(:,2:1001),2);
            xvals = bintimes(l.psthbins(dvec)+shuffle);
            xvals2 = bintimes(l.psthbins+shuffle);
            if(rem(sidx,2))
                plot(xvals,l.psth1(dvec,1),'.-')
                if(setholdflag)
                    % turn hold on here instead of outside the loop since
                    % we want to be able to use the Overplot option from
                    % nptdata/plot
                    hold on
                    setholdflag = 0;
                end
                plot(xvals,l.psth2(dvec,1),'r.-')
				% plot the surrogate means
				plot(xvals2,pmean1,'c.-')
				plot(xvals2,pmean2,'m.-')
            else % if(rem(sidx,2))
				plot(xvals,l.psth1(dvec,1),'o-')
                if(setholdflag)
                    % turn hold on here instead of outside the loop since
                    % we want to be able to use the Overplot option from
                    % nptdata/plot. Need to have a check here since sidx=1
                    % might be empty so we might not always do the code in
                    % the if statement above.
                    hold on
                    setholdflag = 0;
                end
				plot(xvals,l.psth2(dvec,1),'ro-')
				% plot the surrogate means
				plot(xvals2,pmean1,'co-')
				plot(xvals2,pmean2,'mo-')
            end % if(rem(sidx,2))
        end % if(~isempty(l.sptrains))
	end % for sidx = 1:numwindows
else
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
	xc = obj.data.DataSyncSpikes(lagi,n);
    spxc = obj.data.ShiftPSyncSpikes(lagi,n);
	% find the non-nan values in xc which are the same as lagi
	% xc2 = xc(lagi);
    % spxc2 = spxc(lagi);
	xtime = lags2+(mean(diff(lags2))/2);
	% xtime = xtime(~isnan(xtime));
	imagesc(xtime,0:size(sxc2,1),sxc(sxci,:));
	hold on
	% find points with z-score larger than the threshold
	zi = abs(obj.data.DataZScore(lagi,n))>Args.ZScoreThresh;
    % find the points below threshold
    nzi = ~zi;
    plot(xtime,spxc,'go')
	plot(xtime(nzi),xc(nzi),'k*')
	plot(xtime(zi),xc(zi),'r*')
	hold off
	if(~Args.LabelsOff)
		xlabel('Time (ms)')
		ylabel('Number of synchronous spikes')
	end
end
if(~Args.NoTitle)
    title(getDataDirs('ShortName','DirString',obj.nptdata.SessionDirs{n}))
end

zoom on
