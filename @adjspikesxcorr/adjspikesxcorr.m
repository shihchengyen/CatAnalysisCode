function obj = adjspikesxcorr(varargin)
%@adjspikesxcorr Constructor function for adjspikesxcorr class
%   OBJ = adjspikesxcorr('auto') attempts to create a adjspikesxcorr 
%   object by ...
% 
%   When the Sliding option is used, these are how the data fields are used:
%      numSets = 0; corresponds to number in nptdata
%      setNames = ''; corresponds to SessionDirs in nptdata
%      binSize = []; bin size used in dividing windows for synchrony calculation
%      lags = []; time points for windows in synchrony calculation
%      xc = []; number of synchronous spikes in each window
%      sxc = []; histogram of synchronous spikes in the surrogates for each window
%      ecf = []; z score for each window of the data compared to the surrogates
%      kspvalue = []; mean value of the surrogates for each window of the data
%      ksstat = []; standard deviation of the surrogates for each window of the data
%      legendstr = {};
%      fitgoodness = createEmptyFitGoodness;
%      fresult = {};
%      xcwidth = []; number of central bins used in computing synchronous spikes
%      pciCounts = []; keeps track of the frames in sxc that belong to each pair

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0,'ClusterDirs',{''}, ...
	'MaxLag',100,'BinSize',1,'KSAlpha',0.05,'NumCentralBins',1,'AutoCorr',0, ...
    'Sliding',0, ...
    'SubWindows',3,'SubWindowBins','frame/3','OverlapIncrBins','frame/2', ...
    'WaitForQsub',0,'QsubWaitTime',10,'JobWindows',20,'QsubProgram','asxcwrapper', ...
    'NumSurrogates',0,'SurrFilePrefix','asxcsurr','SurrHistMin',0, ...
    'SurrHistStep',1,'DataFile','asxcdata.mat');
Args.flags = {'Auto','AutoCorr','Sliding','WaitForQsub'};
[Args,varargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'adjspikesxcorr';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'adjxc';

if(nargin==0)
	% create empty object
	obj = createEmptyObject(Args);
elseif( (nargin==1) & isa(varargin{1},Args.classname))
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
		% check for saved object
		if(ispresent(Args.matname,'file','CaseInsensitive') ...
			& (Args.RedoLevels==0))
			fprintf('Loading saved %s object...\n',Args.classname);
			l = load(Args.matname);
			obj = eval(['l.' Args.matvarname]);
		else
			% no saved object so we will try to create one
			% pass varargin in case createObject needs to instantiate
			% other objects that take optional input arguments
			obj = createObject(Args,varargin{:});
		end
	end
end

function obj = createObject(Args,varargin)

if(Args.AutoCorr)
	adjs1 = adjspikes('auto',varargin{:});
	adjs2 = adjs1;
elseif(isempty(Args.ClusterDirs))
    % try finding ClusterDirs using current directory
    Args.ClusterDirs = getDataDirs('GetClusterDirs');
    if(~isempty(Args.ClusterDirs))
		% check that there are only 2 directories
		if(size(Args.ClusterDirs,2)>2)
			fprintf('Warning: %s works only on cell pairs! Using the first two cells...\n', ...
				Args.classname);
		end
		% get current directory
		cwd = pwd;
		% change to first directory
		cd(Args.ClusterDirs{1})
		% instantiate the adjspikes object instead of the ispikes object 
		% since we need to be able shift the repetitions with respect to 
		% each other to compute the significance values for the xcorr
		% pass varargin so things like RedoLevels and SaveLevels will work
		adjs1 = adjspikes('auto',varargin{:});
		% return to previous directory in case directory paths are relative
		cd(cwd);
		% change to second directory
		cd(Args.ClusterDirs{2});
		adjs2 = adjspikes('auto',varargin{:});
		% return to previous directory even if any of the objects are empty
		cd(cwd);
	else
		% no cluster dirs so create empty object
		obj = createEmptyObject(Args);
		return
	end
end

% check if any of the required objects are empty
if( isempty(adjs1) || isempty(adjs2) )
	% create empty object
	obj = createEmptyObject(Args);
else
	% this should be a valid object
	data.numSets = 1;
	data.setNames{1} = pwd;
	% save BinSize used to compute data
	data.binSize = Args.BinSize;
	if(Args.Sliding)
		% set some variables which will be used whether we are generating 
		% surrogates or not
		% get number of repetitions
		reps = size(adjs1.data.raster,2);
		% get frame points
		framepoints = adjs1.data.adjFramePoints;
		% get total number of frames over all repetitions
		frames = length(framepoints) - 1;
		% get number of frames in a repetition
		repframes = frames/reps;
		% first get frame points from the first repetition
		% need to subtract framepoints(1) since the first
		% framepoint could be non-zero
		fp2 = framepoints(1:(repframes+1)) - framepoints(1);
		% divide frame points by number of surrogate windows
		[subbins,binSize,nsbins] = divideBins(fp2,'SubBinSize', ...
			Args.BinSize);
		subwinbins = Args.SubWindowBins;
		nsbinstr = num2str(nsbins);
		if(ischar(subwinbins))
			% replace frame with nsbins
			sw2 = strrep(subwinbins,'frame',nsbinstr);
			% evaluate the string to get subwindow size and call ceil since
			% subwindow size has to be an integer and we want to make
			% sure we cover the entire frame
			SubWindowBins = ceil(eval(sw2));
		else
			SubWindowBins = ceil(subwinbins);
		end
		% get number of subwindows
		subwindows = Args.SubWindows;
		% get the length of the window
		windowlength = (subwindows * SubWindowBins);
		windowdiff = windowlength - 1;
		% length of bins used to compute xcorr
		lagslength = (2 * windowdiff) + 1;
		% get the bins corresponding to the beginning of a frame
		framestartbins = (0:(repframes-1))*nsbins + 1;
		overlap = Args.OverlapIncrBins;
		if(ischar(overlap))
			% replace frame with nsbins
			ov2 = strrep(overlap,'frame',nsbinstr);
			% evaluate the string to get overlap
			OverlapIncrBins = round(eval(ov2));
		else
			OverlapIncrBins = round(overlap);
		end
		% get the start bins for the overlaps
		overlapbins = 0:OverlapIncrBins:(nsbins-1);
		% get the start bins for all windows
		startbins0 = [ones(length(overlapbins),1) overlapbins'] ...
			* [framestartbins; ones(1,length(framestartbins))];
		startbins = startbins0(:);

		% check to see if there are surrogate files present
		surrfiles = nptDir([Args.SurrFilePrefix '*.mat']);
		if(isempty(surrfiles))
			% need to generate surrogates
			% use raster field since spike times are already separated
			% into repetitions so we can access the indices for each
			% repetition easily
			% concatenate raster field in both adjspikes objects
			raster12 = concatenate(adjs1.data.raster,adjs2.data.raster, ...
				'Columnwise');
			% take histogram			
			[binsc,binidx] = histcie(raster12,subbins,'DropLast','DataCols');
			% calculate 2*reps so the asxcSurrSynchrony function won't have
			% to compute it over and over again
			tworeps = reps * 2;
			% get the end bins
			endbins = startbins + windowdiff;
			% calculate total number of windows
			numwindows = length(startbins);
			% bins used to compute xcorr
			lagbins = -windowdiff:1:windowdiff;
			% get vector of incrementing subwindow bins, e.g. 1:12
			swbvector = (1:SubWindowBins)';
			% matrix of nan's used to store window bin numbers with spikes, 
			% e.g. size is 36x200
			wssize = [windowlength tworeps];
			windowspikes = repmat(nan,wssize);
			% matrix of nan's used to store subwindow bin numbers with spikes, 
			% e.g. size is 12x200x3
			swssize = [SubWindowBins tworeps subwindows];
			subwindowspikes = repmat(nan,swssize);
			% matrix of subwindow indices, e.g. size is [1:12]'x200x3
			subwindowindices = repmat(swbvector,[1 tworeps subwindows]);
			% matrix 1 used to compute xcorr, e.g. size is 36x2
			m1 = -ones(windowlength,2);
			% matrix 2 used to compute xcorr, e.g. size is 2x36
			m2 = ones(2,windowlength);
			% first column of xc matrix in output of asxcSurrSynchrony is
			% the xcorr of the actual data followed by NumSurrogates
			% columns of surrogate data
			dnssize = Args.NumSurrogates + 1;
			% matrix to store xcorr for data and surrogates, e.g. size is (35*2+1)x1001
			xctemp = zeros(lagslength,dnssize);
			% vector to store histogram of one particular xcorr
			nxctemp = zeros(lagslength,1);
			% vector contain indices from 2 to SurrSets+1 used to index into xc
			setnvec = 2:dnssize;
			% vector containing repetition indices for second cells, e.g. 101:200
			repnvec2 = (1:reps) + reps;
			% save spike counts so the function that generates surrogates
			% can access it
			save(Args.DataFile,'binidx','reps','tworeps','windowlength','SubWindowBins','subwindows', ...
				'startbins','endbins','lagbins','lagslength','swbvector','windowspikes','wssize', ...
				'subwindowspikes','swssize','subwindowindices','m1','m2','xctemp','dnssize', ...
				'nxctemp','setnvec','repnvec2')
				
			% save asxcdata binsc2 binidx2 winbinsize syncbins nsbins nsbins2 frames reps maxlag
			% call shell script to generate surrogates
			% create system command
			syscmd = [Args.QsubProgram ' ' num2str(numwindows) ' ' ...
				num2str(Args.JobWindows) ' ' Args.DataFile ' ' Args.SurrFilePrefix];
			% check if we should wait for the surrogates to be generated
			% don't try to use nohup of to background the job since we want
			% to take a look at the output to make sure the jobs were 
			% submitted properly
			[s,w] = system(syscmd);
			% make sure shell script ran without problems
			if(s~=0)
				error([Args.classname ': Error creating surrogates!'])
			else
				% display output
				fprintf('%s\n',w);
			end
			if(Args.WaitForQsub)
				% get last surrogate filename 
				lastsurr = [Args.SurrDataFile '.mat'];
				% check if last surrogate data has been generated
				while(~ispresent(lastsurr,'file'))
					% wait for 10 seconds
					pause(Args.QsubWaitTime)
				end
				% load surrogate data and fill in object's data structure
				surrfiles = nptDir([Args.SurrFilePrefix '*.mat']);
				data.lags = subbins(startbins);
				[data.xc,data.kspvalue,data.ksstat,data.sxc,data.ecf] ...
					= computeSurrData(Args,surrfiles,lagslength);
				% compute statistics on surrogates
			else
				% call system command with & so it gets backgrounded and 
				% returns immediately
				% [s,w] = system([syscmd]);
				data.sxc = [];
				data.ecf = [];
				data.kspvalue = [];
				data.ksstat = [];
			end
		else % if(isempty(surrfiles))
			% surrogates created now compute data
			data.lags = subbins(startbins);
			[data.xc,data.kspvalue,data.ksstat,data.sxc,data.ecf] ...
				= computeSurrData(Args,surrfiles,lagslength);
		end % if(isempty(surrfiles))			
		data.legendstr = {};
		data.fitgoodness = createEmptyFitGoodness;
		data.fresult = {};
		data.xcwidth = Args.NumCentralBins;
		data.pciCounts = [];
	else % if(Args.Sliding)
		% use histogram of spike time differences instead of xcorr
		% get rasters
		sp1 = adjs1.data.raster;
		% get number of repetitions and maximum number of spikes
		% for each raster. reps should be the same for both
		[nspikes1,reps] = size(sp1);
		% take the transpose of the second raster so we don't have to
		% take transpose in the for loop
		sp2 = adjs2.data.raster';
		nspikes2 = size(sp2,2);
		% get repetition indices
		rindices = 1:reps;
		% get shifted repetition indices
		shiftrindices = circshift(rindices,[0 -1]);
		% create matrices that will be used to compute all spike
		% time differences so we don't have to allocate memory
		% inside for-loop
		m1 = -ones(nspikes1,2);
		m2 = ones(2,nspikes2);
		% adjust MaxLag according to BinSize
		maxlags = ceil(Args.MaxLag/Args.BinSize);
		% set up bins
		lags = -maxlags:Args.BinSize:maxlags;
		% create empty arrays to store spike time differences within
		% the range of interest. Since we don't know how many points
		% there will be, we have no choice but to change array size
		% inside the for-loop
		xcvals = [];
		sxcvals = [];
		% initialize cell array
		for repn = 1:reps
			% replace data in matrices with spike times from each
			% repetition
			m1(:,1) = sp1(:,repn);
			m2(2,:) = sp2(repn,:);
			mxc = m1 * m2;
			% create matrices to compute shift predictor
			m2(2,:) = sp2(shiftrindices(repn),:);
			msxc = m1 * m2;
			% find values between -maxlags and maxlags
			mxci = find( (mxc>=-maxlags) & (mxc<=maxlags) );
			% find values between -maxlags and maxlags
			msxci = find( (msxc>=-maxlags) & (msxc<=maxlags) );
			% save those values so we can do histcie and kstest2 
			% later
			xcvals = [xcvals; mxc(mxci)];
			sxcvals = [sxcvals; msxc(msxci)];
		end
		% compute histcie on xc and sxc at the same time
		hcounts = histcie(concatenate(xcvals,sxcvals,'Columnwise'), ...
			lags,'DataCols');
		% separate xc and sxc
		xc = hcounts(:,1);
		sxc = hcounts(:,2);
		% get total number of spikes in 1st cell
		nspikes1 = sum(sum(~isnan(adjs1.data.raster)));
		% get total number of spikes in 2nd cell
		nspikes2 = sum(sum(~isnan(adjs2.data.raster)));
		% get total number of spikes fired by 
		nspikes = nspikes1 + nspikes2;
		if(Args.AutoCorr)
			% if computing auto-correlation, subtract out the 0 lag from the 
			% appropriate bin
			axcbin = maxlags + 1;
			xc(axcbin) = xc(axcbin) - nspikes1;
		end
		% get number of counts in bins -1:0 and 0:1
		centralEnd = maxlags + Args.NumCentralBins;
		centralBins = (maxlags-Args.NumCentralBins + 1):centralEnd;
		xccount = sum(xc(centralBins));
		sxccount = sum(sxc(centralBins));
		data.lags = lags;
		data.xc = xc;
		data.sxc = sxc;
		data.ecf = (xccount - sxccount)/nspikes;
		
		% do Kolmogorov-Smirnov test to see if the two distributions
		% are significantly different
		[ksh,data.kspvalue,data.ksstat] = kstest2(xcvals,sxcvals,Args.KSAlpha);
		
		% find width of correlation
		if(ksh==1)
			% distributions are different so take difference between xc 
			% and sxc
			dxc = xc - sxc;
			% find bins that are significant
			[phat,pci] = binofit(xc,nspikes);
			% multiply lower confidence limit by nspikes to get counts
			pciCounts = pci(:,1) * nspikes;
			data.pciCounts = pciCounts;
			% find bins where sxc is lower than pciCounts
			xclarge = sxc < pciCounts;
			% find difference in xclarge
			xcdlarge = diff(xclarge);
			% find positive changes
			xcdplus = find(xcdlarge==1);
			% find negative changes
			xcdminus = find(xcdlarge==-1);
			% find first value in xcdplus that is larger than the 
			% centralBin
			xcdiplus = find(xcdplus>centralEnd);
			% initialize emptyfit so if xcdiplus is empty things won't
			% crash
			emptyfit = 1;
			if(~isempty(xcdiplus) && xcdiplus(1)~=1)
				% set starting point for the gaussian fit to the value 
				% before xcdiplus and then add 1 since we took the diff and 
				% we need to shift the points forward by one to get the 
				% original indices
				gstart = xcdplus(xcdiplus(1)-1)+1;
				% find first value in xcdminus that is larger than the 
				% centralBin
				xcdiminus = find(xcdminus>centralEnd);
				% set end point for the gaussian fit to the point before the 
				% first value in xcdiminus but we also have to add 1 because 
				% we used diff so we just use the value in xcdminus
				gend = xcdminus(xcdiminus(1));
				% extract data between gstart and gend from the difference
				% between xc and sxc
				fitpts = gstart:gend;
				% need at least 3 points to do Gaussian fit with 3 parameters
				if(length(fitpts)>2)
					% fit dxc to gaussian
					% save fitgoodness in separate variable instead of directly
					% in data so we can keep the field order in data the same
					[cfobj,fitgoodness] = fit(data.lags(fitpts)', ...
						dxc(fitpts),'gauss1');
					% check the center of the gaussian to make sure it is
					% inside fitpts otherwise, it is probably not a valid
					% fit
					% get the center of the gaussian
					gcenter = cfobj.b1;
					if(gcenter>data.lags(gstart) && gcenter<data.lags(gend))
						gwidth = cfobj.c1;
						data.legendstr = {sprintf('Gaussian Width %f, Center: %f',gwidth,cfobj.b1)};
						% store fitgoodness
						data.fitgoodness = fitgoodness;
						data.fresult = {cfobj};
						data.xcwidth = gwidth;
						emptyfit = 0;
					else % if(gcenter>data.lags(gstart) && gcenter<data.lags(gend))
						data.legendstr = {'Gaussian center not within fit points!'};
					end % if(gcenter>data.lags(gstart) && gcenter<data.lags(gend))
				else % if(length(fitpts)>2)
					data.legendstr = {'Not enough fit points!'};
				end % if(length(fitpts)>2)
			else % if(~isempty(xcdiplus))
				data.legendstr = {'No fit points found!'};
			end % if(~isempty(xcdiplus))
			if(emptyfit)
				% fitdata is empty or center of gaussian is outside 
				% fitpts so set fresult to empty cfit object and widths 
				% to NaN. Do this here instead of outside the if(ksh==1)
				% block to keep pciCounts.
				data.fitgoodness = createEmptyFitGoodness;
				data.fresult = {cfit};
				data.xcwidth = NaN;
			end
		else % if(ksh==1)
			% distributions are not different so set fresult to empty 
			% cfit object and widths to NaN
			data.pciCounts = [];
			data.legendstr = {sprintf('Distributions not different!')};
			data.fitgoodness = createEmptyFitGoodness;
			data.fresult = {cfit};
			data.xcwidth = NaN;
		end % if(ksh==1)
	end % if(Args.Sliding)
	% create nptdata
	n = nptdata(data.numSets,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	if(Args.SaveLevels)
		fprintf('Saving %s object...\n',Args.classname);
		eval([Args.matvarname ' = obj;']);
		% save object
		eval(['save ' Args.matname ' ' Args.matvarname]);
	end
end % if( isempty(adjs1) | isempty(adjs2) )

function obj = createEmptyObject(Args)

data.numSets = 0;
data.setNames = '';
data.binSize = [];
data.lags = [];
data.xc = [];
data.sxc = [];
data.ecf = [];
data.kspvalue = [];
data.ksstat = [];
data.legendstr = {};
data.fitgoodness = createEmptyFitGoodness;
data.fresult = {};
data.xcwidth = [];
data.pciCounts = [];
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);

function r = createEmptyFitGoodness

r.sse = [];
r.rsquare = [];
r.dfe = [];
r.adjrsquare = [];
r.rmse = [];

function [data,surrmean,surrstd,surrhist,zscore] = computeSurrData(Args,surrlist, ...
	lagslength)

% get number of files to open
surrnum = length(surrlist);
% allocate memory for output arguments
data = zeros(surrnum,1);
sdata = zeros(Args.NumSurrogates,surrnum);
% create vector for summing appropriate bins in xcorr vector
% create vector using windowlength
sumvec = zeros(1,lagslength);
% find the center of lagslength
llcenter = ceil(lagslength/2);
centerbins = (Args.NumCentralBins-1)/2;
% create ones in places that correspond to NumCentralBins
sumvec((llcenter-centerbins):(llcenter+centerbins)) = 1;
for sidx = 1:surrnum
	l = load(surrlist(sidx).name);
	% sum counts within maxlags for data and surrogates by creating appropriate
	% vector
	nsyncspikes = sumvec * l.xc;
	% separate the surrogates from the data
	data(sidx) = nsyncspikes(1);
	sdata(:,sidx) = nsyncspikes(2:end)';
end
% compute histogram of surrogates
surrmean = mean(sdata)';
surrstd = std(sdata)';
% set bins so it covers the max of the data as well as of the surrogates
maxn = max([sdata(:); data(:)]);
bins = Args.SurrHistMin:Args.SurrHistStep:(ceil(maxn/Args.SurrHistStep) ...
	* Args.SurrHistStep);
% do hist since values are integers
surrhist = hist(sdata,bins);
% find bins with std of 0
si = find(surrstd);
% create z-score vector consisting of nans
zscore = repmat(nan,size(surrmean));
% fill bins that don't have std of 0
zscore(si) = (data(si) - surrmean(si)) ./ surrstd(si);
