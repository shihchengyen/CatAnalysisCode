% function [cf,rfd,edges,hcounts,maxi,pf,pfpts] = spikeRecovery(intervals,varargin)
function [cf,rfd,edges,hcounts,maxi] = spikeRecovery(intervals,varargin)
%spikeRecovery Compute spike recovery function
%
%   [CF,RFD,EDGES,HCOUNTS,MAXI,PF] = spikeRecovery(INTERVALS,VARARGIN) 
%   computes the recovery function after a spike. The inter-spike 
%   interval (ISI) histogram is first computed and the peak of the 
%   histogram is found. The histogram counts after the peak are 
%   fitted to an exponential decay function to determine q, time 
%   constant of the exponential. The recovery function is computed 
%   for values leading up to the peak of the histogram by first 
%   converting the histogram values to probabilities by dividing 
%   through by the sum of the counts at all intervals. The recovery 
%   function at time t is then the probability of intervals at time t 
%   divided by q*(1-sum(ISI(1:peak))). This function is normalized by 
%   dividing all values by the value at the peak of the ISI so that 
%   the recovery function ends up at the value 1. For more details, 
%   see Berry&Meister,98. The recovery function is then fitted with a 
%   sigmoid function: 1./(1+exp(-a*(x-b))), with the first point in
%   the recovery function that crosses 0.5 as the starting point for
%   b and 1 as the starting point for a.
%
%   The optional arguments are as follows:
%      isibinsize - bin size of the histogram in ms (default 0.2)
%
%      pfpts - number of points after the peak to use to compute q
%              (default 4)
%
%      expfitmax - maximum interval used to find peak of histogram
%                  in ms (default 15)
%
%   These are the output arguments:
%      CF - curve fit object. Use getFunctionValues(cf) to return 
%           function values (length is not neccessarily equal to MAXI).
%      RFD - the recovery data (length is equal to MAXI)
%      EDGES - the histogram bins
%      HCOUNTS - the ISI histogram
%      MAXI - index of histogram bin corresponding to the peak
%      PF - array containing the slope of the exponential fit, followed
%           by the intercept in log units. The fit can be plotted
%           using: pfstart = maxi; pfend = maxi+pfpts; pfind =
%           pfstart:pfend;plot(edges(pfind),exp(pf(2)).*exp(pf(1)
%           *edges(pfind))
%
%   [CF,RFD,EDGES,HCOUNTS,MAXI,PF] = spikeRecovery(INTERVALS,VARARGIN) 

% default values for optional arguments
Args = struct('isibinsize',0.2, ... % step size of histogram
			  'pfpts',4, ... % number of points to use for polyfit
			  'expfitmax',15, ... % max of exponential fit in ms
			  'RunningFit',0, ... % do fits on sliding window
			  'FitSearchRange',[0 25], ... % look for FitDuration Window in this range
			  'FitDuration',5, ... % duration of window to do fit in ms
			  'Display',0,'Pause',0,'FitWindow',[],'MaxEndThreshold',0.3, ...
              'MaxPercent',0.95,'MinFitDurations',2,'MinSpikes',100); 
Args.flags = {'RunningFit','Display','Pause'};
% get optional arguments
Args = getOptArgs(varargin,Args);

isibinsize = Args.isibinsize;
pfpts = Args.pfpts;
expfitmax = Args.expfitmax;

edges = (0:isibinsize:max(intervals))';
% do histogram for all intervals since we need to do successive fits if the
% RunningFit option is specified and also when the meratio is below
% threshold in the default case
hcounts = histcie(intervals,edges);

% find the peak of the distribution within expfitmax
expfitmaxpts = round(expfitmax/isibinsize) + 1;
expfitpts = 1:expfitmaxpts;
pfrangedata = hcounts(expfitpts);

% check total number of spikes within the fit range
% tspikes = sum(pfrangedata);
% if(tspikes<Args.MinSpikes)
%     cf = [];
%     rfd = [];
%     maxi = [];
%     pf = [];
%     pfpts = [];
%     return;
% end

% if(Args.RunningFit)
% 	% find best fit to straight line starting from Args.FitSearchRange(1)
% 	% convert ms to data points
% 	fsRange = Args.FitSearchRange/Args.isibinsize;
% 	fsDuration = Args.FitDuration/Args.isibinsize;
% 	fitstart = fsRange(1) + 1;
% 	fitstop = fsRange(2) - fsDuration;
% 	% get logical vector indicating which values are non-zero
% 	hclogical = (hcounts~=0);
% 	% take the natural log of hcounts and not the log to base 10 since we
% 	% are plugging the value back into an exponential function. We would
% 	% have to scale the log10 value by log(10) in order for the fit to line
% 	% up with the data when we are plotting both at the same time
% 	hcounts2 = zerolog(hcounts);
% 	% initialize bestadjrs to inf
% 	bestadjrs = inf;
% 	for ind = fitstart:fitstop
% 		% get the indices
% 		hcind = ind:(ind+fsDuration-1);
% 		% get the edges
% 		edges2 = edges(hcind);
% 		% get the hcounts2
% 		hcounts3 = hcounts2(hcind);
% 		% use the hclogical to pull out just the values that are non-zero
% 		hcindvalid = hclogical(hcind);
% 		[fitmodel,goodness] = fit(edges2(hcindvalid),hcounts3(hcindvalid),'poly1');
% 		if(Args.Display)
% 			hold on
% 			plot(edges2(hcindvalid),hcounts3(hcindvalid),'.-')
% 			plot(fitmodel)
% 			legend off
% 			if(Args.Pause)
% 				key = input('RETURN - Next; q - Continue to End: ','s');
% 				n = str2num(key);
% 				if strcmp(key,'q')
% 					Args.Pause = 0;
% 				end
% 			end
% 		end
% 		if(goodness.adjrsquare<bestadjrs)
% 			% keep parameters for this last fit
% 			bestadjrs = goodness.adjrsquare;
% 			pf = fitmodel.p1;
% 			maxi = ind;
% 		end
% 	end
% elseif(~isempty(Args.FitWindow))
%     pfwind = Args.FitWindow/Args.isibinsize + 1;
%     pfind = pfwind(1):pfwind(2);
%     maxi = pfind(1);
%     pf = polyfit(edges(pfind),zerolog(hcounts(pfind)),1);
%     pfpts = pfind(2) - pfind(1) + 1;
% else
%     % take the log of the counts
% 	% take the natural log of hcounts and not the log to base 10 since we
% 	% are plugging the value back into an exponential function. We would
% 	% have to scale the log10 value by log(10) in order for the fit to line
% 	% up with the data when we are plotting both at the same time
    hcounts2 = zerolog(hcounts);    
	[maxh,maxi] = max(hcounts2(expfitpts));
% 	% find the difference between the values at the max and at expfitmax
%     % get a value at expfitmax that is not NaN, which is possible since we
%     % are taking the log. Take the next point until it is not NaN.
%     endindex = expfitmaxpts;
%     hcend = hcounts2(endindex);
%     while(isnan(hcend))
%         endindex = endindex - 1;
%         hcend = hcounts2(endindex);
%     end
% 	meratio = (maxh-hcend)/maxh;
%     if(isnan(meratio))
%         % return empty objects
%         cf = [];
%         rfd = [];
%         maxi = [];
%         pf = [];
%         pfpts = [];
%         return;    
%     elseif(meratio<Args.MaxEndThreshold)
%         maxperc = find(hcounts2>(Args.MaxPercent*maxh));
%         maxi = maxperc(1);
%         pfstart = maxi;
%         % do fit starting from maxi to 2 times the FitDuration (i.e. 10 ms, 
%         % 15 ms, etc.) until the exponential is negative
%         pfstep = Args.FitDuration/Args.isibinsize;
%         % set pfend to the larger of the next multiple of FitDuration 
%         % from the maximum or the minimum set in MinFitDurations
%         nfd = max(ceil(maxi/pfstep),Args.MinFitDurations);
%         pfend = nfd * pfstep + 1;
%         % make sure there is at least pfpts+1 points from the maximum to
%         % pfend
%         if((pfend-maxi)<pfpts)
%             % go to the next multiple of FitDuration to make sure we have
%             % enough points to do a fit
%             pfend = pfend + pfstep;
%         end
%         pf = [1 1];
%         while(pf(1)>0)
%             pfind = pfstart:pfend;
%             % don't add 1 to pfpts here because we compute pfend later in
%             % refractory/plot by doing pfstart+pfpts
%             pfpts = pfend - pfstart;
%             pf = polyfit(edges(pfind),hcounts2(pfind),1);
%             pfend = pfend + pfstep;
%         end
% 	else
% 		pfstart = maxi;
% 		pfend = maxi+pfpts;
%         pfind = pfstart:pfend;
% 		% take the natural log of hcounts and not the log to base 10 since we
% 		% are plugging the value back into an exponential function. We would
% 		% have to scale the log10 value by log(10) in order for the fit to line
% 		% up with the data when we are plotting both at the same time
% 		pf = polyfit(edges(pfind),hcounts2(pfind),1);
% 	end
% end

% compute recovery function using method of Berry & Meister, 98
% get the sum of hcounts
hs = sum(hcounts);
% convert hcounts to probabilty by dividing by hs
hprob = hcounts/hs;

% compute cummulative probablities
maxind = 1:maxi;
hcprob = tril(ones(maxi))*hprob(maxind);

% check to make sure pf(1) is not zero
% if(pf(1)==0 || isnan(pf(1)))
% 	cf = [];
%     rfd = [];
%     return;
% else
	% compute recovery function
	% pf was computed using time in milliseconds so need to multiply by
	% 1000 to get seconds and thus into units of Hz (counts per second)
	% wr = -1/(pf(1)*1000)*(hprob(1:maxi)./(1-hcprob));
    wr = hprob(maxind)./(1-hcprob);
	% normalize values so recovery function goes to 1
	rfd = wr/wr(end);
	
	% create fit for rfd
	ft = fittype('1./(1+exp(-a*(x-b)))','dependent',{'y'}, ...
					'independent',{'x'},'coefficients',{'a','b'});
	
	% find first point that crosses 0.5 to use as starting point in fitting
	% sigmoid curve to wr
	rfdi = find(rfd>0.5);
	% starting points for fit
	st = [1 edges(rfdi(1))];
	% do curve fit
	cf = fit(edges(maxind),rfd,ft,'Startpoint',st);
    if(Args.Display)
        plot(edges(maxind),rfd,'.')
        hold on
        plot(cf)
        hold off
        beep
        pause
    end
% end
