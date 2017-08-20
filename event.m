% load data
load jpdisplayevents

% set cell number
celln = @;

% number of reps
reps = 60;

% set psth bin size
psthbinsize = 3;

% number of points used to fit a gaussian to the central points of the 
% average xcross of the spike trains
ncpts = 60;

% find the maximum spike time for all repetitions
mi = zeros(reps,1);
sptimes = [];
for i = 1:reps
	% get max spike time
	spt = displayevents(celln).raster{i};
	mst = max(spt);
	% make sure max time is not empty
	if ~isempty(mst)
		mi(i) = mst;
	end
	sptimes = [sptimes; spt];
end
% make sure mimax is a multiple of bins size
mimax = ceil(max(mi)/psthbinsize) * psthbinsize;
% sort the spike times
sptimes = sort(sptimes);

% get bins for histc
bins = 0:1:mimax;
% get length of spike train
spl = mimax + 1;

% get spike trains for all repetitions
strain = zeros(spl,reps);
for i = 1:reps
	% make sure that there is at least 1 spike in each repetition
	if ~isempty(displayevents(celln).raster{i})
		s = histc(displayevents(celln).raster{i},bins);
		% make sure s is a column vector
		if (size(s,1) == 1)
			strain(:,i) = s';
		else
			strain(:,i) = s;
		end
	end
end

% get xcorr between all different pairs
% [stxc,lags] = xcorr(strain,'coeff');
stxc = zeros(2*spl-1,1);
for i = 1:reps
	for j = 1:reps
		if (i~=j)
			[c,lags] = xcorr(strain(:,i),strain(:,j));
			stxc = [stxc + c];
		end
		fprintf('%d,%d\n',i,j);
	end
end

% pick out center points
hncpts = ncpts/2;
cpoints = (spl-hncpts):(spl+hncpts);
xdata = lags(cpoints);
ydata = stxc(cpoints);

% pick out difference between center and peripheral value as starting 
% point for scale parameter in the curve fitting
yscale = ydata(hncpts+1)-ydata(1);
% use peripheral value for starting point for pedestal parameter in
% curve fitting
ypedestal = ydata(1);

% define function for curve fitting
f = inline('x(4) + x(3) * exp(-((xdata-x(1)).^2)/(2*x(2)^2))','x','xdata');

% do curve fit
x = lsqcurvefit(f,[0 6 yscale ypedestal],xdata,ydata);

% get gaussian values from the parameters of the curve fit
% divide standard deviation by sqrt(2)
y = gaussian([0 x(2)/sqrt(2) 1 0],xdata);

% subsample according to the psth resoution
gfilter = y(1:psthbinsize:(ncpts+1));

% normalize filter so that the sum is 1
gfilter = gfilter/sum(gfilter);

% get psth
% psth3 = jpcellspsth(1).psth{celln}.psth;
% compute the psth, divide by the number of reps and then convert to
% spikes per second
psth = histc(sptimes,0:psthbinsize:mimax) / (reps*psthbinsize) * 1000;
p = conv(psth,gfilter);

% get rid of points padded in the front and in the back
endpts = hncpts/psthbinsize;
ps = p((endpts+1):(end-endpts));

% take derivative
psd = diff(ps);

% find threshold crossings: ps3dx.falling are the peaks, and ps3dx.rising
% are the valleys
psdx = nptThresholdCrossings(psd,0,'separate','ignorefirst');

% set up time scale matching the 3 ms psth
pl = length(psth);
xt = (psthbinsize/2):psthbinsize:(pl*psthbinsize);

% get the peak values
pp = ps(psdx.falling);
pv = ps(psdx.rising);

% set values of valleys below 1 to 1 to prevent artificially inflating the pv ratio
pv1 = pv;
pv1(find(pv<1),1) = 1;
pv = sqrt(pp(1:(end-1)) .* pp(2:end)) ./ pv1;

% get mean of psth
pm = mean(ps);
% find peaks that are actually higher than the mean
pp1 = pp>pm;
% find consecutive peaks that are lower than the mean and take them out
% of pv to try to get pv distribution of real events
% shift p3p1 and add together. values that are zero represent two 
% consecutive peaks below the mean
pp2 = pp1(1:(end-1)) + pp1(2:end);
% change to binary so that 0 are two consecutive peaks below the mean
% and 1 where at least 1 of the peaks is above the mean
pp3 = pp2>0;
pv1 = pp3 .* pv;
% do ttest until we find mean that is definitely larger than pv distribution
% at p=0.05
% ttest(pv1,2.2,0.05,-1)
% get 95% confidence limit
[muhat,sigmahat,muci,sigmaci] = normfit(pv1(find((pv1>0) & (pv1<2.6))),0.05);
% find times for event boundaries
boundsi = psdx.rising(find(pv1>1.7));
bounds = xt(boundsi);

% do plots
stem(xt,psth,'.')
hold on
plot(xt,ps,'r')
plot(xt(psdx.falling),pp,'bs')
plot(xt(psdx.rising),pv,'ms')
% get max of psth
psthmax = max(psth);
% draw boundaries
line([bounds; bounds],[zeros(size(bounds)); ones(size(bounds))*psthmax])

% label pv ratios
% text( (xt3(ps3dx.rising(1:(end-1))))', (p3v(1:(end-1))+5)',num2str(pv','%.1f'))
text( xt(psdx.rising)', pv+5,num2str(pv1,'%.1f'))

% plot mean
line([0 10000],[pm pm],'Color',[0 0 0])

% plot spikes
% set offset
os = psthmax - reps;
for i = 1:reps
	plot(displayevents(celln).raster{i},os+i,'g.');
end

% subtract the mean from the smoothed psth to compute area between
% smoothed psth and the mean
pspm = ps - pm;
% remove the negative values
pspm(find(pspm<0)) = 0;

% for each event, sum the ps3pm
bl = length(boundsi) - 1;
% add 1 to the start of boundsi
boundsi = [1; boundsi];
% create memory for ps3pmArea
pspmArea = zeros(bl,1);
for i = 1:bl
	pspmArea(i) = sum(pspm(boundsi(i):boundsi(i+1)));
	% find the mean spike times in each event
	meanspt(i) = mean(sptimes(find((sptimes>xt(boundsi(i)) ) & (sptimes<=xt(boundsi(i+1))))));
end

% plot area on graph
text( (xt(boundsi(1:(end-1)))+1)', ones(bl,1)*(psthmax/2),num2str(pspmArea,'%.1f'))

% find areas that are larger than 70
eventt = meanspt(find(pspmArea>70));
% compute event intervals
eint = diff(eventt);
