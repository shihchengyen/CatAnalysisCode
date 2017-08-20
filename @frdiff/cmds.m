% compute objects
[ndintra,frdintra] = ProcessDirs(ndintra,'Object','frdiff');
[ndinter,frdinter] = ProcessDirs(ndinter,'Object','frdiff');

[ndintra,frdintra10] = ProcessDirs(ndintra,'Object','frdiff','AdjSpikes','redo','BinSize',10);
[ndinter,frdinter10] = ProcessDirs(ndinter,'Object','frdiff','AdjSpikes','redo','BinSize',10);

[ndintra,frdintra50] = ProcessDirs(ndintra,'Object','frdiff','AdjSpikes','redo','BinSize',50);
[ndinter,frdinter50] = ProcessDirs(ndinter,'Object','frdiff','AdjSpikes','redo','BinSize',50);

[ndintra,frdintra100] = ProcessDirs(ndintra,'Object','frdiff','AdjSpikes','redo','BinSize',100);
[ndinter,frdinter100] = ProcessDirs(ndinter,'Object','frdiff','AdjSpikes','redo','BinSize',100);

frdintra = frdintra100;
frdinter = frdinter100;
figure
plot(frdintra,'Hist','Corr','HistBins',-0.4:0.1:0.9);
figure
plot(frdinter,'Hist','Corr','HistBins',-0.4:0.1:0.9);
std(frdintra.data.corrcoef)
std(frdinter.data.corrcoef)
[r,p] = lillietest(frdintra.data.corrcoef)
[r,p] = lillietest(frdinter.data.corrcoef)
[r,p] = ttest2(frdintra.data.corrcoef,frdinter.data.corrcoef)

% check if distributions are different
p = anova1(concat(frdintra10.data.corrcoef,frdintra.data.corrcoef,frdintra50.data.corrcoef,frdintra100.data.corrcoef,'Columnwise'))
p = anova1(concat(frdinter10.data.corrcoef,frdinter.data.corrcoef,frdinter50.data.corrcoef,frdinter100.data.corrcoef,'Columnwise'))

p = kruskalwallis(concat(frdintra10.data.vang,frdintra.data.vang,frdintra50.data.vang,frdintra100.data.vang,'Columnwise'),[],'off')
p = kruskalwallis(concat(frdinter10.data.vang,frdinter.data.vang,frdinter50.data.vang,frdinter100.data.vang,'Columnwise'),[],'off')

% compute values for Poisson control
[frdintra,frdiffs,frdiffctrl] = plot(frdintra,'DiffThreshBins','ReturnVars', ...
    {'frdiffs','frdiffctrl'});

[frdintra,argsout] = plot(frdintra,'DiffThreshBins','ReturnVars', ...
    {'frdiffs','frdiffctrl','frthresh'});
[frdintra,argsout] = plot(frdintra,'ReturnVars', ...
    {'frdiffs','frdiffctrl','frthresh'});
frdiffsnt = argsout{2};
frthresh = argsout{6};

[frdintra,frdiffs,frdiffctrl] = plot(frdintra,'ReturnVars', ...
    {'frdiffs','frdiffctrl'});

% original calculation which is just the square root of the mean
frdc2 = frdiffctrl;
% frdiffctrl is the square root of the mean, i.e. sd
% but we need sqrt(2) * sd so multiply by sqrt2
frdc2 = frdiffctrl * sqrt(2);

% load matrix from poissonctrl calculation which contains the population
% mean for each site
load poissonctrl6
% take square root of the mean to get std and then multiple by sqrt2 since
% we are subtracting two Poisson processes which doubles the mean
pcmfr = sqrt(pc.data.mfr) * sqrt(2);

% load data from new poissonctrl calculation which uses the max fano for
% each frame to determine the variability
load poissonctrlobjs
% mfr fields contains mean and std in alternating columns
pcmfr = pcframemax.data.mfr(:,2:2:26) * sqrt(2);
% index into pcmfr to create appropriate comparison for frdiffs
frdc2 = pcmfr(:,[1 2 2 2 3 4 5 5 5 6 6 6 7 7 7 8 8 8 8 8 8 9 9 9 10 11 11 11 12 12 12 13]);

frd3sd = frdc2 * 3;
nb1 = frdiffs > frd3sd;
% find bins that are not NaN by cell
totalcell = sum(~isnan(frdiffs));
% find all bins that are not NaN
totalbins = sum(totalcell)
% find bins that are larger than 3 SD by cell
sdcell = sum(nb1);
% find all bins larger than 3 SD
sdbins = sum(sdcell)
% find percentage of bins
sdbins/totalbins
% find percentage of bins for each cell
percentcell = sdcell ./ totalcell;
% percentiles
prctile(percentcell,25)
prctile(percentcell,50)
prctile(percentcell,75)
% divide differences by std to get mean "z-score"
nb2 = frdiffs ./ frdc2;
nanmean(nb2(:))
nb3 = frdiffsnt ./ frdc2;
nanmean(nb3(:))

% mean percentage
mpc = mean(percentcell)
spc = std(percentcell)

% try out intra-inter correlation differences
% pairs sorted in IntraInterPairs.xls
% indices for intra and inter pairs are stored in first and second columns
% of a respectively
intracc = frdintra.data.corrcoef(a(:,1));
intercc = frdinter.data.corrcoef(a(:,2));
[h,p] = ttest(intracc,intercc);
% h = 0, p = 0.3398

% compute psth-cc for different window sizes
[ndintra2,frdintra10] = ProcessDirs(ndintra2,'Object','frdiff','AdjSpikes','redo','BinSize',10);
[ndinter2,frdinter10] = ProcessDirs(ndinter2,'Object','frdiff','AdjSpikes','redo','BinSize',10);

[ndintra2,frdintra50] = ProcessDirs(ndintra2,'Object','frdiff','AdjSpikes','redo','BinSize',50);
[ndinter2,frdinter50] = ProcessDirs(ndinter2,'Object','frdiff','AdjSpikes','redo','BinSize',50);

[ndintra2,frdintra100] = ProcessDirs(ndintra2,'Object','frdiff','AdjSpikes','redo','BinSize',100);
[ndinter2,frdinter100] = ProcessDirs(ndinter2,'Object','frdiff','AdjSpikes','redo','BinSize',100);


