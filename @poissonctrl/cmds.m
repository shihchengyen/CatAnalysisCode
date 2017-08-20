pc = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','VAng');

pcc = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','VAng','CellsVar','Events');

pcframe = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','VAng','FrameVar','Events');

pcframemax = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','VAng','FrameVar','Events','FanoMax');

pcframesum = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','VAng','FrameVar','Events','SumVar');

pcframesum1 = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','VAng','FrameVar','Events','SumVar', ...
    'StdMinFR',1,'RedoLevels',2,'Repetitions');

pcstd3 = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','VAng','StdMultiple',3,'Events');

pcstd6 = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','VAng','StdMultiple',6,'Events');

pcevents = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','Events','NoCorr');

pcpsparse = ProcessDays(poissonctrl,'IntraGroup','Cells',ndresp.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup','PSparse','NoCorr','Abs','UseNCells', ...
    'FrameVar','FanoMax');

[ndintra,pccpv] = ProcessDirs(ndintra,'Object','poissonctrl','FrameVar', ...
    'CellPopVar','RedoLevels',2,'Repetitions');

[ndintra2,pccpv] = ProcessDirs(ndintra2,'Object','poissonctrl','FrameVar', ...
    'CellPopVar','RedoLevels',1);

InspectGUI(pc)

figure; plot(pc)
pcobj = pcframemax;
pcobj = pcframesum;
pcobj = pcframesum1;
pcobj = pcstd6;

% plot heterogeneity of corr coef
load frdobjs
figure
plot(frdintra,'Corr','GroupedIntra');
hold on
plot(pcobj,'XShift',-0.5);
hold off

% plot vector separation of corr coef
figure
plot(frdintra,'VAng','GroupedIntra','NoSig');
hold on
plot(pcobj,'VAng','XShift',-0.5);
hold off

% plot heterogeneity of joint event probability 
load jointeventprob
figure
plot(jp);
hold on
plot(pcobj,'Events','XShift',-0.5);
hold off

a = sum(pcpsparse.data.pszscore>2)
b = sum(~isnan(pcpsparse.data.mfr(1:2:26)))
d = a ./ b
mean(d)
prctile(d,50)
e = prctile(pc.data.mfr,75)
m1 = pc.data.mfr;
m1(m1==0) = nan;
f = prctile(m1,75)
g = repmat(f,850,1);
h = m1 > g;
hs = sum((pc.data.pszscore>2) .* h);
mean(js ./ hs)

% plot poisson control data on top of ascorr data
plot(ascintra);
xvals = [1 3 5 6 8 11 14 18.5 23 25 27 30 32]; 
hold on
errorbar(xvals,pc.data.meancorr,pc.data.stdcorr)
hold off

% plot poisson control data on top of ascorr vang data
plot(ascintra,'VAng');
hold on
errorbar(xvals,pc.data.meanvang,pc.data.stdvang);
hold off

% pszscore field now stores percentage of surrogates with equal or smaller
% population sparseness values
% compute percentage of bins with percentage larger than 99
a = sum(pcpsparse.data.pszscore > 0.99) ./ sum(~isnan(pcpsparse.data.pszscore));
mean(a)
% find bins with means that are non-zero 
nzbins = pcpsparse.data.mfr(:,1:2:26)>0;
c = sum(pcpsparse.data.pszscore(nzbins) > 0.99) ./ sum()
% find average percentage of surrogates with PS values smaller than or 
% equal to PS value of the data
nanmean(pcpsparse.data.pszscore(:))
% find number of bins
sum(~isnan(pcpsparse.data.pszscore(:)))
% find bins in which at least one cell's firing rate is above the 75th percentile
load frdobjs
[frdintra,argsout] = plot(frdintra,'DiffThreshBins','ReturnVars', ...
    {'frdiffs','frdiffctrl','frthresh'});
frdiffs = argsout{2};
frthresh = argsout{6};
psp = logical(zeros(size(pcpsparse.data.pszscore)));
psp(:,1) = sum(frthresh(:,1:2),2)>0;
psp(:,2) = sum(frthresh(:,3:8),2)>0;
psp(:,3) = sum(frthresh(:,9:10),2)>0;
psp(:,4) = sum(frthresh(:,11:12),2)>0;
psp(:,5) = sum(frthresh(:,13:18),2)>0;
psp(:,6) = sum(frthresh(:,19:24),2)>0;
psp(:,7) = sum(frthresh(:,25:30),2)>0;
psp(:,8) = sum(frthresh(:,31:42),2)>0;
psp(:,9) = sum(frthresh(:,43:48),2)>0;
psp(:,10) = sum(frthresh(:,49:50),2)>0;
psp(:,11) = sum(frthresh(:,51:56),2)>0;
psp(:,12) = sum(frthresh(:,57:62),2)>0;
psp(:,13) = sum(frthresh(:,63:64),2)>0;
d = pcpsparse.data.pszscore(psp);
mean(d)

% mean percentage of bins with at least 1 cell with lower fano than
% the combined fano: 0.9942
mean(pccpv.data.mediancorr)

% mean ratio of combined fano to minimum fano: 1.3579
mean(pccpv.data.lqvang ./ pccpv.data.medianvang)

% find number of pairs with 1 autocorr lower than xcorr: 3
sum(pccpv.data.uqvang==1)
% find number of pairs with both autocorr lower than xcorr: 29
sum(pccpv.data.uqvang==2)

mfano = max(avgfano);
% find bins that are above 75th percentile
m1 = mean(smat(1:100,:));
m2 = mean(smat(101:200,:));
m1p = m1;
m2p = m2;
% replace 0's with nan
m1p(m1==0) = nan;
m2p(m2==0) = nan;
p1 = prctile(m1p,75);
p2 = prctile(m2p,75);
rbins = (m1>p1) | (m2>p1);
sum(rbins)
sum(popfano(rbins)>mfano(rbins))
rbins2 = (m1>p1) & (m2>p1);
sum(rbins2)
sum(popfano(rbins2)>mfano(rbins2))
s1 = std(smat(1:100,:));
s2 = std(smat(101:200,:));
stairs(popfano)
hold on
stairs(avgfano(1,:),'r')
stairs(avgfano(2,:),'g')
help stairs
stairs(m2>p2,'g:')
stairs((m1>p1)-0.5,'r:')
mfano2 = min(avgfano);
sum(mfano2(rbins2)<popfano(rbins2))
% check distributions of correlation coefficients between cells against
% distributions of correlation coefficients within cell
