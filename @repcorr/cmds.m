cd /Users/syen/Documents/Work/Data/Neural/Cat/a4/site05/session25
g25 = ProcessSession(nptdata,'UnitType','s')
ng25 = ProcessSession(nptgroup,'IntraGroup','Cells',g25.SessionDirs,'AnalysisLevel','AllIntraGroup');
g25intra = ProcessSession(nptdata,'AnalysisLevel','pairs','Intragroup','Cells',g25.SessionDirs);
cd ../session28
g28 = ProcessSession(nptdata,'UnitType','s')
ng28 = ProcessSession(nptgroup,'IntraGroup','Cells',g28.SessionDirs,'AnalysisLevel','AllIntraGroup');
g28intra = ProcessSession(nptdata,'AnalysisLevel','pairs','Intragroup','Cells',g28.SessionDirs);
cd ../session29
g29 = ProcessSession(nptdata,'UnitType','s')
ng29 = ProcessSession(nptgroup,'IntraGroup','Cells',g29.SessionDirs,'AnalysisLevel','AllIntraGroup');
g29intra = ProcessSession(nptdata,'AnalysisLevel','pairs','Intragroup','Cells',g29.SessionDirs);
cd ../session31
g31 = ProcessSession(nptdata,'UnitType','s')
ng31 = ProcessSession(nptgroup,'IntraGroup','Cells',g31.SessionDirs,'AnalysisLevel','AllIntraGroup');
g31intra = ProcessSession(nptdata,'AnalysisLevel','pairs','Intragroup','Cells',g31.SessionDirs);

InspectGUI(ng28,'Object',{'adjspikes',{'Fast'}},'GroupEvent','GroupPlotSep','Vertical')
InspectGUI(ng29,'Object',{'adjspikes',{'Fast'}},'GroupEvent','GroupPlotSep','Vertical')
InspectGUI(ng31,'Object',{'adjspikes',{'Fast'}},'GroupEvent','GroupPlotSep','Vertical')

[g25intra,rc25] = ProcessDirs(g25intra,'Object','repcorr');
[g28intra,rc28] = ProcessDirs(g28intra,'Object','repcorr','Grating');
[g29intra,rc29] = ProcessDirs(g29intra,'Object','repcorr','Grating');

% compute intra-cell correlation coefficients as well as inter-cell
% correlation coefficients at the repetition level
[ndintra2,rcintra2] = ProcessDirs(ndintra2,'Object','repcorr','IntraInter');
% find the median of each distribution of correlation coefficients
md1 = prctile(rcintra2.data.repcorr,50);
% grab the intra-cell medians
md1a = md1(1:3:360);
md1b = md1(2:3:360);
% grab the inter-cell medians
md1c = md1(3:3:360);
% subtract the intra-cell medians from the inter-cell median
md1ac = md1a - md1c;
md1bc = md1b - md1c;
% plot the two differences
plot(md1ac,'.-')
hold on
plot(md1bc,'r.-')
% add the two differences to find the pairs that have the highest
% difference between the medians for the intra and inter distributions
md1abc = md1ac + md1bc;
[md1sabc,md1si] = sort(md1abc,2,'descend');
% create the proper indices to plot
ma = (md1si(1:10)-1)*3+1;
mb = [ma; ma+1; ma+2];
boxplot(data(:,mb(:)))
