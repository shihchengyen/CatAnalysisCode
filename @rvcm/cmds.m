% process data using original method
[ndresp,rvcmresp] = ProcessDirs(ndresp,'Object','rvcm','Original');
InspectGUI(rvcmresp,'Diff')
plot(rvcmresp,'MeanStd');
% process data using Jonathan's method
[ndresp,rvcmresp2] = ProcessDirs(ndresp,'Object','rvcm');
plot(rvcmresp2,'MeanStd');

figure
subplot(2,1,1)
% create boxplots of correlations computed using Yao's method
boxplot(rvcmobj.data.rvcmc)
subplot(2,1,2)
% create boxplots of correlations for all pairwise repetitions
boxplot(rvcmobj.data.repcc)