[ndintra,ac] = ProcessDirs(ndintra,'Object','ascorr','redo');

InspectGUI(ac)

figure; plot(ac)
% show PSTH correlations
hold on
load frdobjs
plot(frdintra.data.corrcoef,'ms')
hold off

figure; plot(ac,'VAng')
% show PSTH correlations
hold on
load frdobjs
plot(frdintra.data.vang,'ms')
hold off

% create list in unix shell
cut -d "/" -f 1-2,4-5 movie-responsive.txt > moviecells.txt
cut -d "/" -f 1-2,4-5 mseq.txt > mseqcells.txt
comm -12 moviecells.txt mseqcells.txt > moviemseqcells.txt
% edit to insert session number for movie to create moviecellsmseq.txt
% edit to insert session number for mseq to create mseqcellsmovie.txt
% read in files
s = textread('moviecellsmseq.txt','%s');
% set up absolute paths
s1 = strcat([pwd '/'],s);
ndmoviemseq = nptdata('SessionDirs',s1);
s = textread('mseqcellsmovie.txt','%s');
% set up absolute paths
s1 = strcat([pwd '/'],s);
ndmseqmovie = nptdata('SessionDirs',s1);

ndmoviemseqintra = ProcessDays(nptdata,'Cells',ndmoviemseq.SessionDirs,'AnalysisLevel','Pairs','Intragroup');
ndmoviemseqinter = ProcessDays(nptdata,'Cells',ndmoviemseq.SessionDirs,'AnalysisLevel','Pairs','Intergroup');
ndmseqmovieintra = ProcessDays(nptdata,'Cells',ndmseqmovie.SessionDirs,'AnalysisLevel','Pairs','Intragroup');
ndmseqmovieinter = ProcessDays(nptdata,'Cells',ndmseqmovie.SessionDirs,'AnalysisLevel','Pairs','Intergroup');
ndmoviesparseintra = ProcessDays(nptdata,'Cells',ndmoviesparse.SessionDirs,'AnalysisLevel','Pairs','Intragroup');
ndsparsemovieintra = ProcessDays(nptdata,'Cells',ndsparsemovie.SessionDirs,'AnalysisLevel','Pairs','Intragroup');

[ndmoviemseqintra,asmoviemseqintra] = ProcessDirs(ndmoviemseqintra,'Object','ascorr','iSpikes','redo');
[ndmseqmovieintra,asmseqmovieintra] = ProcessDirs(ndmseqmovieintra,'Object','ascorr','iSpikes','redo');
[ndmoviemseqintra,asmoviemseqintra10] = ProcessDirs(ndmoviemseqintra,'Object','ascorr','iSpikes','redo','WindowSize',10);
[ndmseqmovieintra,asmseqmovieintra10] = ProcessDirs(ndmseqmovieintra,'Object','ascorr','iSpikes','redo','WindowSize',10);
[ndmoviemseqintra,asmoviemseqintra50] = ProcessDirs(ndmoviemseqintra,'Object','ascorr','iSpikes','redo','WindowSize',50);
[ndmseqmovieintra,asmseqmovieintra50] = ProcessDirs(ndmseqmovieintra,'Object','ascorr','iSpikes','redo','WindowSize',50);
[ndmoviemseqintra,asmoviemseqintra100] = ProcessDirs(ndmoviemseqintra,'Object','ascorr','iSpikes','redo','WindowSize',100);
[ndmseqmovieintra,asmseqmovieintra100] = ProcessDirs(ndmseqmovieintra,'Object','ascorr','iSpikes','redo','WindowSize',100);
[ndmoviemseqintra,asmoviemseqintrancc] = ProcessDirs(ndmoviemseqintra,'Object','ascorr','iSpikes','redo','NoCorrCoef');
[ndmseqmovieintra,asmseqmovieintrancc] = ProcessDirs(ndmseqmovieintra,'Object','ascorr','iSpikes','redo','NoCorrCoef');
[ndsparsemoviemseqintra,assparsemoviemseqintra] = ProcessDirs(ndsparsemoviemseqintra,'Object','ascorr','iSpikes','redo');
[ndmoviemseqinter,asmoviemseqinter] = ProcessDirs(ndmoviemseqinter,'Object','ascorr','iSpikes','redo');
[ndmseqmovieinter,asmseqmovieinter] = ProcessDirs(ndmseqmovieinter,'Object','ascorr','iSpikes','redo');
[ndsparsemovieintra,assparsemovieintra] = ProcessDirs(ndsparsemovieintra,'Object','ascorr','iSpikes','redo');
[ndmoviesparseintra,asmoviesparseintra] = ProcessDirs(ndmoviesparseintra,'Object','ascorr','iSpikes','redo');

a = asmoviemseqintra.data.contcoefs;
b = asmseqmovieintra.data.contcoefs;
figure
plot(a,b,'.')
hold on
line([-0.1 0.25],[-0.1 0.25],'Color','r')
hold off
[r,p] = ttest(a,b)
% r = 0; p = 0.4090
figure
hist(a-b)

a = asmoviemseqintra10.data.contcoefs;
b = asmseqmovieintra10.data.contcoefs;
figure
plot(a,b,'.')
hold on
line([-0.1 0.3],[-0.1 0.3],'Color','r')
hold off
[r,p] = ttest(a,b)
% r = 0; p = 0.2036
figure
hist(a-b)

a = asmoviemseqintra50.data.contcoefs;
b = asmseqmovieintra50.data.contcoefs;
figure
plot(a,b,'.')
hold on
line([-0.1 0.3],[-0.1 0.3],'Color','r')
hold off
[r,p] = ttest(a,b)
% r = 0; p = 0.3653
figure
hist(a-b)

a = asmoviemseqintra100.data.contcoefs;
b = asmseqmovieintra100.data.contcoefs;
figure
plot(a,b,'.')
hold on
line([-0.1 0.4],[-0.1 0.4],'Color','r')
hold off
[r,p] = ttest(a,b)
% r = 0; p = 0.4000
figure
hist(a-b)

% turns out to be exactly the same as using corrcoef
a = asmoviemseqintrancc.data.contcoefs;
b = asmseqmovieintrancc.data.contcoefs;
figure
plot(a,b,'.')
hold on
line([-0.1 0.25],[-0.1 0.25],'Color','r')
hold off
[r,p] = ttest(a,b)
% r = 0; p = 0.1980
figure
hist(a-b)

% get correlations for sparse noise stimuli
moviemseqidx = [1 2 3 20 21 22];
a = asmoviemseqintra.data.contcoefs(moviemseqidx);
b = asmseqmovieintra.data.contcoefs(moviemseqidx);
c = assparsemoviemseqintra.data.contcoefs;
plot(a,b,'.')
hold on
plot(a,c,'ro')
line([-0.02 0.2],[-0.02 0.2],'Color','g')
hold off
legend('M-sequence','Sparse Noise',0)
[r,p] = ttest(a,c)
% r = 0, p = 0.4503

% check spike count std
[ndmseqmovieintra,asstd] = ProcessDirs(ndmseqmovieintra,'Object','ascorr','iSpikes','redo','CheckStd');

a = asmoviemseqinter.data.contcoefs;
b = asmseqmovieinter.data.contcoefs;
figure
plot(a,b,'.')
hold on
line([-0.1 0.25],[-0.1 0.25],'Color','r')
hold off
[r,p] = ttest(a,b)
% r = 0; p = 0.3421
figure
hist(a-b)

% check correlations for pairs with visible kernel structure
idx = logical(zeros(23,1));
xidx = 1:23;
idx([3 4 5 6 13 14 15 17 18 19 20]) = 1;
nidx = ~idx;
a = asmoviemseqintra.data.contcoefs;
b = asmseqmovieintra.data.contcoefs;
figure
plot(xidx(idx),b(idx),'.')
hold on
plot(xidx(nidx),b(nidx),'rx')
hold off
figure
plot(a(idx),b(idx),'.')
hold on
plot(a(nidx),b(nidx),'m.')
hold off
[r,p] = ttest(a(idx),b(idx))
% r = 0; p = 0.5217

a = asmoviesparseintra.data.contcoefs;
b = assparsemovieintra.data.contcoefs;
figure
plot(a,b,'.')
hold on
line([-0.1 0.5],[-0.1 0.5],'Color','r')
hold off
[r,p] = ttest(a,b)
% r = 0; p = 0.5965
figure
hist(a-b)

[ndmoviemseqintra,acmoviemseq] = ProcessDirs(ndmoviemseqintra,'Object','ascorr','redo');
p = repmat(nan,200,46);
p(:,1:2:46) = acmoviemseq.data.ascorr;
p(1:102,2:2:46) = asmseqreps.data.ascorr;
boxplot(p,1)
% set breakpoint in nptdata/GroupedIntraPlot and do bpx = bpx*2 - 0.5
% do kstest2 to see which stimuli gives lower correlations
ksr = repmat(nan,23,6);
for idx = 1:23
    [ksr(idx,1) ksp(idx,2)] = kstest2(asmseqreps.data.ascorr(:,idx),acmoviemseq.data.ascorr(:,idx),0.05,1);
    [ksr(idx,3) ksp(idx,4)] = kstest2(asmseqreps.data.ascorr(:,idx),acmoviemseq.data.ascorr(:,idx));
    [ksr(idx,5) ksp(idx,6)] = kstest2(asmseqreps.data.ascorr(:,idx),acmoviemseq.data.ascorr(:,idx),0.05,-1);
end
% 19 out of 23 pairs are different, 
% mseq values are smaller in pairs 1 5 6 8 9 11 12 16 20 21 22
% movie values are smaller in pairs 2 3 4 7 15 16 17 19 23
% range of correlations for pair 16 is much larger for mseq than for movie
% whch may be why it is simultaneously smaller and larger than movie
% according to the kstest2

