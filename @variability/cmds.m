Variability calculations and surrogate generation:
% create nptdata object containing list of cells to process
ndcells = nptdata('SessionsFile','l5.txt');
% gather some statistics from ISI histogram
% tabluated in Refractory.xls
[ndcells,idata] = ProcessDirs(ndcells,'nptDirCmd','iobj = isi(''auto''); a(1) = struct(''type'',''.'',''subs'',''data''); a(2) = struct(''type'',''.'',''subs'',''isi''); intervals = subsref(iobj,a); edges = 0:0.2:50; hcounts = histcie(intervals,edges); data = [data; sum(hcounts)];');
[ndcells,idata] = ProcessDirs(ndcells,'nptDirCmd','iobj = isi(''auto''); a(1) = struct(''type'',''.'',''subs'',''data''); a(2) = struct(''type'',''.'',''subs'',''isi''); intervals = subsref(iobj,a); [cf,rfd,edges,hcounts,maxi,pf,pfpts] = spikeRecovery(intervals); if(~isempty(pf)) data = [data; pf(1)]; else data = [data; nan]; end');
% exclude cells with max spike counts in 50 ms lower than 50 (n=28)
% copy list from Refractory.xls and paste into variable a
s = ndcells.SessionDirs;    
s3 = {s{a}};
s4 = setdiff(s,s3);
% left with 60 cells
ndcb = nptdata('SessionDirs',s4);
[ndcb,isicb] = ProcessDirs(ndcb,'Object','isi');
plot(isicb,'SPXTickLabelsOff','LabelsOff','SPTitleOff', ...
    'XScale','linear','YScale','log','XMax',25,'Color','k');
% create refractory objects
% ProcessDirs(ndcells,'nptDirCmd','refractory(''auto'',''SaveLevels'',3);');
% ProcessDirs(ndcells,'Object','refractory','SaveLevels',3,'frameavg');
% create refractory objects with free firing rates changing at the frame rate
ProcessDirs(ndcb,'Object','refractory','Redo','SaveLevels',3,'frameavg','expfitmax',50,'AdjSpikes');
% create refractory objects using free firing rates changing at 1/10th of a movie frame
ProcessDirs(ndcb,'Object','refractory','Redo','SaveLevels',3,'FrameSplit',10,'expfitmax',50,'AdjSpikes');
InspectGUI(ndcb,'addObjs',{ndcb,ndcb}, ...
    'optArgs',{{'Objects',{{'isi',{'XScale','linear','YScale','linear','XMax',50}}}}, ...
               {'Objects',{{'refractory',{'recovery','xMax',50}}}}, ...
               {'Objects',{{'adjspikes',{'Fast'}}}}},'SP',[3 1])
% check for refresh artifacts
[ndcb,latndcb] = ProcessDirs(ndcb,'Object','latency');
InspectGUI(latndcb,'addObjs',{latndcb},'optArgs',{{},{'FFT'}},'SP',[2 1])
figure
plot(latndcb,'SPTitleOff','LabelsOff','BarColor','k','LineColor',[0.5 0.5 0.5]);
[latndcb,npit] = plot(latndcb,'FFT','ReturnVars',{'npit'});
% find the cells with no refresh artifacts
% left with 51 cells
b = find(npit{2});
ndcc = nptdata('SessionDirs',{s4{b}});
% reorder the figures with subplots so we can group the cells with
% refresh artifacts at the end
c = find(~npit{2});
d = {s4{b} s4{c}};
ndcd = nptdata('SessionDirs',d);
[ndcd,isicd] = ProcessDirs(ndcd,'Object','isi');
figure
plot(isicd,'SPXTickLabelsOff','LabelsOff','SPTitleOff', ...
    'XScale','linear','YScale','log','XMax',25,'Color','k');
[ndcd,latndcd] = ProcessDirs(ndcd,'Object','latency');
figure
plot(latndcd,'SPTitleOff','LabelsOff','BarColor',[0.5 0.5 0.5],'LineColor','k');
% reorder the figures so that cells with similar ISI histograms are grouped
% together
e = [1 16 25 34 39 40];
f = setdiff(1:60,e);
g = {d{e} d{f}};
ndce = nptdata('SessionDirs',g);
[ndce,isice] = ProcessDirs(ndce,'Object','isi');
figure
plot(isice,'SPXTickLabelsOff','LabelsOff','SPTitleOff', ...
    'XScale','linear','YScale','log','XMax',25,'Color','k');
[ndce,latndce] = ProcessDirs(ndce,'Object','latency');
figure
plot(latndce,'SPTitleOff','LabelsOff','BarColor',[0.5 0.5 0.5],'LineColor','k');

% create variability objects using the new surrogates
[ndce,va] = ProcessDirs(ndce,'Object','variability','redo', ...
    'SurrogateFF','refsga','SurrogateTRE','refsgaTRE10');

% reorder so that high refresh rates are at the beginning and low refresh
% rates are at the end so we can exclude the 9 cells with refresh artifacts
s = ndcd.SessionDirs';
s1 = {s{21:45} s{1:20} s{46:60}};
ndcf = nptdata('SessionDirs',s1);
[ndcf,latndcf] = ProcessDirs(ndcf,'Object','latency');
figure
plot(latndcf,'SPTitleOff','LabelsOff','BarColor',[0.5 0.5 0.5],'LineColor','k');
[ndcf,isicf] = ProcessDirs(ndcf,'Object','isi');
figure
plot(isicf,'SPXTickLabelsOff','LabelsOff','SPTitleOff', ...
'XScale','linear','YScale','log','XMax',25,'Color','k');
hbins = 0:0.2:25;
nhist = histcie(isicf.data.isi,hbins);
[nhm,nhmi] = max(nhist);
maxisi = hbins(nhmi)';
[(1:60)' sortrows([(1:60)' maxisi],2)]

% create firing rate fits
load froreps

% create subplot with real and surrogate data
subplot(2,1,1)
adj = adjspikes('auto');
plot(adj,'Fast','Color','k');
set(gca,'YTick',[0 50 100],'YTickLabel',[0 50 100])
set(gca,'XTick',[],'XTickLabel','')
xlabel('')
title(getDataDirs('ShortName'))
subplot(2,1,2)
sb = readSurrogateBin('refsga5.bin');
spt = cell2array(sb{50});
yi = repmat(1:100,size(spt,1),1);
plot(spt,yi,'k.')
set(gca,'YTick',[0 50 100],'YTickLabel',[0 50 100])
ylabel('Repetitions')
xlabel('Time (ms)')
title('Surrogate Spike Trains refsga5-50')

% create subplots with real data, precision, and two different surrogate
% time windows
cd ~/data/variability/
load nprecision
cd /opt/data/cat/
load ndobjs
cd ~u0509335/data
load variabilitya
global nptDataDir
nptDataDir = pwd;
InspectGUI(ndcc,'addObjs',{npc,va,va}, ...
    'optArgs',{{'Objects',{'adjspikes',{'Fast','Color','k'}}},{},{'ShowSurrogateRasters',{'refsga5.bin',50}},{'ShowSurrogateRasters',{'refsgj5.bin',50}}}, ...
    'SP',[4 1],'LinkedZoom')
InspectGUI(ndcc,'addObjs',{va,va}, ...
    'optArgs',{{'Objects',{'adjspikes',{'Fast','Color','k'}}},{'ShowSurrogateRasters',{'refsga5.bin',50}},{'ShowSurrogateRasters',{'refsgj5.bin',50}}}, ...
    'SP',[3 1],'LinkedZoom')

% create subplot with rasters and Fano Factors
cd variability
load variabilitya
InspectGUI(ndcc,'addObjs',{varefsga},'optArgs',...
    {{'Objects',{'adjspikes',{'Fast'},{}}}, ...
     {'Fano','ShowSurrogateData','XTime','SurrogateFF','refsgaFF.mat'}}, ...
     'SP',[2 1],'LinkedZoom')

% create subplot with rasters, entropy, and precision
InspectGUI(ndcc,'addObjs',{varefsga,p},'optArgs',...
    {{'Objects',{'adjspikes',{'Fast'},{}}}, ...
     {'Entropy','ShowSurrogateData','XTime','SurrogateTRE','refsgaTRE10.mat'}, ...
     {}}, ...
     'SP',[3 1],'LinkedZoom')

% create subplot with rasters and precision
cd variability
load precision
InspectGUI(ndcc,'addObjs',{npc},'optArgs',...
    {{'Objects',{'adjspikes',{'Fast'},{}}}, ...
     {}}, ...
     'SP',[2 1],'LinkedZoom')
 
% create histogram of entropy
histedges = 0:0.05:1;
n1 = histcie(varegsfa.data.entropySurrPercent,histedges);
msc1 = varefsga.data.scmean>=1;
n2 = histcie(varefsga.data.entropySurrPercent(msc1),histedges);
hbar1 = bar(histedges,n1,'histc');
hb1 = get(hbar1,'Children');
v = get(hb1,'Vertices');
v(v(:,2)==0,2) = 1;
set(hb1,'Vertices',v)
set(gca,'YScale','log')
ylim([1 1e5])
hbar2 = bar(histedges,n2,'histc');
hb2 = get(hbar2,'Children');
v2 = get(hb2,'Vertices');
v2(v2(:,2)==0,2) = 1;
set(hb2,'Vertices',v2)
xlim([0 1])
print -depsc2 -r300 -cmyk -noui EntropyHist.eps

histedgesz = -25:2:10;
z1 = histcie(varefsga.data.entropySurrZScores,histedgesz);
hzbar1 = bar(histedgesz,z1,'histc');
hz1 = get(hzbar1,'Children');
vz1 = get(hz1,'Vertices');
vz1(vz1(:,2)==0,2)=1;
set(hz1,'Vertices',vz1);
set(gca,'YScale','log')
hold on
z2 = histcie(varefsga.data.entropySurrZScores(msc1),histedgesz);
hzbar2 = bar(histedgesz,z2,'histc');
hz2 = get(hzbar2,'Children');
vz2 = get(hz2,'Vertices');
vz2(vz2(:,2)==0,2)=1;
set(hz2,'Vertices',vz2);
xlim([-23 5])
set(gca,'XTick',-21:6:5)

% create figure with rasters, surrogate rasters at 2 FrameSplit settings
% and the precision plot
InspectGUI(ndcc,'addObjs',{varefsga,varefsga,p},'optArgs',...
    {{'Objects',{'adjspikes',{'Fast'},{}}}, ...
     {'ShowSurrogateRasters',{'refsga5.bin',50}}, ...
     {'ShowSurrogateRasters',{'refsgh5.bin',50}}, ...
     {}}, ...
     'SP',[4 1],'LinkedZoom')
 
% create histogram of precision windows
% find indices that are not 0
% n0idx = p.data.trepre>0;
histbinsall = 0:8;
pnall = hist(p.data.trepre,histbinsall);
no0idx = 2:9;
pn = pnall(no0idx);
histbins = histbinsall(no0idx);
bar(histbins,pn);
figure
pie(pn,{'1','2','3','4','5','6','7','8'})
sum(pn)
% remove cell 4 and recreate pie chart
no4idx = [p.data.index(1):(p.data.index(4)-1) p.data.index(5):(p.data.index(end)-1)];
% n0idx4 = p.data.trepre(no4idx)>0;
pn4all = hist(p.data.trepre(no4idx),histbinsall);
pn4 = pn4all(no0idx);
figure
pie(pn4,{'1','2','3','4','5','6','7','8'})
sum(pn4)
% make sure everything adds up by doing histogram just for cell 4
pn4allb = hist(p.data.trepre(p.data.index(4):(p.data.index(5)-1)),histbinsall);
pn4b = pn4allb(no0idx);
pn4check = [pn4; pn4b];
pn4csum = sum(pn4check);
[pn4check; pn4csum; pn]
% another pie chart with cell 4 and all other cells separated

% see which cells have the high precision windows
p8idx = find(p.data.trepre==8);
hpn = histcie(p8idx,p.data.index,'DropLast');
figure
bar(1:51,hpn)
% rearrange according to order in ndce
% regspike = [1 16 25 34 39 40];
% regspike = [1 16 25 34 39 40];
% fastspike = setdiff(1:51,regspike);
% reorderhistidx = [regspike fastspike];
% rearrange according to order in ndcf
reorderhistidx = [21:45 1:20 46:51];
figure
hphbar = bar(1:51,hpn(reorderhistidx),1);
hph = get(hphbar,'Children');
vph = get(hph,'Vertices');
vph(vph(:,2)==0,2)=0.1;
set(hph,'Vertices',vph);
set(gca,'YScale','log')
ylim([0.1 1000])
xlim([0.5 51.5])

% check to make sure order is correct
ndces = ndce.SessionDirs';
ndce2 = nptdata('SessionDirs',{ndces{1:51}});
InspectGUI(ndce2,'Objects',{'adjspikes',{'Fast'},{}})


% write list of cells to text file refcells.txt
% in shell
cat refcells.txt | xargs ccwrapper 100 10 refsga
% once above jobs have finished, i.e. surrogates have been created
cat refcells.txt | xargs trewrapper 10 refsgaTRE10

% commands for creating figures
InspectGUI(va,'Fano','ShowSurrogateData');
InspectGUI(va,'Entropy','ShowSurrogateData')
plot(va,'MeanFano');
plot(va,'Entropy');

% code to compare surrogate generation methods 
% create surrogates for cells that were in the previous list of 32 cells
% (ndrnosync in ndobjs)
s1 = ndcc.SessionDirs';
s2 = ndrnosync.SessionDirs';
s3 = setdiff(s2,s1);
s4 = strrep(s3,[pwd filesep],'');
% save list of cells to text file refbcells.txt and call ccwrapper and
% trewrapper as above
% compare summary numbers for new and old variability objects
va = vaA0; % contains duplicate spike times
va = vaA1; % generated by Xiao Xiao, checked to have no duplicate spike times
va = vaA2; % simplified recovery function, checked to have no duplicate spike times

% these measurements are not expected to change in any of the above objects
% number of total windows
vlength = length(va.data.fano);
% number of windows with FF smaller than 1
ffidx = va.data.fano<1;
sum(ffidx)/vlength
% number of windows with mean spike count greater or equal to 1
msc1idx = va.data.scmean>=1;
msc1 = sum(msc1idx)
% number of above windows with FF smaller than 1
msc1ffidx = ffidx & msc1idx;
msc1ff = sum(msc1ffidx)
msc1ff/msc1
% get number of cells
ncells = length(va.data.framebins)
% get number of cells with msc1ff windows
cellsmsc1ffidx = unique(va.data.cellid(msc1ffidx));
ncellsmsc1ff = length(cellsmsc1ffidx)
ncellsmsc1ff/ncells
% find correlation coefficient between spike count and fano
[r,p] = corrcoef(va.data.scmean(msc1idx),va.data.fano(msc1idx))
[r,p] = spcorr(va.data.scmean(msc1idx),va.data.fano(msc1idx))
[r,p] = corrcoef(va.data.scmean(msc1ffidx),va.data.fano(msc1ffidx))
[r,p] = spcorr(va.data.scmean(msc1ffidx),va.data.fano(msc1ffidx))

% these measurements might be different for the different objects
% number of windows with FF lower than 95% of the surrogates
fspidx = va.data.fanoSurrPercent>0.95;
fsp = sum(fspidx)
% fsp for the various objects: 19, 18, 8
% find windows with spike count >= 1 and FanoSurrPercent greater than 95%
msc1fspidx = msc1idx & fspidx;
msc1fsp = sum(msc1fspidx)
% msc1fsp: 9, 9, 6
% find number of cells
cellsmsc1fspidx = unique(va.data.cellid(msc1fspidx));
ncellsmsc1fsp = length(cellsmsc1fspidx)
% ncellsmsc1fsp: 7, 7, 5
ncellsmsc1fsp/ncells
% Percent: 21.9%, 21.9%, 15.6%
% windows with TR-Entropy lower than 95% of the surrogates
nspidx = va.data.entropySurrPercent>0.95;
nsp = sum(nspidx)
% nsp: 2689, 2689, 2674
% windows with spike count >= 1 and TR-Entropy lower than 95% of the surr
msc1nspidx = msc1idx & nspidx;
msc1nsp = sum(msc1nspidx)
% msc1nsp: 236, 236, 234
% cells with above conditions satisfied
cellsmsc1nspidx = unique(va.data.cellid(msc1nspidx));
ncellsmsc1nsp = length(cellsmsc1nspidx)
% ncellsmsc1nsp: 25, 25, 24
ncellsmsc1nsp/ncells
% Percent: 78.1%, 78.1%, 75%
% check correlation coefficient between FanoSurrZScores and
% entropySurrZScores
[r,p] = corrcoef(va.data.fanoSurrZScores(msc1idx), ...
    va.data.entropySurrZScores(msc1idx))
% R (p): -0.0949 (0.0307), -0.1287 (0.0033), -0.1025 (0.0196)
[r,p] = corrcoef(va.data.fanoSurrZScores(msc1ffidx), ...
    va.data.entropySurrZScores(msc1ffidx))
% R (p): 0.4172 (0.0001), 0.3656 (0.0006), 0.4312 (0)

% windows with Z-scores smaller than 2 
% use 2 std because 95% of the data will be between -2 and +2 std
nsp2idx = va.data.entropySurrZScores<-2;
nsp2 = sum(nsp2idx)
% windows with spike count >= 1 and TR-Zscore smaller than -2
msc1nsp2idx = msc1idx & nsp2idx;
msc1nsp2 = sum(msc1nsp2idx)
% cells with above conditions satisfied
cellsmsc1nsp2idx = unique(va.data.cellid(msc1nsp2idx));
ncellsmsc1nsp2 = length(cellsmsc1nsp2idx)
ncellsmsc1nsp2/ncells

% in shell
find $(cat l5.txt) -name "refractory.mat" | xargs tar -cvzf rf3.tgz
scp rf3.tgz ogier:~/home/cat/newcatdata3
% on ogier
cd home/cat/newcatdata3
gzip -dc rf3.tgz | tar -xvf -
% submit pbs jobs
find . -type d -path "*cluster*" | xargs ~/bin/ccwrapper 100 10 framesg
find . -type d -path "*cluster*" | xargs ~/bin/trewrapper 10 framesgTRE10

% check for duplicate spike times
[ndcc,zisi] = ProcessDirs(ndcc,'nptDirCmd','for sfi = 1:10,spt2 = readSurrogateBin([''refsga'' num2str(sfi) ''.bin'']); for sni = 1:100, zi = find(diff(cell2array(spt2{sni}))==0); data = [data; [repmat([sfi sni],size(zi,1),1) zi]]; end; end;');

[ndresp,sgisi] = ProcessDirs(ndresp,'nptDirCmd','a1 = readSurrogateBin(''framesg1.bin''); sgn = length(a1); repn = length(a1{1}); repsd = zeros(repn,1); for repi = 1:repn, repd = zeros(sgn,1); for sgi = 1:sgn, repd(sgi) = ~isempty(find(diff(a1{sgi}{repi})==0)); end; repsd(repi) = sum(repd); end; isireps = find(repsd); data = [data; [repmat(i,length(isireps),1) isireps repsd(isireps)]];');

% explanation of above command
% read surrogate file
a1 = readSurrogateBin('framesg1.bin'); 
% get number of surrogates in each file
sgn = length(a1); 
% get number of repetitions
repn = length(a1{1}); 
for fi = 1:10
	% allocate memory
	repsd = zeros(repn,1); 
	% loop over repetitions
	for repi = 1:repn, 
		% allocate memory to store 1 if surrogate contains 0 isi
		repd = zeros(sgn,1); 
		% loop over surrogates
		for sgi = 1:sgn, 
			% find inter-spike-intervals that are 0
			% set repd(sgi) to 1 if there are 0 isi
			repd(sgi) = ~isempty(find(diff(a1{sgi}{repi})==0)); 
		end; % end loop over number of surrogates 
		% find number of surrogates with 0 isi for this repetition
		repsd(repi) = sum(repd); 
	end; % end loop over number of repetitions
	% find repetitions that have isi's that are 0
	isireps = find(repsd); 
	% return cell number, repetition number, and number of isi's that are 0
	data = [data; [repmat(i,length(isireps),1) isireps repsd(isireps)]];
	if(fi<10)
		a1 = readSurrogateBin(['framesg' num2str(fi+1) '.bin']); 
	end
end
% find unique cell numbers
us = unique(sgisi(:,1));
% get session directories
s = ndresp.SessionDirs';
% find unique session directories
usd = s(us);
% remove path prefix
usd2 = strrep(usd,'/Users/syen/Documents/ShihCheng/Data/Neural/Cat/newcatdata/','');
