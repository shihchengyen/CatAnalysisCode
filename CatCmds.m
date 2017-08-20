Variability calculations and surrogate generation:
	% create nptdata object containing list of cells to process
	ndcells = nptdata('SessionsFile','l5.txt');
	% create refractory objects
	% ProcessDirs(ndcells,'nptDirCmd','refractory(''auto'',''SaveLevels'',3);');
	ProcessDirs(ndcells,'Object','refractory','SaveLevels',3,'frameavg');

	% in shell
	find $(cat l5.txt) -name "refractory.mat" | xargs tar -cvzf rf3.tgz
	scp rf3.tgz ogier:~/home/cat/newcatdata3
	% on ogier
	cd home/cat/newcatdata3
	gzip -dc rf3.tgz | tar -xvf -
	% submit pbs jobs
	find . -type d -path "*cluster*" | xargs ~/bin/ccwrapper 100 10 framesg
	find . -type d -path "*cluster*" | xargs ~/bin/trewrapper 10 framesgTRE10
	
s = textread('movie-cells.txt','%s');
s1 = strcat([pwd filesep],s);
ndcells = nptdata('SessionDirs',s1);
ngsite = ProcessDays(nptgroup,'IntraGroup','Cells',ndcells.SessionDirs, ...
    'AnalysisLevel','AllIntraGroup');
InspectGUI(ndcells,'Objects',{{'adjspikes',{'Fast'}}})
InspectGUI(ngsite2,'Object',{'adjspikes',{'Fast'}},'GroupEvent', ...
    'GroupPlotSep','Vertical')

ndintra = ProcessDays(nptdata,'AnalysisLevel','pairs','Intragroup','Cells',textread('movie-responsive.txt','%s'));
ndinter = ProcessDays(nptdata,'AnalysisLevel','pairs','Intergroup','Cells',textread('movie-responsive.txt','%s'));

ProcessDirs(ndintra,'nptDirCmd','adjspikesxcorr(''auto'',''Sliding'',''redo'',''NumSurrogates'',1000,''OverlapIncrBins'',''frame/2'',''JobWindows'',100,''DataFile'',''asxcdata1.mat'',''SurrFilePrefix'',''asxcsurrA'');');
	
InspectGUI(ndinter,'Objects', ...
	{{'nptgroup',{'Object',{'adjspikes',{'Fast','Interleaved'}}}, ...
            {'GetClusterDirs'}}; ...
	{'nptgroup',{'Object',{'firingrate',{'PSTH'},{'Repetitions','Rate','redo'}}}, ...
		{'GetClusterDirs'}}; ...
	{'adjspikesxcorr',{'SurrHist','Sliding'}}},'SP',[3 1],'LinkedZoom')
InspectGUI(ndintra,'Objects', ...
	{{'nptgroup',{'Object',{'adjspikes',{'Fast','Interleaved'}}}, ...
            {'GetClusterDirs'}}; ...
	{'nptgroup',{'Object',{'firingrate',{'PSTH'},{'Repetitions','Rate','redo'}}}, ...
		{'GetClusterDirs'}}},'SP',[2 1],'LinkedZoom')
InspectGUI(ndinter,'Objects', ...
	{{'adjspikesxcorr',{'SurrHist','Sliding'}}; ...
    {'adjspikesxcorr',{'ZScore','Sliding'}}; ...
	{'nptgroup',{'Object',{'latency'}},{'GetClusterDirs'}}},'SP',[3 1],'LinkedZoom')
InspectGUI(ndintra,'addObjs',{axcintra3,axcintra3,axcintraA3,axcintraA3}, ...
    'optArgs',{ {'Objects',{{'nptgroup',{'Object',{'adjspikes',{'Fast','Interleaved'}}}, ...
                    {'GetClusterDirs'}}}}; ...
                {'SurrHist','Sliding'}; ...
                {'SurrMeanStd','Sliding'}; ...
                {'SurrHist','Sliding'}; ...
                {'SurrMeanStd','Sliding'}},'SP',[5 1],'LinkedZoom')
InspectGUI(axcintra,'addObjs',{ngintra},'optArgs', ...
	{{'NoSig'}, ...
	{'GroupEvent','Object',{'adjspikesxcorr',{'SeparateGroupPlot','NoSig'},{'AutoCorr'}}}}, ...
	'SP',[2 1])
InspectGUI(ngintra,'addObjs',{axcintra},'optArgs', ...
	{{'GroupEvent','Object',{'adjspikes',{'Fast','Interleaved'},{'GetClusterDirs'}}},{'NoSig'}}, ...
	'SP',[2 1])
InspectGUI(axcintra,'addObjs',{ngintra,ngintra,ngintra},'optArgs', ...
	{{'NoSig'}, ...
	{'GroupEvent','Object',{'adjspikesxcorr',{'SeparateGroupPlot','NoSig'}}}, ...
	{'GroupEvent','Object',{'isi',{'XScale','linear'},{'AdjSpikes','redo','NoBurst'}},'GroupPlotSep','Horizontal'}, ...
	{'Object',{'latency'},'GroupEvent'}}, ...
	'SP',[4 1])
InspectGUI(axcinter,'addObjs',{nginter,nginter,nginter},'optArgs', ...
	{{'NoSig'}, ...
	{'GroupEvent','Object',{'adjspikesxcorr',{'SeparateGroupPlot','NoSig'}}}, ...
	{'GroupEvent','Object',{'isi',{'XScale','linear'},{'AdjSpikes','redo','NoBurst'}},'GroupPlotSep','Horizontal'}, ...
	{'Object',{'latency'},'GroupEvent'}}, ...
	'SP',[4 1])
InspectGUI(axcintra,'addObjs',{ngintra,ngintra},'optArgs', ...
	{{'NoSig'}, ...
	{'GroupEvent','Object',{'adjspikesxcorr',{'SeparateGroupPlot','NoSig'}}}, ...
	{'Object',{'latency'},'GroupEvent'}}, ...
	'SP',[3 1])
rp = 75;
InspectGUI(ndintra,'addObjs',{frdintra,frdintra},'optArgs', ...
    {{'Objects',{{'nptgroup',{'Object',{'adjspikes',{'Fast','Interleaved'}}},{'GetClusterDirs'}}}}, ...
	 {'FRMean','RatePercentile',rp}, ...
     {'WinNoiseCorr','RatePercentile',rp}},'SP',[3 1],'LinkedZoom',1)
InspectGUI(ndinter,'addObjs',{frdinter,frdinter},'optArgs', ...
    {{'Objects',{{'nptgroup',{'Object',{'adjspikes',{'Fast','Interleaved'}}},{'GetClusterDirs'}}}}, ...
	 {'FRMean','RatePercentile',rp}, ...
     {'WinNoiseCorr','RatePercentile',rp}},'SP',[3 1],'LinkedZoom',1)
InspectGUI(frdintra,'addObjs',{frdintra,frdintra},'optArgs',{{'FRMean'},{'WinNoiseCorr'},{'HistWinFRNoiseCorr','RatePercentile',rp,'HistBins',rp:5:100,'HistBinsY',-1:0.1:1}},'SP',[3 1]);
InspectGUI(frdintra,'WinFRNoiseScatter','RatePercentile',rp);
InspectGUI(ndintra,'addObjs',{frdintra,frdintra,frdintra},'optArgs', ...
    {{'Objects',{{'nptgroup',{'Object',{'adjspikes',{'Fast','Interleaved'}}},{'GetClusterDirs'}}}}, ...
	 {'FRMean','RatePercentile',rp}, ...
	 {'CV','RatePercentile',rp}, ...
     {'WinNoiseCorr','RatePercentile',rp}},'SP',[4 1],'LinkedZoom',1)
% create nptgroup object to plot spike trains of all responsive cells in a group
ngsite = ProcessDays(nptgroup,'IntraGroup','Cells',ndresp.SessionDirs,'AnalysisLevel','AllIntraGroup');
ngsite2 = ProcessDays(nptgroup,'IntraGroup','Cells',ndcells.SessionDirs,'AnalysisLevel','AllIntraGroup');
% create nptgroup object to plot all pairs of responsive cells in a group
ngsitepairs = ProcessDays(nptgroup,'IntraGroup','Cells',ndresp.SessionDirs,'AnalysisLevel','AllPairs');
ngsitepairs2 = ProcessDays(nptgroup,'IntraGroup','Cells',ndcells.SessionDirs,'AnalysisLevel','AllPairs');
InspectGUI(ngsite,'addObjs',{ngsitepairs},'optArgs',{{'Object',{'adjspikes',{'Fast','Interleaved'}},'GroupEvent'}, ...
        {'Object',{'frdiff',{'WinNoiseCorr'}},'GroupEvent','GroupPlotSep','Vertical'}},'SP',[2 1])