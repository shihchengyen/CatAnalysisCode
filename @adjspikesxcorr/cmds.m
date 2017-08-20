[ndintra,axcintra] = ProcessDirs(ndintra,'Object','adjspikesxcorr','redo');
[ndinter,axcinter] = ProcessDirs(ndinter,'Object','adjspikesxcorr','redo');
InspectGUI(axcintra,'NoSig')
InspectGUI(axcintra,'addObjs',{ngintra},'optArgs', ...
	{{'NoSig'}, ...
	 {'GroupEvent','Object', ...
	 	{'adjspikesxcorr',{'SeparateGroupPlot','NoSig'},{'AutoCorr'}}}}, ...
	'SP',[2 1])
InspectGUI(axcintra,'addObjs',{ngintra,ngintra,ngintra},'optArgs', ...
	{{'NoSig'}, ...
	{'GroupEvent','Object', ...
		{'adjspikesxcorr',{'SeparateGroupPlot','NoSig'},{'AutoCorr'}}}, ...
	{'GroupEvent','Object', ...
		{'isi',{'XScale','linear'},{'AdjSpikes','redo','NoBurst'}},'GroupPlotSep','Horizontal'}, ...
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
