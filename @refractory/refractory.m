function obj = refractory(varargin)
%refractory Object used to compute spike recovery function
%   OBJ = refractory('auto',varargin) attempts to load a mat file with
%   name '*g*c*.mat' to create the REFRACTORY object.
%
%   OBJ = refractory(CELLSTRUCT,'CellName',CNAME,varargin) computes the 
%   spike recovery function using the data structure CELLSTRUCT and uses
%   CNAME as the name of the cell since the name of the cell is not
%   stored in CELLSTRUCT.
%
%   The optional input arguments are:
%      isibinsize - bin size of the histogram in ms (default 0.2).
%
%      pfpts - number of points after the peak to use to compute q
%              (default 4).
%
%      expfitmax - maximum interval used to find peak of histogram
%                  in ms (default 15).
%
%      qtbinsize - bin size of free firing rate as well as spike
%                  probability function and PSTH in ms (default 0.2).
%
%      duration - duration of repetition in ms (default 30000).
%
%      binsize - changes isibinsize and qtbinsize at the same time.
%
%   The object fields are:
%      OBJ.data.isibinsize - step size of histogram (default is 0.2 ms)
%      OBJ.data.expfitmax - max time in ms to look for max of histogram
%                      (default is 15 ms).
%      OBJ.data.pfpts - number of points past the peak to use to fit 
%                  exponential function (default is 4).
%      OBJ.data.qtbinsize() - bin size of q(t), W(t) & r(t) in ms.
%      OBJ.data.duration - duration of repetition.
%      OBJ.data.edges() - ISI histogram bin values (one vector for all 
%                       cells).
%      OBJ.data.nCells - number of cells in object.
%      OBJ.data.nEdges - maximum number of bins in the ISI histogram.
%      OBJ.data.rfd{} - recovery function data for each cell.
%      OBJ.data.hcounts{} - ISI histogram counts for each cell.
%      OBJ.data.bins() - number of histogram bins.
%      OBJ.data.maxi() - index corresponding to max of ISI histogram for 
%                      each cell.
%      OBJ.data.pf() - array containing polyfit values of the exponential 
%               decay function in log units. pf(1) is the slope
%               and pf(2) is the intercept.
%      OBJ.data.qt{} - free firing rate for each cell.
%      OBJ.data.Wt{} - spike probability functionn for each cell.
%      OBJ.data.rt{} - PSTH for each cell.
%      OBJ.data.rtEdges{} - bin edges used for the PSTH. Useful for using 
%                       the function bar(,rtEdges,'histc').
%      OBJ.data.repetitions() - the number of stimulus repetitions for
%                             each cell.
%      OBJ.data.wt{} - recovery function for each cell.
%
%   OBJ = REFRACTORY(jonscells,'isibinsize',0.2,'pfpts',4,'expfitmax',15,)

% default values for optional arguments
Args = struct('isibinsize',0.2, ... % bin size of ISI histogram
			  'expfitmax',15, ... % max of exponential fit in ms
			  'pfpts',4, ... % number of points to use for polyfit
			  'qtbinsize',0.2, ... % bin size of qt, Wt & rt
			  'duration',30000, ... % duration of repetition
			  'CellName','', ...
			  'RedoLevels',0,'SaveLevels',0, ...
			  'Auto',0);
			  
[Args,varargin2] = getOptArgs(varargin,Args, ...
	'aliases',{'binsize' {'isibinsize','qtbinsize'}}, ...
	'flags',{'Auto'}, ...
	'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
	'subtract',{'RedoLevels','SaveLevels'});

if nargin==0
	% create empty object
	obj = createEmtpyObject;
elseif( (nargin==1) & (isa(varargin{1},'refractory')) )
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
		% check for saved object
		if(ispresent('refractory.mat','file','CaseInsensitive') & (Args.RedoLevels==0))
			fprintf('Loading saved refractory object...\n');
			% load saved object and exit
			l = load('refractory.mat');
			obj = l.rf;
			return
		else
			% try instantiating isi object
			ispisi = isi('auto',varargin2{:});
			if(~isempty(ispisi))
				% get number of cells
				nCells = 1;
				n = nptdata(nCells,0,pwd);
				data.cellname = pwd;
				data.nCells = nCells;
				data.nEdges = 0;
				% create fields so the field order in the structure is the same
				% use cell arrays since these vectors might not be the same
				% length if the refresh rate of the stimulus was different
				data.isibinsize = Args.isibinsize;
				data.pfpts = [];
				data.expfitmax = Args.expfitmax;
				data.duration = Args.duration;
				data.qtbinsize = Args.qtbinsize;
				data.rfd = {};
				data.hcounts = {};
				data.maxi = 0;
				data.pf = 0;
				data.bins = 0;
				data.qt = {};
				data.Wt = {};
				data.rt = {};
				data.rtEdges = {};
				data.repetitions = 0;
				data.wt = {};
	
				% try instantiating stiminfo object
% 				[cf,data.rfd{1},edges,data.hcounts{1},data.maxi, ...
% 					data.pf,data.pfpts] = spikeRecovery(ispisi.data.isi,varargin2{:});
				[cf,data.rfd{1},edges,data.hcounts{1},data.maxi] ...
                    = spikeRecovery(ispisi.data.isi,varargin2{:});
				% check if edges is longer than data.nEdges
				lEdges = length(edges);
				data.bins = lEdges;
				data.edges = edges;
				data.nEdges = lEdges;
                
                % check if there were enough spikes to do a proper fit
                if(isempty(cf))
                    % return empty refractory object
                    d.data = data;
                    obj = class(d,'refractory',n);
                    return
                end
                
				% instantiate stiminfo object
				st = stiminfo('auto',varargin2{:});
				% instantiate adjspikes object
				adjs = adjspikes('auto',varargin2{:});
				jcells.stimulus_info.adjusted_frames_vector = adjs.data.adjFramePoints;
				jcells.stimulus_info.start_frame = st.data.iniInfo.start_frame;
				jcells.stimulus_info.end_frame = st.data.iniInfo.end_frame;
				jcells.cell_info.adjusted_spiketrain = adjs.data.adjSpiketrain;
				
				% get free firing rate, spike probability function, and PSTH
				[data.qt{1},data.Wt{1},data.rt{1},data.rtEdges{1}, ...
					data.repetitions,data.qtbinsize,data.wt{1},data.qtframebins] = ...
					getFreeFiringRate2(jcells,data.rfd{1},'curvefitobj',cf, ...
					'ISIBinSize',Args.isibinsize,varargin2{:});
				d.data = data;
				obj = class(d,'refractory',n);
				if(Args.SaveLevels>0)
					fprintf('Saving refractory object...\n');
					rf = obj;
					save refractory rf
				end
                return
			else
				% no saved object so try to create one			
				% try to find name of .mat file
				a = nptDir('*g*c*.mat');
				if(isempty(a))
					% create empty object and exit from function
					obj = createEmptyObject;
					return
				else
					mat = load(a(1).name);
					jcells = mat.data;
					% try to set the name of the cell
					[fpath,fname] = nptFileParts(a(1).name);
					data.cellname = fname;
				end
			end
		end
	else
		% hopefully we were given a data structure
		% create object using data structure
		jcells = varargin{1};
		% set name to argument if present
		data.cellname = Args.CellName;
	end
	% figure which version of the structure we got
	if(isfield(jcells,'cell_info'))
		% this is version 2
		structversion = 2;
	else
		structversion = 1;
	end
	% get number of cells
	nCells = length(jcells);
	n = nptdata(nCells,0);
	data.nCells = nCells;
	data.nEdges = 0;
	% create fields so the field order in the structure is the same
	% use cell arrays since these vectors might not be the same
	% length if the refresh rate of the stimulus was different
	data.isibinsize = Args.isibinsize;
	data.pfpts = [];
	data.expfitmax = Args.expfitmax;
	data.duration = Args.duration;
	data.qtbinsize = repmat(Args.qtbinsize,nCells,1);
	data.rfd = cell(1,nCells);
	data.hcounts = cell(1,nCells);
	data.maxi = zeros(nCells,1);
	data.pf = zeros(nCells,2);
	data.bins = zeros(nCells,1);
	data.qt = cell(1,nCells);
	data.Wt = cell(1,nCells);
	data.rt = cell(1,nCells);
	data.rtEdges = cell(1,nCells);
	data.repetitions = zeros(nCells,1);
	data.wt = cell(1,nCells);

	fprintf('Calculating cell ');
	if(structversion==2)
		for i = 1:nCells
			fprintf('%i ',i);
			[cf,data.rfd{i},edges,data.hcounts{i},data.maxi(i), ...
				data.pf(i,1:2),data.pfpts] = spikeRecovery(diff(jcells(i).cell_info.original_spiketrain), ...
				varargin{:});
			% check if edges is longer than data.nEdges
			lEdges = length(edges);
			data.bins(i) = lEdges;
			if (lEdges>data.nEdges)
				data.edges = edges;
				data.nEdges = lEdges;
			end
			% get free firing rate, spike probability function, and PSTH
			[data.qt{i},data.Wt{i},data.rt{i},data.rtEdges{i}, ...
				data.repetitions(i),data.qtbinsize(i),data.wt{i},data.qtframebins] = ...
				getFreeFiringRate2(jcells(i),'curvefitobj',cf, ...
				varargin{:});
		end
	elseif(structversion==1)
		for i = 1:nCells
			fprintf('%i ',i);
			[cf,data.rfd{i},edges,data.hcounts{i},data.maxi(i), ...
				data.pf(i,1:2),data.pfpts] = spikeRecovery(jcells(i).isi,varargin{:});
			% check if edges is longer than data.nEdges
			lEdges = length(edges);
			data.bins(i) = lEdges;
			if (lEdges>data.nEdges)
				data.edges = edges;
				data.nEdges = lEdges;
			end
			% get values for recovery function
			data.wt{i} = getFunctionValues(cf,'binsize',data.qtbinsize,varargin{:});
			% get free firing rate, spike probability function, and PSTH
			[data.qt{i},data.Wt{i},data.rt{i},data.rtEdges{i}, ...
				data.repetitions(i)] = getFreeFiringRate(jcells(i).raster,...
											data.wt{i},varargin{:});
		end
	end
	fprintf('\n');
	d.data = data;
	obj = class(d,'refractory',n);
	if(Args.SaveLevels>0)
		fprintf('Saving refractory object...\n');
		rf = obj;
		save refractory rf
	end
end

function obj = createEmtpyObject

n = nptdata(0,1);
data.nCells = 0;
data.nEdges = 0;
d.data = data;
obj = class(d,'refractory',n);
