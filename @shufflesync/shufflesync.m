function obj = shufflesync(varargin)
%@shufflesync Constructor function for shufflesync class
%   OBJ = shufflesync('auto') attempts to create a shufflesync 
%   object by ...
% 
%   ClusterDirs - Cell array containing the two cluster directories that
%                 are to be used in this calculation. If not specified, 
%                 getDataDirs is called and it will try to determine the
%                 two cluster directories using the current directory name,
%                 assuming the current directory is a combination 
%                 directory (default: {}).
%   BinSize - Size of bins in ms used to bin spikes. All calculations are 
%             performed in indices of these bins (default: 1).
%   WindowBins - Size of each window in number of bins or in terms of 
%                'frame' (default: 'frame').
%   WindowIncrBins - Size of window increment in number of bins or in terms 
%                    of 'frame' (default: 'frame/2').
%   Shuffle - Width of uniform distribution in ms used to jitter the spikes
%             (default: 2).
%   NumSurrogates - Number of surrogates to be generated (default: 1000).
%   DataFile - Name of data file (default: 'shsydata.mat'). This file 
%              contains variables computed entirely from the data so it can
%              be shared between runs with different values for Shuffle.
%              The creation of this file will be skipped if it is already
%              present, unless the RedoDataFile is specified (see below).
%   SurrData' - Name of surrogate data file (default: 'shsysdata.mat'). 
%               This file contains variables specific to the value for 
%               Shuffle. This file is always created, even if it is already
%               present. 
%   SurrFilePrefix - Suffix of surrogate files that are created (default: 
%                    'shsysurr').
%   ResultFile - Name of result file (default: 'shsyresult.mat').
%   RedoDataFile - Flag that specifies that the data file should be 
%                  generated even if it already exists.
%   RecalcData - Flag that returns an empty object after calculating the 
%                data files. 
%   JobWindows - Number of windows in each job submitted to the queue
%                (default: 20).
%   QsubProgram - Name of shell script used to submit jobs to the queue
%                 (default: 'shsywrapper').
%   NoQsub - Flag that indicates that all computations should be performed
%            on the local machine. 
%   WaitForQsub - Flag that indicates that this function should wait till 
%                 the result file is created. 
%   QsubWaitTime - Time in ms to wait between checks to see if the result 
%                  file has been created (default: 10).
%   SkipNoResult - Flag that creates an empty object if the result file is 
%                  not found.
%   SurrHistMin - Minimum value for the histogram of spike pairs at each 
%                 time lag (default: 0). This value is saved to the 
%                 surrogate data file and used in the computation of the 
%                 result file.
%   SurrHistStep - Step size in ms of histogram of spike pairs at each time 
%                  lag (default: 1). This value is saved to the surrogate 
%                  data file and used in the computation of the result 
%                  file.
%   NumCentralBins - Number of central bins that considered synchronous
%                    (default: 3). This value is saved to the surrogate 
%                    data file and used in the computation of the result 
%                    file. 
%   PSTHWindow - Number of bins over which the firing rate is averaged to 
%                compute the smoothed PSTH (default: 11). This value is
%                saved to the surrogate data file and is used in the
%                computation of the result file. 

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0,'ClusterDirs',{''}, ...
	'BinSize',1,'WindowBins','frame','WindowIncrBins','frame/2','Shuffle',2, ...
    'NoQsub',0,'WaitForQsub',0,'QsubWaitTime',10,'JobWindows',20, ...
    'QsubProgram','shsywrapper', 'NumSurrogates',1000, ...
    'SurrFilePrefix','shsysurr','DataFile','shsydata.mat', ...
    'ResultFile','shsyresult.mat','SkipNoResult',0, ...
    'SurrHistMin',0,'SurrHistStep',1,'NumCentralBins',3,'RecalcData',0, ...
    'SurrData','shsysdata.mat','RedoDataFile',0,'PSTHWindow',11);
Args.flags = {'Auto','WaitForQsub','NoQsub','SkipNoResult','RecalcData', ...
	'RedoDataFile'};
[Args,varargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'shufflesync';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'shsy';

if(nargin==0)
	% create empty object
	obj = createEmptyObject(Args);
elseif( (nargin==1) & isa(varargin{1},Args.classname))
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
		% check for saved object
		if(~isempty(nptDir(Args.matname,'CaseInsensitive')) ...
			&& (Args.RedoLevels==0))
			fprintf('Loading saved %s object...\n',Args.classname);
			l = load(Args.matname);
			obj = eval(['l.' Args.matvarname]);
		else
			% no saved object so we will try to create one
			% pass varargin in case createObject needs to instantiate
			% other objects that take optional input arguments
			obj = createObject(Args,varargin{:});
		end
	end
end

function obj = createObject(Args,varargin)

% check to see if the result file is present
resfile = nptDir(Args.ResultFile);
if(Args.RecalcData || isempty(resfile))
	if(Args.SkipNoResult)
		% create empty object
		obj = createEmptyObject(Args);
	else
		% create object
		% check if ClusterDirs is empty
		if(isempty(Args.ClusterDirs))
			% try finding ClusterDirs using current directory
			Args.ClusterDirs = getDataDirs('GetClusterDirs');
		end
		
		if(~isempty(Args.ClusterDirs))
			% check that there are only 2 directories
			if(size(Args.ClusterDirs,2)>2)
				fprintf('Warning: %s works only on cell pairs! Using the first two cells...\n', ...
					Args.classname);
			end
			% get current directory
			cwd = pwd;
			% change to first directory
			cd(Args.ClusterDirs{1})
			% instantiate the adjspikes object instead of the ispikes object 
			% since we need to be able shift the repetitions with respect to 
			% each other to compute the significance values for the xcorr
			% pass varargin so things like RedoLevels and SaveLevels will work
			adjs1 = adjspikes('auto',varargin{:});
			% return to previous directory in case directory paths are relative
			cd(cwd);
			% change to second directory
			cd(Args.ClusterDirs{2});
			adjs2 = adjspikes('auto',varargin{:});
			% return to previous directory even if any of the objects are empty
			cd(cwd);
			% check if any of the required objects are empty
			if( isempty(adjs1) | isempty(adjs2) )
				% create empty object
				obj = createEmptyObject(Args);
			else % if( isempty(adjs1) | isempty(adjs2) )
				% set some variables which will be used whether we are generating 
				% surrogates or not
				% get number of repetitions
				reps = size(adjs1.data.raster,2);
				% get frame points
				framepoints = adjs1.data.adjFramePoints;
				% get total number of frames over all repetitions
				frames = length(framepoints) - 1;
				% get number of frames in a repetition
				repframes = frames/reps;
				% first get frame points from the first repetition
				% need to subtract framepoints(1) since the first
				% framepoint could be non-zero
				fp2 = framepoints(1:(repframes+1)) - framepoints(1);
				% divide frame points by number of surrogate windows
				[subbins,binSize,nsbins] = divideBins(fp2,'SubBinSize', ...
					Args.BinSize);
				winbins = Args.WindowBins;
				nsbinstr = num2str(nsbins);
				if(ischar(winbins))
					% replace frame with nsbins
					sw2 = strrep(winbins,'frame',nsbinstr);
					% evaluate the string to get subwindow size and call ceil since
					% subwindow size has to be an integer and we want to make
					% sure we cover the entire frame
					windowlength = ceil(eval(sw2));
				else
					windowlength = ceil(winbins);
				end
				windowdiff = windowlength - 1;

				% check if the data file exists
				if(isempty(dir(Args.DataFile)) || Args.RedoDataFile)
					% get the bins corresponding to the beginning of a frame
					framestartbins = (0:(repframes-1))*nsbins + 1;
					overlap = Args.WindowIncrBins;
					if(ischar(overlap))
						% replace frame with nsbins
						ov2 = strrep(overlap,'frame',nsbinstr);
						% evaluate the string to get overlap
						WindowIncrBins = round(eval(ov2));
					else
						WindowIncrBins = round(overlap);
					end
					% get the start bins for the overlaps
					overlapbins = 0:WindowIncrBins:(nsbins-1);
					% get the start bins for all windows
					startbins0 = [ones(length(overlapbins),1) overlapbins'] ...
						* [framestartbins; ones(1,length(framestartbins))];
					startbins = startbins0(:);
					% calculate total number of windows
					numwindows = length(startbins);
					% extend startbins by window length so as to eliminate 
					% boundary effects in the xcorr calculation
					% replace negative values with 1 since binidx will 
                    % contain 0's for all values not in the range of bin
                    % values. This will cause all the 0's to be found to
                    % belong to windows with negative bin limits
					startbins2 = max([startbins - windowlength ones(numwindows,1)],[],2);
			
					% need to generate surrogates
					% use raster field since spike times are already separated
					% into repetitions so we can access the indices for each
					% repetition easily
					% concatenate raster field in both adjspikes objects
					raster12 = concatenate(adjs1.data.raster,adjs2.data.raster, ...
						'Columnwise');
					% take histogram			
					[binsc,binidx] = histcie(raster12,subbins,'DropLast','DataCols');
					% get size of binidx
					binidxsize = size(binidx);
					% get the end bins
					endbins = startbins + windowdiff;
					% extend startbins by window length so as to eliminate 
					% boundary effects in the xcorr calculation
					endbins2 = endbins + windowlength;
					cellvec = [1 2];
					needrefcheck = zeros(2,1);
					spiketrains = cell(2,1);
					% save start times for each window which will be used by 
					% the makeObject function
					winstarttime = subbins(startbins);
					% create variable used when no surrogates are generated
					% because there were no spikes or because there were no
					% reps with spikes in both cells
					% psthtemp = sparse(zeros(lagsmax+1,dnssize));
					% convert binsc and binidx to sparse matrices to save 
					% disk space
					binsc = sparse(binsc);
					binidx = sparse(binidx);
					% create matrix to store limits
					tmpstartlimits = ones(reps,4);
                    % don't use sparse matrices since they don't work
                    % properly when working with empty matrices (e.g. 
                    % endlimits(repnums1,CEN1COL) = cell1limits does not
                    % work if repnums1 and cell1limits are both empty. The
                    % above code works fine with full matrices.
					tmpendlimits = zeros(reps,4);
											
					% save variables so the function that generates surrogates
					% won't have to recompute them over and over again
					% save numwindows and winstarttime which will be used by 
					% the makeObject function
					save(Args.DataFile,'binsc','binidx','binidxsize','reps', ...
						'windowlength','startbins','endbins','startbins2', ...
						'endbins2','cellvec','needrefcheck','spiketrains', ...
						'numwindows','winstarttime','subbins','tmpstartlimits', ...
						'tmpendlimits','subbins','binSize');
                else
                    % load data file since some of the variables (e.g. numwindows) 
                    % are needed later
                    load(Args.DataFile);
				end % if(isempty(dir(Args.DataFile)) || Args.RedoDataFile)

				shuffle = Args.Shuffle;
				% length of bins used to compute xcorr
				lagsmax = windowdiff + (2*shuffle);
				lagslength = (2*lagsmax) + 1;
				numSurrogates = Args.NumSurrogates;
				surrvec = 1:numSurrogates;
				% first column of xc matrix in output of asxcSurrSynchrony is
				% the xcorr of the actual data followed by the xcorr of the
				% shift predictor and then followed by NumSurrogates
				% columns of surrogate data. So we need to add 2 to
				% numSurrogates and surrvec.
				numSurrogates2 = numSurrogates + 2;
				% first column of the sptrains matrix is the data followed
				% by numSurrogates columns of surrogate spike trains. So
				% we need to add 1 to numSurrogates and surrvec.
				surrvec2 = surrvec + 2;
				% matrix to store xcorr for data and surrogates, e.g. size is (35*2+1)x1001
				numSurrogates1 = numSurrogates + 1;
				surrvec1 = surrvec + 1;
				xctemp = zeros(lagslength,numSurrogates2);
				shufflevecl = (2 * shuffle) + 1;
				shufflevec = [ones(1,shufflevecl); -shuffle:shuffle];
				% bins used to compute xcorr
				lagbins = -lagsmax:lagsmax;
                % create lagbins2 which will allow us to ignore values
                % outside the range of interest using histcie
                lagbins2 = [lagbins-0.5 lagsmax+0.5];
				% set variables so we can save them easily
				cellpairdatafile = Args.DataFile;
				surrprefix = Args.SurrFilePrefix;
				resultfile = Args.ResultFile;
				NumCentralBins = Args.NumCentralBins;
				SurrHistMin = Args.SurrHistMin;
				SurrHistStep = Args.SurrHistStep;
				psthwin = Args.PSTHWindow;
				
				save(Args.SurrData,'shuffle','shufflevecl','shufflevec', ...
					'lagbins','lagbins2','surrvec','numSurrogates','xctemp', ...
					'numSurrogates1','numSurrogates2','surrvec1', ...
					'surrvec2','cellpairdatafile','surrprefix','resultfile', ...
					'NumCentralBins','SurrHistMin','SurrHistStep', ...
					'lagslength','psthwin');
					
                if(Args.RecalcData)
                    obj = createEmptyObject(Args);
                elseif(Args.NoQsub)
					% run the surrogate generation program in Matlab
					shufflesyncsurr(1,numwindows,Args.SurrData);
					computeSurrData(Args.SurrData);
					% create data strcture and instantiate object
					obj = makeObject(Args);
				else % if(Args.NoQsub)
					% save asxcdata binsc2 binidx2 winbinsize syncbins nsbins nsbins2 frames reps maxlag
					% call shell script to generate surrogates
					% create system command
					syscmd = [Args.QsubProgram ' ' num2str(numwindows) ' ' ...
						num2str(Args.JobWindows) ' ' Args.SurrDataFile];
					% check if we should wait for the surrogates to be generated
					% don't try to use nohup of to background the job since we want
					% to take a look at the output to make sure the jobs were 
					% submitted properly
					[s,w] = system(syscmd);
					% make sure shell script ran without problems
					if(s~=0)
						error([Args.classname ': Error creating surrogates!'])
					else
						% display output
						fprintf('%s\n',w);
					end
					if(Args.WaitForQsub)
						fprintf('Waiting for qsub ...');
						% check if result file has been generated
						while(isempty(nptDir(Args.ResultFile)))
							% wait for 10 seconds
							pause(Args.QsubWaitTime)
							fprintf('.');
						end
						fprintf('\n');
						% create data structure and instantiate object
						obj = makeObject(Args);
					else % if(Args.WaitForQsub)
						% call system command with & so it gets backgrounded and 
						% returns immediately
						% [s,w] = system([syscmd]);
						% create empty object
						obj = createEmptyObject(Args);
					end % if(Args.WaitForQsub)
				end % if(Args.NoQsub)
			end % if( isempty(adjs1) | isempty(adjs2) )
		else % if(~isempty(Args.ClusterDirs))
			% create empty object
			obj = createEmptyObject(Args);
		end % if(~isempty(Args.ClusterDirs))
	end % if(Args.SkipNoResult)
else % if(isempty(resfile))
	% create data structure and instantiate object
	obj = makeObject(Args);
end

function obj = createEmptyObject(Args)

data.numSets = 0;
data.Args = [];
data.binSize = [];
data.WinStartTime = [];
data.DataSyncSpikes = [];
data.ShiftPSyncSpikes = [];
data.SurrSyncMean = [];
data.SurrSyncStd = [];
data.SurrSyncHist = [];
data.DataZScore = [];
data.WindowIndex = [];
data.RasterMarkMatrices = [];
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);

function obj = makeObject(Args)

% load data and result file
load(Args.DataFile);
load(Args.ResultFile);
objdata.numSets = 1;
objdata.Args = Args;
objdata.binSize = binSize;
objdata.reps = reps;
objdata.WinStartTime = winstarttime;
objdata.DataSyncSpikes = data;
objdata.ShiftPSyncSpikes = shiftp;
objdata.SurrSyncMean = surrmean;
objdata.SurrSyncStd = surrstd;
objdata.SurrSyncHist = surrhist;
objdata.DataZScore = dzscore;
% save number of windows
objdata.WindowIndex = [0; numwindows];				
% save the synchronous spikes for the raster display
objdata.RasterMarkMatrices = { { {logical(markmat1) logical(spmarkmat1)} ...
			{logical(markmat2) logical(spmarkmat2)} } };
% create nptdata
n = nptdata(objdata.numSets,0,pwd);
d.data = objdata;
obj = class(d,Args.classname,n);
if(Args.SaveLevels)
	fprintf('Saving %s object...\n',Args.classname);
	eval([Args.matvarname ' = obj;']);
	% save object
	eval(['save ' Args.matname ' ' Args.matvarname]);
end
