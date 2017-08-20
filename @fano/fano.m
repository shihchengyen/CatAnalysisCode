function obj = fano(varargin)
%@fano/fano Constructor function for FANO object
%   OBJ = fano(VARARGIN)

Args = struct('RedoLevels',0,'SaveLevels',0, ...
			  'Auto',0,'CellName','');
			  
[Args,varargin] = getOptArgs(varargin,Args, ...
	'flags',{'Auto'}, ...
	'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
	'subtract',{'RedoLevels','SaveLevels'});

if nargin==0
	% create empty object
	obj = createEmptyObject;
elseif( (nargin==1) & (isa(varargin{1},'fano')) )
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
		% check for saved object
		if(ispresent('fano.mat','file','CaseInsensitive') & (Args.RedoLevels==0))
			fprintf('Loading saved fano object...\n');
			% load saved object and exit
			l = load('fano.mat');
			obj = l.ff;
		else
			% no saved object so try to create one
			% try to find name of .mat file
			a = nptDir('*g*c*.mat');
			if(isempty(a))
				% create empty object
				obj = createEmptyObject;
			else
				mat = load(a(1).name);
				[fpath,fname] = nptFileParts(a(1).name);
				obj = createObject(mat.data,Args.SaveLevels,fname);
			end
		end
	else
		% hopefully we were given a data structure
		% create object using data structure
		obj = createObject(varargin{1},Args.SaveLevels,Args.CellName);
	end
end

function obj = createEmptyObject

n = nptdata(0,0);
d.data.cellname = '';
d.data.scmean = [];
d.data.scstd = [];
d.data.cellid = [];
d.data.fano = [];
d.data.surrPercent = [];
d.data.surrZScores = [];
obj = class(d,'fano',n);

function obj = createObject(data,savelevels,cname)

d.data.cellname = {cname};
sc = histcie(data.cell_info.adjusted_spiketrain, ...
	data.stimulus_info.adjusted_frames_vector,'DropLast');
% compute number of frames presented
nframes = data.stimulus_info.end_frame - data.stimulus_info.start_frame + 1;
% reshape into repetitions and then transpose so we can take the mean and std easily
sc1 = reshape(sc,nframes,[]);
d.data.scmean = mean(sc1,2);
d.data.scstd = std(sc1,0,2);
warning off MATLAB:divideByZero
d.data.fano = (d.data.scstd.^2) ./ d.data.scmean;
% keep track of cell id which will make it easier to identify cells
% when Fano objects are added together
d.data.cellid = ones(nframes,1);
% load the data from the surrogates and compute percentage of surrogates
% that have higher Fano
load framesgFF
% compute Fano for surrogates, dimensions are (numframes x surrsets)
surrFano = (sgSC.scstd.^2) ./ sgSC.scmean;
% get number of sets of surrogates
surrSets = size(surrFano,2);
% replicate data Fano
dataFano = repmat(d.data.fano,1,surrSets);
% compute percent of surrogates with higher Fano than data for each frame
d.data.surrPercent = sum(dataFano<surrFano,2) / surrSets;
% get mean and std of surrogate data
surrMean = mean(surrFano,2);
surrStd = std(surrFano,0,2);
d.data.surrZScores = (d.data.fano - surrMean) ./ surrStd;
n = nptdata(1,0);
obj = class(d,'fano',n);
if(savelevels>0)
	fprintf('Saving fano object...\n');
	ff = obj;
	save fano ff
end
warning on MATLAB:divideByZero