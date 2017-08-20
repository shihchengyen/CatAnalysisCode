function obj = trentropy(varargin)
%@trentropy/trentropy Constructor function for TRENTROPY object
%   OBJ = trentropy(VARARGIN)

Args = struct('RedoLevels',0,'SaveLevels',0, ...
			  'Auto',0,'CellName','','FrameBins',10);
			  
[Args,varargin] = getOptArgs(varargin,Args, ...
	'flags',{'Auto'}, ...
	'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
	'subtract',{'RedoLevels','SaveLevels'});

if nargin==0
	% create empty object
	obj = createEmptyObject;
elseif( (nargin==1) & (isa(varargin{1},'trentropy')) )
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
		% check for saved object
		if(ispresent('trentropy.mat','file','CaseInsensitive') ...
		  & (Args.RedoLevels==0))
			fprintf('Loading saved trentropy object...\n');
			% load saved object and exit
			l = load('trentropy.mat');
			obj = l.tre;
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
				obj = createObject(mat.data,Args.SaveLevels,fname, ...
						Args.FrameBins);
			end
		end
	else
		% hopefully we were given a data structure
		% create object using data structure
		obj = createObject(varargin{1},Args.SaveLevels,Args.CellName, ...
				Args.FrameBins);
	end
end

function obj = createEmptyObject

n = nptdata(0,0);
d.data.cellname = '';
d.data.framebins = [];
d.data.entropy = [];
d.data.surrPercent = [];
d.data.surrZScores = [];
obj = class(d,'trentropy',n);


function obj = createObject(data,savelevels,cname,framebins)

d.data.cellname = {cname};
d.data.framebins = framebins;
% compute number of frames presented
nframes = data.stimulus_info.end_frame-data.stimulus_info.start_frame+1;
% get bin limits for histcie
afv = vecc(data.stimulus_info.adjusted_frames_vector);
% get length of afv
afvl = length(afv);
dafv = diff(afv);
binsizes = dafv/framebins;
mat1 = tril(ones(framebins));
mat2 = [afv(1:(afvl-1))'; repmat(binsizes',(framebins-1),1)];
blimits = mat1 * mat2;
binlimits = [reshape(blimits,[],1); afv(afvl)];
sc = histcie(data.cell_info.adjusted_spiketrain, ...
	binlimits,'DropLast');
% reshape into repetitions and then transpose so we can take the mean 
% and std easily
sc1 = reshape(sc,framebins,[]);
dentropy = getTREntropy(sc1,nframes);
% keep track of cell id which will make it easier to identify cells
% when trentropy objects are added together
d.data.cellid = ones(nframes,1);
% load the data from the surrogates and compute percentage of surrogates
% that have higher entropy
load(['framesgTRE' num2str(framebins)])
% surrogate entropies in matrix with dimensions surrsets x nframes
surrSets = size(entropy,1);
% replicate data entropy
dataentropy = repmat(dentropy,surrSets,1);
d.data.entropy = dentropy';
% compute percent of surrogates with higher entropy than data for each 
% frame
d.data.surrPercent = (sum(dataentropy<entropy) / surrSets)';
% compute mean and std of surrogate data
surrMean = mean(entropy);
surrStd = std(entropy);
warning off MATLAB:divideByZero
d.data.surrZScores = ((dentropy - surrMean) ./ surrStd)';
warning on MATLAB:divideByZero
n = nptdata(1,0);
obj = class(d,'trentropy',n);

if(savelevels>0)
	fprintf('Saving trentropy object...\n');
	tre = obj;
	save trentropy tre
end
