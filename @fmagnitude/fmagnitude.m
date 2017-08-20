function obj = fmagnitude(varargin)
%fmagnitude Constructor for fmagnitude objects
%   OBJ = fmagnitude('auto',VARARGIN) computes the amplitude in the 
%   frequency spectrum of the spike trains contained in the MAT file 
%   stored in the local directory. 
%
%   OBJ = fmagnitude(MAT,'CellName',NAME,VARARGIN) computes the FFT 
%   amplitude of the data structure MAT. It uses the NAME argument as 
%   the name of the cell.
%
%   The following are the optional input arguments:
%      RepsPerStep - number of repetitions that are passed to the FFT
%                    function at one time (default: 6).

Args = struct('RepsPerStep',6,'RedoLevels',0,'SaveLevels',0, ...
			  'Auto',0,'CellName','');
[Args,modvarargin] = getOptArgs(varargin,Args,'flags',{'Auto'}, ...
	'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'remove',{'Auto','RepsPerStep','CellName'});

data.nCells = 0;
data.cellname = '';
data.mean = [];
data.stdev = [];
data.f = [];
data.fftMag = [];
data.colIndex = [];
data.timeVec = [];

numArgin = nargin;
if(numArgin==0)
	n = nptdata(0,0);
	d.data = data;
	obj = class(d,'fmagnitude',n);
elseif( (numArgin==1) & (isa(varargin{1},'fmagnitude')) )
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
		% check for saved object
		if(ispresent('fmagnitude.mat','file','CaseInsensitive') ...
			& (Args.RedoLevels==0))
			fprintf('Loading saved fmagnitude object...\n');
			% load saved object and exit
			l = load('fmagnitude.mat');
			obj = l.fm;
		else
			% no saved object so try to create one
			% first try to instantiate ispikes object
            isp = adjspikes('auto',modvarargin{:});
            if(~isempty(isp))
                obj = createObject(Args,isp,modvarargin{:});
            else
				% try to find name of .mat file
				a = nptDir('*g*c*.mat');
                if(isempty(a))
                    % try data.mat instead
                    a = nptDir('data.mat');
                end
				if(isempty(a))
					% create empty object
					n = nptdata(0,0);
					d.data = data;
					obj = class(d,'fmagnitude',n);
				else
					mat = load(a(1).name);
                    if(strcmp(a(1).name,'data.mat'))
                        fname = strrep(pwd,'/Users/syen/Documents/ShihCheng/Data/Neural/Cat/catdata/','');
                    else
        				[fpath,fname] = nptFileParts(a(1).name);
                    end
					Args.CellName = fname;
					obj = createObject(Args,mat.data,modvarargin{:});
				end
            end
		end
	else
		% create using input data structure
		obj = createObject(Args.NumericArguments{1},Args,modvarargin);
	end
end

function obj = createObject(Args,data,varargin)

d.data.nCells = 1;

if(isa(data,'adjspikes'))
    d.data.cellname = {pwd};
    [mx,f,t] = nptSpikeTimesFFT(data.data.raster,varargin{:});
    d.data.mean = mean(mx,2);
    d.data.stdev = std(mx,0,2);
    d.data.f = f;
    % get length of t vector
    tlength = length(t);
    d.data.fftMag = mean(reshape(mx,length(f),tlength,[]),3);
    d.data.colIndex = [0; tlength];
    d.data.timeVec = t;
else        
	d.data.cellname = Args.CellName;
	% get number of repetitions
	reps = length(data.cell_info.raster);
	% get rep numbers for each step
	stepreps = 0:Args.RepsPerStep:reps;
	% if last number is not reps, add it
	if(stepreps(end)~=reps)
		stepreps = [stepreps reps];
	end
	mx = [];
	% concatenate reps together before calling nptSpikeTimes
	for j = 1:(length(stepreps)-1)
		b = [];
		for k = (stepreps(j)+1):stepreps(j+1)
            if(~isempty(data.cell_info.raster{k}))
        		b = concatenate(b,data.cell_info.raster{k});
            else
                b = concatenate(b,NaN);
            end
		end
		% make sure b is not empty, if so skip these reps
		if(~isempty(b))
			% compute FFT for each column
			[mx2,f,tvec] = nptSpikeTimesFFT(b',varargin{:});
			mx = [mx mx2];
		end
	end
	d.data.mean = mean(mx,2);
	d.data.stdev = std(mx,2);
	d.data.f = f;
    d.data.fftMag = mx;
    d.data.colIndex = [0; size(mx,2)];
    d.data.timeVec = tvec;
end
n = nptdata(1,0,pwd);
obj = class(d,'fmagnitude',n);

if(Args.SaveLevels>0)
	fprintf('Saving fmagnitude object...\n');
	fm = obj;
	save fmagnitude fm
end
