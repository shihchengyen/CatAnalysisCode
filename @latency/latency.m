function obj = dirfiles(varargin)
%@latency Constructor function for LATENCY class
%   OBJ = latency('auto') attempts to create a LATENCY object by ...

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0,'SubBinSize',1, ...
    'NumFrames',1);
Args.flags = {'Auto'};
[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'latency';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'la';

numArgin = nargin;
if(numArgin==0)
	% create empty object
	obj = createEmptyObject(Args);
elseif( (numArgin==1) & isa(varargin{1},Args.classname))
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
		% check for saved object
		if(ispresent(Args.matname,'file','CaseInsensitive') ...
			& (Args.RedoLevels==0))
			fprintf('Loading saved %s object...\n',Args.classname);
			l = load(Args.matname);
			obj = eval(['l.' Args.matvarname]);
		else
			% no saved object so we will try to create one
			% pass varargin in case createObject needs to instantiate
			% other objects that take optional input arguments
			obj = createObject(Args,modvarargin{:});
		end
	end
end

function obj = createObject(Args,varargin)

% instantiate adjspikes object since we want to make sure all the frame 
% windows are the same length
asp = adjspikes('auto',varargin{:});
% get path to session dir
sdir = getDataDirs('session','Relative');
% save current directory
cwd = pwd;
% change to session dir
cd(sdir);
% instantiate stiminfo object
stm = stiminfo('auto',varargin{:});
% move back to original directory
cd(cwd);
% check if ispikes object is empy
if( ~isempty(asp) || ~isempty(stm) )
	data.SubBinSize = Args.SubBinSize;
	% divide frame vector into requested size
	[subbins,data.ActualSubBinSize,nsbins] = divideBins( ...
		asp.data.adjFramePoints,'SubBinSize',Args.SubBinSize);
	% take histogram using requested bin size
	binsc = histcie(asp.data.adjSpiketrain,subbins,'DropLast');
	% reshape into nframes*nsbins x floor(frames*reps/nframes)
    % figure out indices in binsc to use
    % nframes = Args.NumFrames;
    % binscend = nframes*nsbins*(floor(nframes*reps/
	binsc1 = reshape(binsc,Args.NumFrames*nsbins,[]);
	% take the sum of each bin and append 0 at the end so we can use the
	% 'histc' option in bar
	data.BinSums{1} = [sum(binsc1,2); 0];
	% store refresh rate
	data.refreshRate = stm.data.catInfo.video_refresh;
	% create nptdata
	n = nptdata(1,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	if(Args.SaveLevels)
		fprintf('Saving %s object...\n',Args.classname);
		eval([Args.matvarname ' = obj;']);
		% save object
		eval(['save ' Args.matname ' ' Args.matvarname]);
	end
else
	% one of the required objects is empty so create empty object
	obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)

data.SubBinSize = [];
data.ActualSubBinSize = [];
data.BinSums = {};
data.refreshRate = [];
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
