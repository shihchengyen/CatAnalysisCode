function obj = psparse(varargin)
%@psparse Constructor function for psparse class
%   OBJ = psparse('auto') attempts to create a psparse object by ...

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0, ...
                'ClusterDirs',{''});
Args.flags = {'Auto'};
[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'psparse';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'ps';

numArgin = nargin;
if(numArgin==0)
	% create empty object
	obj = createEmptyObject(Args);
elseif( (numArgin==1) & isa(varargin{1},Args.classname))
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
        % change to the proper directory
        % [pdir,cdir] = getDataDirs('session','relative','CDNow');
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
			obj = createObject(Args,modvarargin{:});
		end
        % change back to previous directory if necessary
        % if(isempty(cdir))
        %     cd(cdir)
        % end
	end
end

function obj = createObject(Args,varargin)

% if ClusterDirs is empty, try to generate it from the current directory
if(isempty(Args.ClusterDirs))
    Args.ClusterDirs = getDataDirs('GetClusterDirs');
end
% get number of directories
dsize = size(Args.ClusterDirs,2);
if(dsize>1)
	% get current directory
	cwd = pwd;
	for idx = 1:dsize
		% change to cluster directory
		cd(Args.ClusterDirs{idx})
		frobj = firingrate('auto',varargin{:});
		if(idx==1)
			% get size and instantiate appropriate matrix
			fr = frobj.data.firingRate;
            % noisecounts = frobj.data.meancounts;
            sm = frobj.data.spike_matrix;
			[reps,nbins] = size(sm);
			% length of noise counts not the same as nbins since it is 
			% number of repetitions * nbins
			% lnc = length(noisecounts);
			frmat = zeros(dsize,nbins);
			% smat = cell(1,dsize);
            % noisemat = zeros(dsize,lnc);
			% transpose firing rate so we can take difference in firing
			% rates for each bin easily
			frmat(1,:) = fr';
            % noisemat(1,:) = noisecounts';
            % smat{1} = sm;
            timebins = frobj.data.timebins;
		else
			frmat(idx,:) = frobj.data.firingRate';
            % noisemat(idx,:) = frobj.data.meancounts';
            % smat{idx} = frobj.data.spike_matrix;
		end
		cd(cwd);
	end
	% take the sum for each column
	csum = sum(frmat);
	% the sum of the square of each column
	c2sum = sum(frmat.*frmat);
    warning off MATLAB:divideByZero
    ds1 = 1/dsize;
	dfr = (1 - ds1 * ((csum .* csum)./ c2sum))/(1-ds1);
    warning on MATLAB:divideByZero

	% create object data
	% this is a valid object
	% these are fields that are useful for most objects
	data.numSets = 1;
	data.setNames{1} = cwd;
	
	data.psparse = vecc(dfr);
    data.timebins = vecc(timebins);
    data.ncells = dsize;
	
	% create nptdata so we can inherit from it
	n = nptdata(data.numSets,0,cwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	if(Args.SaveLevels)
		fprintf('Saving %s object...\n',Args.classname);
		eval([Args.matvarname ' = obj;']);
		% save object
		eval(['save ' Args.matname ' ' Args.matvarname]);
	end	
else
	fprintf('Warning: %s requires at least 2 cells!\n', ...
		Args.classname);
	obj = createEmptyObject(Args);
end
	
function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;
data.setNames = '';

% these are object specific fields
data.psparse = [];
data.timebins = [];
data.ncells = 0;

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
