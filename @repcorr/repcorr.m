function obj = repcorr(varargin)
%@dirfiles Constructor function for DIRFILES class
%   OBJ = dirfiles('auto') attempts to create a DIRFILES object by ...

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0,'Grating',0, ...
	'ClusterDirs',{''},'HistBinSize',40,'IntraInter',0);
Args.flags = {'Auto','Grating','IntraInter'};
[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'repcorr';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'rc';

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
		if(~isempty(nptDir(Args.matname,'CaseInsensitive')) ...
			&& (Args.RedoLevels==0))
			fprintf('Loading saved %s object...\n',Args.classname);
			l = load(Args.matname);
			obj = eval(['l.' Args.matvarname]);
            % check if the saved object's Args structure matches the
            % current Args structure
            
		else
			% no saved object so we will try to create one
			% pass varargin in case createObject needs to instantiate
			% other objects that take optional input arguments
			obj = createObject(Args,modvarargin{:});
		end
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
	if(Args.Grating)
		adjs = cell(dsize,1);
		for idx = 1:dsize
			% change to cluster directory
			cd(Args.ClusterDirs{idx})
			adjs{idx} = adjspikes('auto',varargin{:});
		end
		reps = size(adjs{1}.data.raster,2);
		repindices = (1:reps)';
		rasters = concat(adjs{1}.data.raster,adjs{2}.data.raster,'Columnwise');
		maxsptime = max(rasters(:));
		histbinsize = Args.HistBinSize;
		edges = 0:histbinsize:(ceil(maxsptime/histbinsize)*histbinsize);
		smat = histcie(rasters,edges);
	else % if(Args.Grating)
		for idx = 1:dsize
			% change to cluster directory
			cd(Args.ClusterDirs{idx})
			frobj = firingrate('auto',varargin{:});
			sm = frobj.data.spike_matrix;
			if(idx==1)
				[reps,nbins] = size(sm);
				% set up indices so we don't have to do this inside the loop
				repindices = (1:reps)';
				% set up array to store reps of all cells
				smat = zeros(reps*dsize,nbins);
				smat(repindices,:) = sm;
			else
				smat((idx-1)*reps+repindices,:) = sm;
			end
		end % for idx = 1:dsize
        % transpose smat so spike counts are in columns
        % this will allow corrcoef to be computed correctly
        smat = smat';
	end % if(Args.Grating)
	r = corrcoef(smat);
    if(Args.IntraInter)
        reps2 = repindices + reps;
        c1 = r(repindices,repindices);
        c2 = r(reps2,reps2);
        c3 = r(reps2,repindices);
        ctrilind = logical(tril(ones(reps,reps),-1));
        c1v = c1(ctrilind);
        c2v = c2(ctrilind);
        c = concat(c1v,c2v,c3(:),'Columnwise');
    else
        % find correlations between corresponding repetitions
        reps2 = 2 * reps;
        c = r(sub2ind([reps2 reps2],repindices,repindices+reps));
    end
    data.numSets = 1;
    data.Args = Args;
    data.repcorr = c;
		
	% create nptdata so we can inherit from it
	n = nptdata(data.numSets,0,pwd);
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
	% create empty object
	obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;
data.Args = [];
data.repcorr = [];

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
