function [obj, varargout] = surrpsth(varargin)
%@surrpsth Constructor function for surrpsth class
%   OBJ = surrpsth(varargin)
%
%   OBJ = surrpsth('auto') attempts to create a SURRPSTH object which
%   compares the surrogate PSTH to the PSTH from an adjspikes object.

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
	'SurrName','refsga');
Args.flags = {'Auto','ArgsOnly'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'surrpsth';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'df';

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
    return;
end

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
		else
			% no saved object so we will try to create one
			% pass varargin in case createObject needs to instantiate
			% other objects that take optional input arguments
			obj = createObject(Args,modvarargin{:});
		end
	end
end

function obj = createObject(Args,varargin)

% instantiate adjspikes object
adj = adjspikes('auto',varargin{:});
% get number of repetitions
ras1 = adj.data.raster;
reps = size(ras1,2);
% get frame points
framepts = adj.data.adjFramePoints;
% get number of frames
frames = (length(framepts) - 1) / reps;
fpts = framepts(1:(frames+1));
% subtract first frame point since it might not be zero
% while the rasters all have been adjusted to start at zero
fp = fpts-fpts(1);
% compute spike counts
sc = histcie(ras1,fp,'DropLast');
% compute mean and std of spike counts
scm = nanmean(sc,2);
scstd = nanstd(sc,0,2);
sem = scstd/sqrt(reps);

% find surrogate files
sflist = nptDir(Args.SurrName);
% get number of files
sfnum = size(sflist,2);
% initialize variables to compute running std
sumx = zeros(frames,reps);
sumx2 = sumx;
for sfi = 1:sfnum
    spt2 = readSurrogateBin('refsga1.bin');
    % get number of surrogates inside each file
    surrnum = size(spt2,2);
    for sni = 1:surrnum
        sc1 = histcie(cell2array(spt2{1}),fp,'DropLast');
        sumx = sumx + sum(sc1,2);
        sumx2 = sumx2 + sum(sc1.^2,2);
    end
end


% example object
dlist = nptDir;
% get entries in directory
dnum = size(dlist,1);

% check if the right conditions were met to create object
if(dnum>0)
	% this is a valid object
	% these are fields that are useful for most objects
	data.numSets = 1;
    data.Args = Args;
	
	% these are object specific fields
	data.dlist = dlist;
	% set index to keep track of which data goes with which directory
	data.setIndex = [0; dnum];
	
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
	% create empty object
	obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;
data.setNames = '';

% these are object specific fields
data.dlist = [];
data.setIndex = [];

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
