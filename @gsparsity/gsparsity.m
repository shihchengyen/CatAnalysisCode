function obj = gsparsity(varargin)
%@gsparsity Constructor function for gsparsity class
%   OBJ = gsparsity('auto') attempts to create a gsparsity object by ...

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0,...
             'stimInfoDir',['..' filesep '..' ],'ClusterDirs',{''},...
             'Normalize',0,'Sum',0,'Repetitions',0,'Synchrony',0,...
             'Binsize',5,'SpikeCounts',0);
         
Args.flags = {'Auto','Repetitions','Normalize','Sum','Synchrony','SpikeCounts'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'gsparsity';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'gs';

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
        % change back to previous directory if necessary
        % if(isempty(cdir))
        %     cd(cdir)
        % end
	end
end

function obj = createObject(Args,varargin)
% get number of directories
dsize = size(Args.ClusterDirs,2);
if(dsize>1)
	% get current directory
	cwd = pwd;
    cd(Args.stimInfoDir) %% Change to get the stiminfo
    st = stiminfo('auto');
	for idx = 1:dsize
		% change to cluster directory
		cd(Args.ClusterDirs{idx})
        sp = adjspikes('auto');
		if(idx==1)
			% get size and instantiate appropriate matrix
			Spikes = sp.data.adjSpiketrain';
            FramePoints = sp.data.adjFramePoints';			
		else
			Spikes = concatenate(Spikes,sp.data.adjSpiketrain','ColumnWise');
        end
		cd(cwd);
	end	    
    [S,Values] = GroupSparseness(Args,Spikes,FramePoints,st);
    data.groupsparsity = S;
    data.groupValues = Values;
	% create object data
	% this is a valid object
	% these are fields that are useful for most objects
	data.numSets = 1;
	data.setNames{1} = cwd;
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
data.gsparsity = [];

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
