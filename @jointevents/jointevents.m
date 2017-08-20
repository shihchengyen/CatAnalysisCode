function obj = jointevents(varargin)
% constructor for jointevents class
%   OBJ = jointevents('auto')
%
%   Must input the two or more ClusterDirs into a cell array 
%   to load the spiketrains and events objects.
%    
%dependencies:  adjspikes,events,stiminfo

Args = struct('RedoLevels',0,'SaveLevels',0, ...
    'Auto',0,'ClusterDirs',{''},'numRand',1000);

[Args,varargin] = getOptArgs(varargin,Args, ...
    'flags',{'Auto'}, ...
    'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
    'subtract',{'RedoLevels','SaveLevels'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'jointevents';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'je';

if(nargin==0)
    % create empty object
    obj = createEmptyObject(Args);
elseif( (nargin==1) & isa(varargin{1},Args.classname))
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
            obj = createObject(Args,varargin{:});
        end
    end
end

%%%%% Create an Empty Object %%%%%%
function obj = createEmptyObject(Args)
data.numSets = 0;
data.setNames = '';
data.jointevents = [];
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);

%%%%% Create the JointEvents Object %%%%
function obj = createObject(Args,varargin)
%%% jointevents object location %%%%%
cwd = pwd;
%%% Go into the directories and get the adjspikes and events objects %%%
for ii = 1:size(Args.ClusterDirs,2)
    cd(Args.ClusterDirs{ii})
    sp = adjspikes('auto');
    if isempty(sp)
        obj=createEmptyObject;
        return
    end
    JointEventsInfo(ii).spiketrain = sp.data.adjSpiketrain;
    JointEventsInfo(ii).framepoints = sp.data.adjFramePoints;
    ev = events('auto');
    if isempty(ev)
        obj=createEmptyObject;
        return
    end
    JointEventsInfo(ii).StartEvent = ev.data.events.StartEvent;
    JointEventsInfo(ii).EndEvent = ev.data.events.EndEvent;
    JointEventsInfo(ii).EventProbability = ev.data.events.EventProbability;
    JointEventsInfo(ii).Binsize = ev.data.events.Binsize;
    JointEventsInfo(ii).PSTH = ev.data.events.PSTH;
    JointEventsInfo(ii).Thresholds = ev.data.events.Thresholds;      
    % Calculate the PSTH at the frame resolution
    counts = mean((reshape(histcie(sp.data.adjSpiketrain,sp.data.adjFramePoints,'DropLast'),ev.data.events.NumFrames,ev.data.events.NumRepetitions))');
    JointEventsInfo(ii).FramePSTH = counts';
end

%%% Get the stiminfo %%
% get session directory
sessiondir = getDataDirs('session','relative');
% change to session directory
cd(sessiondir);
stimInfo = stiminfo('auto',1,'save');
if isempty(stimInfo)
    obj=createEmptyObject;
    return
end

%%% Change back to the jointevents object location %%
cd(cwd)

% this should be a valid object
data.numSets = 1;
data.setNames{1} = pwd;

%%% Private function that creates the jointevents data %%%%%
data.jointevents = JointEventsAnalysis(JointEventsInfo,stimInfo,Args);
d.data = data;

n = nptdata(data.numSets,0,pwd);
obj = class(d,Args.classname,n);
if(Args.SaveLevels)
    fprintf('Saving %s object...\n',Args.classname);
    eval([Args.matvarname ' = obj;']);
    % save object
    eval(['save ' Args.matname ' ' Args.matvarname]);
end