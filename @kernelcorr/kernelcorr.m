function obj = kernelcorr(varargin)
% constructor for kernelcorr class
%   OBJ = kernelcorr('auto')
%
%   Must input the two or more ClusterDirs into a cell array 
%   to load the spiketrains.
%    
% dependencies:  revcorr

Args = struct('RedoLevels',0,'SaveLevels',0, ...
    'Auto',0,'ClusterDirs',{''},'DegreeDiameter',6,...
    'NumFrames',6,'FrameNumber',0,'Plot',1,'Threshold',0);

[Args,varargin] = getOptArgs(varargin,Args, ...
    'flags',{'Auto'}, ...
    'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
    'subtract',{'RedoLevels','SaveLevels'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'kernelcorr';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'kc';

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
data.KernelCorr = [];
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);

%%%%% Create the KernelCorr Object %%%%
function obj = createObject(Args,varargin)
%%% jointevents object location %%%%%
cwd = pwd;
%%% Go into the directories and get the adjspikes
cd(Args.ClusterDirs{1})
rf1 = revcorr('auto');
if isempty(rf1)
    obj=createEmptyObject;
    return
end
rf1 = rf1.data.R;
if size(rf1,3) ~= Args.NumFrames
    rf1 = rf1(:,:,1:Args.NumFrames);
end
if Args.FrameNumber
    rf1 = rf1(:,:,Args.FrameNumber);
end
cd(Args.ClusterDirs{2})
rf2 = revcorr('auto');
if isempty(rf2)
    obj=createEmptyObject;
    return
end
rf2 = rf2.data.R;
if size(rf2,3) ~= Args.NumFrames
    rf2 = rf2(:,:,1:Args.NumFrames);
end
if Args.FrameNumber
    rf2 = rf2(:,:,Args.FrameNumber);
end

%%% Get the stiminfo %%
sessiondir = dirLevel('session','relative');
% change to session directory
cd(sessiondir);
ses_dir = pwd;
session_name = [ses_dir(end-18:end-17) ses_dir(end-1:end)]; % Get the session name
stimInfo = stiminfo('auto');
if isempty(stimInfo)
    obj=createEmptyObject;
    return
end

%%% Change back to the combo directory object location %%
cd(cwd)
data.numSets = 1;
data.setNames{1} = pwd;

%%% Private function that creates the kernelcorr data %%%%%
[data.KernelCorr,data.KernelPVal] = brunokernelcorr(rf1,rf2,session_name,Args);
%[data.KernelCorr,data.KernelPVal] = calckernelcorr(rf1,rf2,stimInfo,Args);
d.data = data;

n = nptdata(data.numSets,0,pwd);
obj = class(d,Args.classname,n);
if(Args.SaveLevels)
    fprintf('Saving %s object...\n',Args.classname);
    eval([Args.matvarname ' = obj;']);
    % save object
    eval(['save ' Args.matname ' ' Args.matvarname]);
end