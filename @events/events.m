function obj = events(varargin)
% constructor for events class
%   OBJ = events(varargin) 
%   %   
%example es = events('save','redo')
% 
%dependencies:  adjspikes,stiminfo

Args = struct('RedoLevels',0,'SaveLevels',0, ...
    'Auto',0,'stimInfoDir',['..' filesep '..' ],....
    'Binsize',10,'ThresholdType','percentile', ...
    'SlideDuration',1,'PercentageReps',33, ...
    'Upper',95,'Lower','mean');

[Args,varargin] = getOptArgs(varargin,Args, ...
    'flags',{'Auto'}, ...
    'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
    'subtract',{'RedoLevels','SaveLevels'});

Args.classname = 'events';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'ev';

if nargin==0
    % create empty object
    obj = createEmptyObject(Args);
elseif( (nargin==1) & (isa(varargin{1},'adjspikes'))  )
    obj = varargin{1};
elseif Args.RedoLevels==0
    % check for saved object
    dirlist = nptDir('events.mat');
    if ~isempty(dirlist)
        fprintf('Loading saved events object...\n');
        % load saved object and exit
        try
            l = load(dirlist(1).name);
            obj = l.ev;
        catch
            fprintf('older format... Recalculating...')
            obj = createObject(Args);
        end
    else
        % no saved object so try to create one
        obj = createObject(Args);
    end
elseif Args.RedoLevels
    obj = createObject(Args);
end

function obj = createEmptyObject(Args)
data.numSets = 0;
data.setNames = '';
data.dlist = [];
data.setIndex = [];
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);

function obj = createObject(Args)
%get ispikes
sp = adjspikes('auto');
if isempty(sp)
    obj=createEmptyObject;
    return
end
%get stiminfo
%need to go up to session directory so...
pdir = pwd;
cd(Args.stimInfoDir)
stimInfo = stiminfo('auto');
cd(pdir)
data.numSets = 1;
data.setNames{1} = pwd;
data.dlist = nptDir;
%%%% Calculate the events analysis %%%%
data.events = EventsAnalysis(sp,stimInfo,Args);
d.data = data;
number = 1;
holdaxis = get(sp,'HoldAxis');
nd = nptdata(number,holdaxis);
obj = class(d,'events',nd);
if(Args.SaveLevels>0)
    fprintf('Saving events object...\n');
    filename = 'events.mat';
    ev = obj;
    save(filename,'ev')
end