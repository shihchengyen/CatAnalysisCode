function obj = sparsity(varargin)
% constructor for sparsity class
%   OBJ = adjISpikes(varargin) 
%      
%example as = sparsity('save','redo')
% 
%dependencies:  adjspikes,stiminfo

Args = struct('RedoLevels',0,'SaveLevels',0, ...
    'Auto',0,'stimInfoDir',['..' filesep '..' ],'Repetitions',0, ...
    'WindowSize','frame','IgnoreFrames',0,'Rate',0,'ISpikes',0);

Args.flags = {'Auto','Repetitions','IgnoreFrames','Rate','ISpikes'};

[Args,varargin] = getOptArgs(varargin,Args,...
    'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
    'subtract',{'RedoLevels','SaveLevels'});

Args.classname = 'sparsity';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 's';

if nargin==0
    % create empty object
    obj = createEmptyObject(Args);
elseif( (nargin==1) & (isa(varargin{1},'adjspikes'))  )
    obj = varargin{1};
elseif Args.RedoLevels==0
    % check for saved object
    dirlist = nptDir('sparsity.mat');
    if ~isempty(dirlist)
        fprintf('Loading saved sparsity object...\n');
        % load saved object and exit
        try
            l = load(dirlist(1).name);
            obj = l.s;
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
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);

function obj = createObject(Args)
%get spikes object
if Args.ISpikes
    sp = ispikes('auto');
else
    sp = adjspikes('auto');
end
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

[S,Values] = Sparseness(sp,stimInfo,Args);
data.sparsity = S;
data.Values = Values;
data.Args = Args;
d.data = data;

number = 1;
holdaxis = get(sp,'HoldAxis');
nd = nptdata(number,holdaxis,pdir);

obj = class(d,'sparsity',nd);
if(Args.SaveLevels>0)
    fprintf('Saving sparsity object...\n');
    filename = 'sparsity.mat';
    s = obj;
    save(filename,'s')
end
