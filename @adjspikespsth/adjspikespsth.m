function obj = adjspikespsth(varargin)
% constructor for adjspikespsth class
%   OBJ = adjspikespsth(varargin) 
%   
%   Calculates the PSTH from the adjusted
%   spike trians.  Only works 
%   on the cell level.
%   
%   example as = adjspikespsth('save','redo')
% 
%dependencies:  adjspikes, stimInfo

Args = struct('RedoLevels',0,'SaveLevels',0, ...
    'Auto',0,'stimInfoDir',['..' filesep '..' ],'binsize','frame','overlap',1);

[Args,varargin] = getOptArgs(varargin,Args, ...
    'flags',{'Auto'}, ...
    'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
    'subtract',{'RedoLevels','SaveLevels'});

Args.classname = 'adjspikespsth';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'ap';

if nargin==0
    % create empty object
    obj = createEmptyObject(Args);
elseif( (nargin==1) & (isa(varargin{1},'adjspikespsth'))  )
    obj = varargin{1};
elseif Args.RedoLevels==0
    % check for saved object
    dirlist = nptDir('adjspikespsth.mat','CaseInsensitive');
    if size(dirlist,1)> 1
        warning('More than one adjspikespsth file found.  Loading %s\n',dirlist(1).name)
    end
    if ~isempty(dirlist)
        fprintf('Loading saved adjspikespsth object...\n');
        % load saved object and exit
        try
            l = load(dirlist(1).name);
            obj = l.ap;
        catch
            fprintf('older format... Recalculating...\n')
            obj = createObject(group,Args);
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

[data.psth,data.plotVector] = adjPSTH(sp,stimInfo,Args.binsize,Args.overlap);
d.data = data;

number = get(sp,'Number');
holdaxis = get(sp,'HoldAxis');
nd = nptdata(number,holdaxis);
obj = class(d,'adjspikespsth',nd);
if(Args.SaveLevels>0)
    fprintf('Saving adjspikespsth object...\n');
    filename = 'adjspikespsth.mat';
    ap = obj;
    save(filename,'ap')
end
