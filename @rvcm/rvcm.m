function obj = rvcm(varargin)
% constructor for rapid visual cortical modification (rvcm) class
%   OBJ = rvcm(varargin) 
%      
%example rvcm = rvcm('auto','redo','save')
% 
%dependencies:  adjspikes

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0,'Significance',0.05,...
    'PlotDist',0,'RandTimeTrial',0,'RandSets',10000,'Original',0);
Args.flags = {'Auto','PlotDist','RandTimeTrial','Original'};

[Args,varargin] = getOptArgs(varargin,Args, ...
    'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
    'subtract',{'RedoLevels','SaveLevels'});

Args.classname = 'rvcm';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'rv';

if nargin==0
    % create empty object
    obj = createEmptyObject(Args);
elseif( (nargin==1) && (isa(varargin{1},Args.classname))  )
    obj = varargin{1};
elseif Args.RedoLevels==0
    % check for saved object
    dirlist = nptDir(Args.matname);
    if ~isempty(dirlist)
        fprintf('Loading saved rvcm object...\n');
        % load saved object and exit
        try
            l = load(dirlist(1).name);
            obj = l.rv;
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
data.rv = nan;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);

function obj = createObject(Args)
sp = adjspikes('auto');
if isempty(sp)
    obj=createEmptyObject;
    return
end

%need to go up to session directory so...

d.data = CalcRVCM(sp,Args);

number = 1;
holdaxis = get(sp,'HoldAxis');
nd = nptdata(number,holdaxis,d.data.setNames{1});

obj = class(d,Args.classname,nd);
if(Args.SaveLevels>0)
    fprintf('Saving rvcm object...\n');
    filename = 'rvcm.mat';
    rv = obj;
    save(filename,'rv')
end=