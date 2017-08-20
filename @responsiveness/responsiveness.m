function obj = responsiveness(varargin)
% constructor for responsiveness class
%   OBJ = adjISpikes(varargin) 
%   
%   
%example as = responsiveness('save','redo')
% 
%dependencies:  adjspikes,stiminfo

Args = struct('RedoLevels',0,'SaveLevels',0, ...
    'Auto',0,'Significance',0.05,'PlotDist',0,'RandTimeTrial',0, ...
    'RandSets',1000);

[Args,varargin] = getOptArgs(varargin,Args, ...
    'flags',{'Auto','PlotDist','RandTimeTrial'}, ...
    'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
    'subtract',{'RedoLevels','SaveLevels'});

Args.classname = 'responsiveness';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'r';

if nargin==0
    % create empty object
    obj = createEmptyObject(Args);
elseif( (nargin==1) & (isa(varargin{1},Args.classname))  )
    obj = varargin{1};
elseif Args.RedoLevels==0
    % check for saved object
    dirlist = nptDir(Args.matname);
    if ~isempty(dirlist)
        fprintf('Loading saved responsiveness object...\n');
        % load saved object and exit
        try
            l = load(dirlist(1).name);
            obj = l.r;
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
data.responsiveness = nan;
data.corrcoef = nan;
data.corrp = nan;
data.ctrlcoeff = nan;
data.ctrlp = nan;
data.distdiffr = nan;
data.distdiffp = nan;
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
data.numSets = 1;
data.setNames{1} = pdir;

[data.responsiveness,data.corrcoef,data.corrp,data.ctrlcoeff, ...
    data.ctrlp,data.distdiffr,data.distdiffp] = CalcResponsiveness(sp,Args);
d.data = data;

number = 1;
holdaxis = get(sp,'HoldAxis');
nd = nptdata(number,holdaxis,pdir);

obj = class(d,Args.classname,nd);
if(Args.SaveLevels>0)
    fprintf('Saving responsiveness object...\n');
    filename = 'responsiveness.mat';
    r = obj;
    save(filename,'r')
end
