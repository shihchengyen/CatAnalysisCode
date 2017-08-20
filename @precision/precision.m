function [obj, varargout] = precision(varargin)

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0,...
    'plottingFano', 0, 'plottingEntropy', 0);
Args.flags = {'Auto','ArgsOnly'};
% The arguments which can be neglected during arguments checking
Args.UnimportantArgs = {'RedoLevels','SaveLevels'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
    'subtract',{'RedoLevels','SaveLevels'}, ...
    'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
    'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'precision';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'a';

% cd ~/data;
if(nargin==0)
    % create empty object
    obj = createEmptyObject(Args);
elseif( (nargin==1) & isa(varargin{1},Args.classname))
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
            obj = createObject(Args,varargin{:});
        end
    end
end

function obj = createObject(Args,varargin)
cd /opt/data/cat;
load ndobjs;
index = 1;
for(i = 1: length(ndcc.SessionDirs))
    cd (ndcc.SessionDirs{i});
    load('refsga1FF.mat');
    index = [index index(end)+size(sg.scmean, 1)];
end
cd ~/data;

%MODIFY here to add additional variability objects
if(ispresent('variabilitya.mat','file')&&ispresent('variabilityb.mat','file')...
        &&ispresent('variabilityc.mat','file')&&ispresent('variabilityd.mat','file')...
        &&ispresent('variabilitye.mat','file')&&ispresent('variabilityf.mat','file')...
        &&ispresent('variabilityg.mat','file')&&ispresent('variabilityh.mat','file')...
        &&ispresent('variabilityi.mat','file')&&ispresent('variabilityj.mat','file'))
    load ('variabilitya.mat');
    load ('variabilityb.mat');
    load ('variabilityc.mat');
    load ('variabilityd.mat');
    load ('variabilitye.mat');
    load ('variabilityf.mat');
    load ('variabilityg.mat');
    load ('variabilityh.mat');
    load ('variabilityi.mat');
    load ('variabilityj.mat');
    objstr = {va vb vc vd ve vf vg vh vi vj};
    %#################################################

    fanomat = zeros(length(vh.data.scmean), length(objstr));
    tremat  = fanomat;
    for i = 1:length(objstr)

        msc1idx = objstr{i}.data.scmean>=1;
        fspidx = objstr{i}.data.fanoSurrPercent>=0.95;
        fanomat(:,i) = msc1idx & fspidx;

        nspidx = objstr{i}.data.entropySurrPercent>=0.95;
        tremat(:,i) = msc1idx & nspidx;
    end
    fanopre = sum(fanomat');
    trepre = sum(tremat');
    data.fanopre = fanopre;
    data.trepre = trepre;
    data.numSets = 51;
    data.index = index;

    % create nptdata so we can inherit from it

    data.Args = Args;
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
data.fanopre = [];
data.trepre = [];

data.numSets = 0;

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
