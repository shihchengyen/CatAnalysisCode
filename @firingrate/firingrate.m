function obj = firingrate(varargin)
% constructor for firingrate class
%   OBJ = firingrate(varargin)
%
% example as = firingrate('save','redo')
%
% dependencies:  adjspikes, stiminfo

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0,...
    'Binsize','frame','Repetitions',0,'Fit',0,'BootStrap',0,...
    'NumBoots',1000,'Surrogates',0,'Probability',0,'TossZeros',0,...
    'numBins',100,'HistRate',5,'Rate',0,'Overlap',0,...
    'Count',0,'AdjR',0,'NoiseCorr',0);

Args.flags = {'Auto','Rate','Repetitions','BootStrap','Fit',...
    'Surrogates','Probability','TossZeros','AdjR','NoiseCorr'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
    'subtract',{'RedoLevels','SaveLevels'}, ...
    'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
    'remove',{'Auto'});

Args.classname = 'firingrate';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'fr';

if nargin==0
    % create empty object
    obj = createEmptyObject(Args);
elseif((nargin==1) & (isa(varargin{1},'adjspikes')))
    obj = varargin{1};
elseif Args.RedoLevels==0
    % check for saved object
    dirlist = nptDir('firingrate.mat','CaseInsensitive');
    if ~isempty(dirlist)
        fprintf('Loading saved firingrate object...\n');
        % load saved object and exit
        try
            l = load(dirlist(1).name);
            obj = l.fr;
        catch
            fprintf('older format... Recalculating...')
            obj = createObject(Args,modvarargin{:});
        end
    else
        % no saved object so try to create one
        obj = createObject(Args,modvarargin{:});
    end
elseif Args.RedoLevels
    obj = createObject(Args,modvarargin{:});
end

function obj = createEmptyObject(Args)
data.numSets = 0;
data.setNames = '';
data.firingRate=[];
data.meancounts = [];
data.spike_matrix = [];
data.timebins = [];
data.binsize = [];
data.Fit=[];
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);

function obj = createObject(Args,varargin)
cdir = pwd; %% Cell Directory
%%%% get adjspikes object
if Args.Repetitions
    sp = adjspikes('Auto');
else
    sp = ispikes('Auto');
end
if isempty(sp)
    obj=createEmptyObject;
    fprintf('No AdjSpikes or Ispikes Objects')
    return
end
%%%%% move to the stiminfo session directory
stimInfo = stiminfo('Auto');
if isempty(stimInfo)
    obj=createEmptyObject;
    fprintf('No StimInfo Object')
    return
end
data.numSets = 1;
data.setNames{1} = cdir;
[data.firingRate,data.meancounts,data.spike_matrix,data.timebins, ...
    data.binsize] = CalcFiringRate(sp,stimInfo,Args);

% If max firing rate is below the histogram binsize, cancel the fit
if Args.Fit
    if max(data.firingRate)<Args.HistRate
        data.Fit.N=[];
        data.Fit.bins=[];
        data.Fit.Surrogate_Mean_R_Values=[];
        data.Fit.Surrogate_H_Value=[];
        data.Fit.Surrogate_P_Value=[];
        data.Fit.BootStrap_R_Values=[];
        data.Fit.BootStrap_Mean_R_Values=[];
        data.Fit.BootStrap_H_Value=[];
        data.Fit.BootStrap_P_Value=[];
        data.Fit.Firing_Rate_R_Values=[];
        data.Fit.FitName1='exp1';
        data.Fit.FitName2='power1';
    else
        data.Fit = fitfiringrate(data,Args,cdir);
    end
else
    data.Fit = [];
end
d.data = data;
number = 1;
holdaxis = get(sp,'HoldAxis');
nd = nptdata(number,holdaxis,cdir);
obj = class(d,'firingrate',nd);
if(Args.SaveLevels>0)
    fprintf('Saving firingrate object...\n');
    filename = 'firingrate.mat';
    fr = obj;
    save(filename,'fr')
end
