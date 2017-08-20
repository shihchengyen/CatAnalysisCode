function obj = catcell(varargin)
%catcell Object used to store data structure from cat experiments.
%   OBJ = catcell('auto',varargin) attempts to load a mat file with
%   name '*g*c*.mat' to create a CATCELL object.

Args = struct('RedoLevels',0,'SaveLevels',0, ...
			  'Auto',0,'CellName','','SubBinSize',1, ...
              'MatSearchString','*g*c*.mat');
			  
[Args,varargin] = getOptArgs(varargin,Args, ...
	'flags',{'Auto'}, ...
	'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
	'subtract',{'RedoLevels','SaveLevels'});

if nargin==0
	% create empty object
	obj = createEmptyObject;
elseif( (nargin==1) & (isa(varargin{1},'catcell')) )
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
		% check for saved object
		if(ispresent('catcell.mat','file','CaseInsensitive') & (Args.RedoLevels==0))
			fprintf('Loading saved catcell object...\n');
			% load saved object and exit
			l = load('catcell.mat');
			obj = l.cc;
		else
			% no saved object so try to create one
			% try to find name of .mat file
			a = nptDir(Args.MatSearchString);
			if(isempty(a))
				% create empty object
				obj = createEmptyObject;
			else
				mat = load(a(1).name);
                % hack to work with current data files
                if(strcmp(Args.MatSearchString,'data.mat'))
                    % get cell name based on directory path
                    cwd = pwd;
                    % get cell name from last part of path
                    [p,cname] = nptFileParts(cwd);
                    % get group name from next part of path
                    [p,gname] = nptFileParts(p);
                    % get session name from next part of path
                    [p,sname] = nptFileParts(p);
                    % get site name from next part of path
                    [p,tname] = nptFileParts(p);
                    % get cat name from next part of path
                    [p,aname] = nptFileParts(p);
                    % reconstruct to get cell name
                    Args.CellName = [aname tname sname gname cname];
                else
					[fpath,fname] = nptFileParts(a(1).name);
					Args.CellName = fname;
				end
				obj = createObject(mat.data,Args);
			end
		end
	else
		% hopefully we were given a data structure
		% create object using data structure
		obj = createObject(varargin{1},Args);
	end
end

function obj = createEmptyObject

n = nptdata(0,0);
d.data.cellname = '';
d.data.scmean = [];
d.data.scstd = [];
obj = class(d,'catcell',n);

function obj = createObject(data,Args)

data.cellname = Args.CellName;
sc = histcie(data.cell_info.adjusted_spiketrain, ...
	data.stimulus_info.adjusted_frames_vector,'DropLast');
% compute number of frames presented
nframes = data.stimulus_info.end_frame - data.stimulus_info.start_frame + 1;
% reshape into repetitions and then transpose so we can take the mean and std easily
sc1 = reshape(sc,nframes,[])';
data.scmean = mean(sc1);
data.scstd = std(sc1);

data.latencySubBinSize = Args.SubBinSize;
% divide frame vector into requested size
[subbins,data.latencyActualSubBinSize,nsbins] = divideBins( ...
	data.stimulus_info.adjusted_frames_vector,'SubBinSize',Args.SubBinSize);
% take histogram using requested bin size
binsc = histcie(data.cell_info.adjusted_spiketrain,subbins,'DropLast');
% reshape into nsbins x (frames*reps)
binsc1 = reshape(binsc,nsbins,[]);
% take the sum of each bin and append 0 at the end so we can use the
% 'histc' option bar
data.latencyBinSums = [sum(binsc1,2); 0];
n = nptdata(1,0);
d.data = data;
obj = class(d,'catcell',n);
if(Args.SaveLevels>0)
	fprintf('Saving catcell object...\n');
	cc = obj;
	save catcell cc
end
