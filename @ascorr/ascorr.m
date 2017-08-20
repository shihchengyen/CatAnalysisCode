function obj = ascorr(varargin)
%@ascorr Constructor function for ascorr class
%   OBJ = ascorr('auto') attempts to create a ascorr object by ...

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0, ...
                'ClusterDirs',{''},'iSpikes',0,'WindowSize','frame', ...
                'NoCorrCoef',0,'RepFrames',[],'CheckStd',0);
Args.flags = {'Auto','iSpikes','NoCorrCoef','CheckStd'};
[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'ascorr';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'ac';

numArgin = nargin;
if(numArgin==0)
	% create empty object
	obj = createEmptyObject(Args);
elseif( (numArgin==1) & isa(varargin{1},Args.classname))
	obj = varargin{1};
else
	% create object using arguments
	if(Args.Auto)
        % change to the proper directory
        % [pdir,cdir] = getDataDirs('session','relative','CDNow');
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
			obj = createObject(Args,modvarargin{:});
		end
        % change back to previous directory if necessary
        % if(isempty(cdir))
        %     cd(cdir)
        % end
	end
end

function obj = createObject(Args,varargin)

% if ClusterDirs is empty, try to generate it from the current directory
if(isempty(Args.ClusterDirs))
    Args.ClusterDirs = getDataDirs('GetClusterDirs');
end
% get number of directories
dsize = size(Args.ClusterDirs,2);
if(dsize==2)
	% get current directory
	cwd = pwd;
	% change to cluster directory
	cd(Args.ClusterDirs{1})
    if(Args.iSpikes)
        isp1 = ispikes('auto',varargin{:});
        % change to second cluster direcctory
		cd(Args.ClusterDirs{2})
		% get adjspikes for second cell
		isp2 = ispikes('auto',varargin{:});
		% return to original directory
		cd(cwd);

        cluster1 = isp1.data.trial.cluster;
        cluster2 = isp2.data.trial.cluster;
        % get max spike time
        % maxsptime = max([cluster1.spikes(cluster1.spikecount) cluster2.spikes(cluster2.spikecount)]);
        % get stiminfo
        st = stiminfo('auto');
        % convert framePoints from data points to time
        % need to subtract 1 since framePoints starts from 1 which is
        % actually time zero
        framePoints = (st.data.framePoints-1)/st.data.catInfo.samplingrate*1000;
        % get size of framePoints
        fpl = length(framePoints)-1;
        wsize = Args.WindowSize;
        if(ischar(wsize) && strcmp(wsize,'frame'))
            % set up historgram bins using framePoints
            histbins = framePoints(1:fpl);
        else
            % set up histogram bins using 1st and second-last framePoints
            % the last point in framePoints is always 0 for some reason
            histbins = framePoints(1):Args.WindowSize:framePoints(fpl);
            fpl = length(histbins);
        end
        hcounts = histcie(concat(isp1.data.trial.cluster.spikes,isp2.data.trial.cluster.spikes)',histbins,'DropLast');

        if(Args.CheckStd)
            % code to check out a couple of measures for determining
            % responsiveness of cells to m-sequence stimuli
            % get std of spike counts
            coefs = std(hcounts);
            % get percentage of frames that have non-zero spike counts
            contcoef = (sum(hcounts>0)/fpl)';
            % get total number of spikes 
			vang = sum(hcounts);
            % get average firing rate
			contvang = vang/(fpl*st.data.catInfo.frame_duration/1000);
        else        
            repFrames = Args.RepFrames;
            if(~isempty(repFrames))
                % figure out how many frames we can fit into a integer multiple
                % of repFrames
                nreps = floor(fpl/repFrames);
                frameidx = 1:(nreps * repFrames);
                % reshape to Args.RepFrames
                hc1 = reshape(hcounts(frameidx,1),repFrames,[]);
                hc2 = reshape(hcounts(frameidx,2),repFrames,[]);
                coefs = zeros(nreps,1);
                for repi = 1:nreps
                    if(Args.NoCorrCoef)
                        % get spike correlation without subtracting the mean before
                        % multiplying the spike counts
                        hctmp = [hc1(:,repi) hc2(:,repi)];
                        coefs(repi) = (mean(prod(hctmp,2))-(prod(mean(hctmp))))/prod(std(hctmp));
                    else
                        r = corrcoef(hc1(:,repi),hc2(:,repi));
                        coefs(repi) = r(1,2);
                    end
                end
            else % if(~isempty(repFrames))
        		coefs = [];
            end
            if(Args.NoCorrCoef)
                % get spike correlation without subtracting the mean before
                % multiplying the spike counts
                contcoef = (mean(prod(hcounts,2))-(prod(mean(hcounts))))/prod(std(hcounts));
            else
                r = corrcoef(hcounts);
        		contcoef = r(1,2);
            end
			vang = [];
			contvang = [];
        end
    else
    	as1 = adjspikes('auto',varargin{:});
		% get rasters
		ras1 = as1.data.raster;
		% get number of repetitions
		reps = size(ras1,2);
		% get frame points
		framepts = as1.data.adjFramePoints;
		% get number of frames
		frames = (length(framepts) - 1) / reps;
		fpts = framepts(1:(frames+1));
        % subtract first frame point since it might not be zero
        % while the rasters all have been adjusted to start at zero
		fp = fpts-fpts(1);
		repcounts1 = histcie(ras1,fp,'DropLast');
		
		% change to second cluster directory
		cd(Args.ClusterDirs{2})
		% get adjspikes for second cell
		as2 = adjspikes('auto',varargin{:});
		repcounts2 = histcie(as2.data.raster,fp,'DropLast');
		% return to original directory
		cd(cwd);
		
		% now compute corrcoef for each rep
		% first allocate memory
		coefs = zeros(reps,1);
        vang = coefs;
		for repi = 1:reps
            rc1 = repcounts1(:,repi);
            rc2 = repcounts2(:,repi);
			r = corrcoef(rc1,rc2);
			coefs(repi) = r(1,2);
			[vmag,vang(repi)] = vecsimilarity([rc1 rc2],varargin{:});
		end
	
        % now compute corrcoef for all reps
        r1 = repcounts1(:);
        r2 = repcounts2(:);
        r = corrcoef(r1,r2);
        contcoef = r(1,2);
        [vmag,contvang] = vecsimilarity([r1 r2],varargin{:});
    end
        
	% create object data
	% this is a valid object
	% these are fields that are useful for most objects
	data.numSets = 1;
	data.setNames{1} = cwd;
	data.ascorr = coefs;
    data.contcoefs = contcoef;
    data.vang = vang;
    data.contvang = contvang;
        
	% create nptdata so we can inherit from it
	n = nptdata(data.numSets,0,cwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	if(Args.SaveLevels)
		fprintf('Saving %s object...\n',Args.classname);
		eval([Args.matvarname ' = obj;']);
		% save object
		eval(['save ' Args.matname ' ' Args.matvarname]);
	end	
else
	fprintf('Warning: %s requires 2 cells!\n', ...
		Args.classname);
	obj = createEmptyObject(Args);
end
	
function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;
data.setNames = '';

% these are object specific fields
data.ascorr = [];
data.contcoefs = [];
data.vang = [];
data.contvang = [];

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
