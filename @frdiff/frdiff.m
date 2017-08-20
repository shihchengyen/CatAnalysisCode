function obj = frdiff(varargin)
%@frdiff Constructor function for frdiff class
%   OBJ = frdiff('auto') attempts to create a frdiff object by ...

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0, ...
                'ClusterDirs',{''},'AdjSpikes',0,'BinSize',10);
Args.flags = {'Auto','AdjSpikes'};
[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'frdiff';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'fd';

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
if(dsize>1)
	% get current directory
	cwd = pwd;
	for idx = 1:dsize
		% change to cluster directory
		cd(Args.ClusterDirs{idx})
        if(~Args.AdjSpikes)
    		frobj = firingrate('auto','NoiseCorr',varargin{:});
        end
		if(idx==1)
            if(Args.AdjSpikes)
                % instantiate adjspikes object
                as1 = adjspikes('auto');
				% get rasters
				ras1 = as1.data.raster;
				% get number of repetitions
				reps = size(ras1,2);
				% get frame points
				framepts = as1.data.adjFramePoints;
				% get number of frames
				frames = (length(framepts) - 1) / reps;
				fpts = framepts(1):Args.BinSize:framepts(frames+1);
                % subtract first frame point since it might not be zero
                % while the rasters all have been adjusted to start at zero
				fp = fpts-fpts(1);
				sm = histcie(ras1,fp,'DropLast')' / Args.BinSize * 1000;
                fr = mean(sm)';
                nc = (sm - repmat(fr',reps,1))';
                noisecounts = nc(:);
                frobj.data.timebins = fp;
            else
				% get size and instantiate appropriate matrix
				fr = frobj.data.firingRate;
                noisecounts = frobj.data.meancounts;
                sm = frobj.data.spike_matrix;
            end
			[reps,nbins] = size(sm);
			% length of noise counts not the same as nbins since it is 
			% number of repetitions * nbins
			lnc = length(noisecounts);
			frmat = zeros(dsize,nbins);
			smat = cell(1,dsize);
            noisemat = zeros(dsize,lnc);
			% transpose firing rate so we can take difference in firing
			% rates for each bin easily
			frmat(1,:) = fr';
            noisemat(1,:) = noisecounts';
            smat{1} = sm;
            timebins = frobj.data.timebins;
		else
            if(Args.AdjSpikes)
                % instantiate adjspikes object
                as1 = adjspikes('auto');
				% get rasters
				ras1 = as1.data.raster;
				sm = histcie(ras1,fp,'DropLast')';
                fr = mean(sm);
                frmat(idx,:) = fr;
                nc = (sm - repmat(fr,reps,1))';
                noisemat(idx,:) = nc(:)';
                smat{idx} = sm;
            else
				frmat(idx,:) = frobj.data.firingRate';
                noisemat(idx,:) = frobj.data.meancounts';
                smat{idx} = frobj.data.spike_matrix;
            end
		end
		cd(cwd);
	end
	if(dsize==2)
		% take absolute difference in firing rates for each bin
		dfr = abs(diff(frmat));
		ndfr = 1;
		% get mean of differences normalized by number of frames
		frdmean = mean(dfr);
		% get max difference
		frdmax = max(dfr);
		[vmag,vang,vangp] = vecsimilarity(frmat',varargin{:});
		% compute the signal correlation coefficient between cells
		[ce1,cep1] = corrcoef(frmat');
		% get the non-redundant values
		ce = ce1(2,1);
		cep = cep1(2,1);
        % compute the Spearman correlation between cells
        [sc1,scp1] = spcorr(frmat');
        sc = sc1(2,1);
        scp = scp1(2,1);
		% compute the noise correlation between cells
		[nc1,ncp1] = corrcoef(noisemat');
        % get the non-redundant values
        nc = nc1(2,1);
        ncp = ncp1(2,1);
        % compute the noise correlations for individual windows
        % initialize with 0's and not nan's since corrcoef will return nan
        % if there is a problem
        winnc = repmat(nan,nbins,1);
        winncp = winnc;
        warning off MATLAB:divideByZero
        for bi = 1:nbins
        	[winnc1 winncp1] = ...
        		corrcoef(smat{1}(:,bi),smat{2}(:,bi));
        	winnc(bi) = winnc1(2,1);
        	winncp(bi) = winncp1(2,1);	
        end
        % convert smat to matrix
        smarray = cell2mat(smat);
        winstd = reshape(std(smarray),nbins,[]);
        warning on MATLAB:divideByZero
	elseif(dsize==3)
		% add first row to the end to make it easy to take the differences
		% between the vectors
		dfr1 = abs(diff([frmat; frmat(1,:)]));
		% rearrange the rows to make it more similar to the case when
		% dsize>3
		dfr = dfr1([1 3 2],:);
		ndfr = 3;
		% get the mean of the pair-wise differences
		dfrs1 = mean(dfr);
		% get the mean of the differences normalized by number of frames
		frdmean = mean(dfrs1);
		% get max difference
		frdmax = max(max(dfr));
		[vmag,vang,vangp] = vecsimilarity(frmat');
		% compute the correlation coefficient between cells
		[ce1,cep1] = corrcoef(frmat');
		% get the non-redundant values in the correlation coefficients
		ce2 = triu(ce1,1);
		% get the mean correlation coefficient
		ce = mean(ce2(:));
		% the max p-value is not really very informative but it is better 
		% than nothing
		cep = max(cep1(:));
		% compute the Spearman correlation between cells
		[sc1,scp1] = spcorr(frmat');
		% get the non-redundant values in the correlation coefficients
		sc2 = triu(sc1,1);
		% get the mean correlation coefficient
		sc = mean(sc2(:));
		% the max p-value is not really very informative but it is better 
		% than nothing
		scp = max(scp1(:));
        % initialize the noise correlation to empty matrices
        nc = [];
        ncp = [];
        winnc = [];
        winncp = [];
        winstd = [];
	else
		% allocate memory to store differences
        % the number of differences is the arithmetic sum of dsize-1
        % where the arithmetic sum is given by: n(first+last)/2
        n = dsize - 1;
		dfr = zeros(n*(1+n)/2,nbins);
		idx = 1;
		for jdx = 1:(dsize-1)
			for kdx = (jdx+1):dsize
				dfr(idx,:) = abs(frmat(kdx,:) - frmat(jdx,:));
				idx = idx + 1;
			end
		end
		ndfr = idx - 1;
		dfrs1 = mean(dfr);
		frdmean = mean(dfrs1);
		% get max difference
		frdmax = max(max(dfr));
		[vmag,vang,vangp] = vecsimilarity(frmat');
		% compute the correlation coefficient between cells
		[ce1,cep1] = corrcoef(frmat');
		% get the non-redundant values in the correlation coefficients
		ce2 = triu(ce1,1);
		% get the mean correlation coefficient
		ce = mean(ce2(:));
		% the max p-value is not really very informative but it is better than
		% nothing
		cep = max(cep1(:));
		% compute the Spearman correlation between cells
		[sc1,scp1] = spcorr(frmat');
		% get the non-redundant values in the correlation coefficients
		sc2 = triu(sc1,1);
		% get the mean correlation coefficient
		sc = mean(sc2(:));
		% the max p-value is not really very informative but it is better 
		% than nothing
		scp = max(scp1(:));
        % initialize the noise correlation to empty matrices
        nc = [];
        ncp = [];
        sc = [];
        scp = [];
        winnc = [];
        winncp = [];
        winstd = [];
	end
	% create object data
	% this is a valid object
	% these are fields that are useful for most objects
	data.numSets = 1;
	data.setNames{1} = cwd;
	
	data.frdiff = dfr;
	data.setIndex = [0; ndfr];
	data.frdmean = frdmean;
	data.vmag = vmag;
	data.vang = vang;
    data.vangp = vangp;
	data.frdmax = frdmax;
    % get mean firing rate
    frmean = mean(frmat);
    frmsqrt = sqrt(frmean);
    data.frdiffctrl = frmsqrt;
    data.frdctrlmean = mean(frmsqrt);
    data.frdctrlmax = max(frmsqrt);
    data.corrcoef = ce;
    data.corrpvalue = cep;
    data.spcorrcoef = sc;
    data.spcorrpvalue = scp;
    data.noisecorr = nc;
    data.noisecorrpvalue = ncp;
    data.winnoisecorr = winnc;
    data.winnoisecorrp = winncp;
    data.timebins = vecc(timebins);
    data.frmat = frmat';
    data.frSetIndex = [0; dsize];
    data.winstd = winstd;
	
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
	fprintf('Warning: %s requires at least 2 cells!\n', ...
		Args.classname);
	obj = createEmptyObject(Args);
end
	
function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;
data.setNames = '';

% these are object specific fields
data.frdiff = [];
data.setIndex = [];
data.frdmean = [];
data.vmag = [];
data.vang = [];
data.vangp = [];
data.frdmax = [];
data.frdiffctrl = [];
data.frdctrlmean = [];
data.frdctrlmax = [];
data.corrcoef = [];
data.corrpvalue = [];
data.spcorrcoef = [];
data.spcorrpvalue = [];
data.noisecorr = [];
data.noisecorrpvalue = [];
data.winnoisecorr = [];
data.winnoisecorrp = [];
data.timebins = [];
data.frmat = [];
data.frSetIndex = [];
data.winstd = [];

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
