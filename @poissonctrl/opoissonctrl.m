function obj = poissonctrl(varargin)
%@poissonctrl Constructor function for poissonctrl class
%   OBJ = poissonctrl('auto') attempts to create a poissonctrl object by ...

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0, ...
                'ClusterDirs',{''},'RNSize',1000,'Events',0,'VAng',0, ...
                'NoCorr',0,'Prctile',75,'PSparse',0,'Abs',0,'UseNCells',0);
Args.flags = {'Auto','Events','VAng','NoCorr','PSparse','Abs','UseNCells'};
[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'poissonctrl';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'pc';

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
		frobj = firingrate('auto',varargin{:});
		if(idx==1)
			% get size and instantiate appropriate matrix
			fr = frobj.data.firingRate;
            % noisecounts = frobj.data.meancounts;
            sm = frobj.data.spike_matrix;
			[reps,nbins] = size(sm);
			% length of noise counts not the same as nbins since it is 
			% number of repetitions * nbins
			% lnc = length(noisecounts);
			frmat = zeros(dsize,nbins);
			% smat = cell(1,dsize);
            % noisemat = zeros(dsize,lnc);
			% transpose firing rate so we can take difference in firing
			% rates for each bin easily
			frmat(1,:) = fr';
            % noisemat(1,:) = noisecounts';
            % smat{1} = sm;
            timebins = frobj.data.timebins;
		else
			frmat(idx,:) = frobj.data.firingRate';
            % noisemat(idx,:) = frobj.data.meancounts';
            % smat{idx} = frobj.data.spike_matrix;
		end
		cd(cwd);
	end
	% take the mean of all the cells
	mfr = mean(frmat);
	% compute the standard deviation which for the Poisson is just the 
	% square root of the mean
    % assume variance is 3 times the mean, i.e. Fano is 3
	sfr = sqrt(mfr*3);
	% generate randn matrix
	rnsize = Args.RNSize;
    % check option for psparse calculation
    if(Args.UseNCells)
        rns2 = rnsize * dsize;
		rmat = randn(rns2,nbins);
    else
        rns2 = rnsize * 2;
		rmat = randn(rns2,nbins);
        orows = 1:2:rns2;
        erows = orows + 1;
    end
	% scale random samples so they have the right mean and standard 
	% deviation
	rmat2 = repmat(mfr,rns2,1) + repmat(sfr,rns2,1) .* rmat;
    
    mr3 = nan;
    sr3 = nan;
    mv = nan;
    sv = nan;
    mjep = nan;
    sjep = nan;
    pszscore = nan;
    kwp = nan;
	if(Args.Events)
		rm2 = rmat2;
		% remove the zeros before finding the 75th percentile
		rm2(rm2==0) = nan;
		% find the 75th percentile of the non-zero firing rate for each surrogate
		thresh = prctile(rm2',Args.Prctile);
		% replicate threshold to create matrix to compare to rm2
		tmat = repmat(thresh',1,nbins);
		% find bins that are above threshold
		tbins = rm2 > tmat;
		% separate into odd and even rows and compute joint events
		jtbins = tbins(orows,:) & tbins(erows,:);
		% calculate event probabilities which is just the sum of each row in jtbins
		% dividied by the total number of bins
		jtprob = sum(jtbins,2)/nbins;
        mjep = mean(jtprob);
        sjep = std(jtprob);
		% take the difference along each row to compute event durations
		% jtd = diff(jtbins,1,2)';
		% find indices where the event started and ended
		% [surrn1,evstart] = find(jtd==1);
		% [surrn2,evend] = find(jtd==-1);
		% check for special cases where the firing rate at t=0 is above threshold
    elseif(Args.PSparse)
		rns1 = Args.RNSize + 1;
		% concatenate data and surrogates so we can do both calculations at the
		% same time
        rmat1 = zeros(rns2+dsize,nbins);
        d2 = 0;
        for didx = 1:dsize
            d1 = d2+1;
            d2 = didx*rnsize;
            rmat1(d1,:) = frmat(didx,:);
            rmat1((d1+1):(d2+1),:) = rmat(d1:d2,:);
        end
		% rmat1 = [frmat(1,:); rmat(1:rnsize,:); frmat(2,:); rmat(rns1:rns2,:)];
		% transpose 2002x850 matrix and then reshape into 850*1001x2
		rm2 = reshape(rmat1',[],dsize);
		% take the sum for each column
		if(Args.Abs)
            csum = sum(abs(rm2),2);
		else
            % this will actually make population more sparse since negative
            % values are set to 0
            rm2(rm2<0) = 0;
            csum = sum(rm2,2);
		end
		% the sum of the square of each column
		c2sum = sum(rm2.*rm2,2);
		warning off MATLAB:divideByZero
		ds1 = 1/dsize;
		dfr = (1 - ds1 * ((csum .* csum)./ c2sum))/(1-ds1);
		% reshape dfr to 850x1000 and then take transpose so we can call nanmean
		% which can only take data in columns
		dfr1 = reshape(dfr,[],rns1)';
		% get columns for surrogates
		sdata = dfr1(2:rns1,:);
		% compute mean and stdev for surrogates
		md = nanmean(sdata);
		sd = nanstd(sdata);
		% compute z-score for data
		pszscore = vecc((dfr1(1,:) - md) ./ sd);
        kwp = kruskalwallis(dfr1',num2str((1:1001)'),'off');
		warning on MATLAB:divideByZero
    end

    if(~Args.NoCorr || Args.VAng)
        r3 = zeros(rnsize,1);
        v = r3;
        tempfrmat = zeros(2,nbins);
        warning off MATLAB:divideByZero
        for idx = 1:rnsize
            fr1 = rmat2(orows(idx),:);
            fr2 = rmat2(erows(idx),:);
            if(~Args.NoCorr)
				tempr = corrcoef(fr1,fr2);
				r3(idx) = tempr(1,2);
			end
			if(Args.VAng)
				tempfrmat(1,:) = fr1;
				tempfrmat(2,:) = fr2;
				[dump,v(idx)] = vecsimilarity(tempfrmat,varargin{:});
			end
        end
        warning on MATLAB:divideByZero
		% compute the corrcoef
		% r = corrcoef(rmat2');
		% grab the correlation coefficients
		% r2 = triu(r,1);
		% r3 = r2(r2>0);
		% get the mean and std for the r values
        if(~Args.NoCorr)
			mr3 = mean(r3);
			sr3 = std(r3);
        end
        if(Args.VAng)
			% get the mean and std for the v values
			mv = mean(v);
			sv = std(v);
        end		
    end
    
	% create object data
	% this is a valid object
	% these are fields that are useful for most objects
	data.numSets = 1;
	data.setNames{1} = cwd;
	
    data.ncells = dsize;
	data.meancorr = mr3;
    data.stdcorr = sr3;
    data.meanvang = mv;
    data.stdvang = sv;
    data.meanjep = mjep;
    data.stdjep = sjep;
    data.pszscore = pszscore;
    data.kwp = kwp;
    data.mfr = mfr';
	
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
data.ncells = 0;
data.meancorr = [];
data.stdcorr = [];
data.meanvang = [];
data.stdvang = [];
data.meanjep = [];
data.stdjep = [];
data.pszscore = [];
data.kwp = [];
data.mfr = [];

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
