function obj = poissonctrl(varargin)
%@poissonctrl Constructor function for poissonctrl class
%   OBJ = poissonctrl('auto') attempts to create a poissonctrl object by ...

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0, ...
                'ClusterDirs',{''},'RNSize',1000,'Events',0,'VAng',0, ...
                'NoCorr',0,'Prctile',75,'PSparse',0,'Abs',0,'UseNCells',0, ...
                'StdMultiple',[],'RepsVar',0,'CellsVar',0,'FrameVar',0, ...
                'FanoMax',0,'SumVar',0,'StdMinFR',[],'CellPopVar',0, ...
                'FRPrctile',75,);
Args.flags = {'Auto','Events','VAng','NoCorr','PSparse','Abs','UseNCells', ...
        'RepsVar','CellsVar','FrameVar','FanoMax','SumVar','CellPopVar'};
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
		sm = frobj.data.spike_matrix;
        minfr = Args.StdMinFR;
        warning off MATLAB:divideByZero
		if(idx==1)
			% get size and instantiate appropriate matrix
			% fr = frobj.data.firingRate;
            fr = mean(sm);
            % noisecounts = frobj.data.meancounts;
			[reps,nbins] = size(sm);
            % set up indices so we don't have to do this inside the loop
            repindices = 1:reps;
			% length of noise counts not the same as nbins since it is 
			% number of repetitions * nbins
			% lnc = length(noisecounts);
			frmat = zeros(dsize,nbins);
            % set up array to store reps of all cells
		    smat = zeros(reps*dsize,nbins);
            smat(repindices,:) = sm;
            if(Args.CellsVar)
                % set up vector to store fano for each cell
                cellfano = zeros(dsize,1);
                % compute variability for each cell
                repspikecounts = sum(sm,2);
                cellfano(1) = (std(repspikecounts)^2) / mean(repspikecounts);
            elseif(Args.FrameVar)
                % set up vector to compute average fano for each frame
                if(Args.SumVar)
                    avgfano = std(sm).^2;
                    if(~isempty(minfr))
                        % set the variance of bins in which the mean is smaller
                        % than 1 to zero
                        avgfano(fr<minfr) = 0;
                    end
                elseif(Args.CellPopVar)
                    avgfano = zeros(dsize,nbins);
                    avgfano(1,:) = std(sm).^2 ./ fr;
                else
                    % avgfano = std(sm).^2 ./ mean(sm);  
                    avgfano = std(sm).^2 ./ fr;  
                end
                % replace nan's with zero
                avgfano(isnan(avgfano)) = 0;
            end
            % noisemat = zeros(dsize,lnc);
			% transpose firing rate so we can take difference in firing
			% rates for each bin easily
			% frmat(1,:) = fr';
            frmat(1,:) = fr;
            % noisemat(1,:) = noisecounts';
            % timebins = frobj.data.timebins;
		else % if(idx==1)
            smat((idx-1)*reps+repindices,:) = sm;
            fr = mean(sm);
            if(Args.CellsVar)
                % compute variability for each cell
                repspikecounts = sum(sm,2);
                cellfano(idx) = (std(repspikecounts)^2) / mean(repspikecounts);
            elseif(Args.FrameVar)
                % set up vector to compute average fano for each frame
                % tfano = std(sm).^2 ./ mean(sm);
                tfano = std(sm).^2 ./ fr;
                % replace nan's with zero
                tfano(isnan(tfano)) = 0;
                if(Args.FanoMax)
                    avgfano = max([avgfano; tfano]);
                elseif(Args.SumVar)
                    tmpstd = std(sm).^2;
                    if(~isempty(minfr))
                        % set the variance of bins in which the mean is smaller
                        % than 1 to zero
                        tmpstd(fr<minfr) = 0;
                    end
                    avgfano = avgfano + tmpstd;
                elseif(Args.CellPopVar)
                    avgfano(idx,:) = std(sm).^2 ./ fr;
                else
                    avgfano = avgfano + tfano;
                end
            end
			% frmat(idx,:) = frobj.data.firingRate';
            frmat(idx,:) = fr;
            % noisemat(idx,:) = frobj.data.meancounts';
            % smat{idx} = frobj.data.spike_matrix;
		end
        warning on MATLAB:divideByZero
		cd(cwd);
	end
    if(Args.RepsVar || Args.CellPopVar)
        % take the mean of all the reps
        mfr = mean(smat);
        % compute the standard deviation for all the reps
        sfr = std(smat);
    else    
		% take the mean of all the cells
		mfr = mean(frmat);
		% compute the standard deviation which for the Poisson is just the 
		% square root of the mean
        % get stdmultiple
        stdmul = Args.StdMultiple;
		if(~isempty(stdmul))
            sfr = sqrt(mfr*stdmul);
        elseif(Args.CellsVar)
            sfr = sqrt(mfr*mean(cellfano));
        elseif(Args.FrameVar)
            if(Args.FanoMax)
                sfr = sqrt(mfr .* avgfano);
            elseif(Args.SumVar)
                sfr = sqrt(avgfano);
            else
                sfr = sqrt(mfr .* avgfano / dsize);
            end
        else
            sfr = sqrt(mfr);
        end
    end
    
    mr3 = nan;
    lqr3 = nan;
    uqr3 = nan;
    mv = nan;
    lqv = nan;
    uqv = nan;
    mjep = nan;
    lqjep = nan;
    uqjep = nan;
    pszscore = nan;
    kwp = nan;        

    if(Args.CellPopVar)
        % compare cell fanos to population fanos
        warning off MATLAB:divideByZero
        popfano = sfr.^2 ./ mfr;
        reps1 = 1:reps;
        reps2 = (reps+1):(2*reps);
        s1 = smat(reps1,:);
        s2 = smat(reps2,:);
		m1 = mean(s1);
		m2 = mean(s2);
		m1p = m1;
		m2p = m2;
		% replace 0's with nan
		m1p(m1==0) = nan;
		m2p(m2==0) = nan;
        frprc = Args.FRPrctile;
		p1 = prctile(m1p,frprc);
		p2 = prctile(m2p,frprc);
		rbins2 = (m1>p1) & (m2>p1);
		jbins = sum(rbins2);
		mfano2 = min(avgfano);
        rb2mfano = mfano2(rbins2)';
        rb2pfano = popfano(rbins2)';
		vbins = sum(rb2mfano<rb2pfano);
        
        [r,p] = corrcoef(smat');
        c1 = r(reps1,reps1);
        c2 = r(reps2,reps2);
        c3 = r(reps2,reps1);
        ctrilind = logical(tril(ones(reps,reps),-1));
        c1v = c1(ctrilind);
        c2v = c2(ctrilind);
        [kwp,kwt,kws] = kruskalwallis(concat(c1v,c2v,c3(:),'Columnwise'),{'1','2','3'},'off');
        c = multcompare(kws,0.05,'off');
        % check if lower confidence limit of col1 - col3 is greater than 0
        if(c(2,3)>0)
            % cross corr is lower than autocorr
            ccac = 1;
        else
            ccac = 0;
        end
        % check if confidence limit of col2 - col3 is greater than 0
        if(c(3,3)>0)
            ccac = ccac + 1;
        end
        mr3 = vbins/jbins;
        warning on MATLAB:divideByZero
    else
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
        
        % vector for computing median and lower and upper quartiles
        prctv = [50 25 75];
	
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
            prc = prctile(jtprob,prctv);
            mjep = prc(1);
            lqjep = prc(2);
            uqjep = prc(3);
	%         mjep = mean(jtprob);
	%         sjep = std(jtprob);
			% take the difference along each row to compute event durations
			% jtd = diff(jtbins,1,2)';
			% find indices where the event started and ended
			% [surrn1,evstart] = find(jtd==1);
			% [surrn2,evend] = find(jtd==-1);
			% check for special cases where the firing rate at t=0 is above threshold
        elseif(Args.PSparse)
			rns1 = rnsize + 1;
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
			% md = nanmean(sdata);
			% sd = nanstd(sdata);
            prc = prctile(sdata,prctv);
            mjep = prc(1);
            lqjep = prc(2);
            uqjep = prc(3);
			% compute z-score for data
			% pszscore = vecc((dfr1(1,:) - md) ./ sd);
            % find percentage of surrogates with lower psparsity
            datamat = repmat(dfr1(1,:),rnsize,1);
            pszscore = (sum(datamat>sdata)/rnsize)';
            kwp = kruskalwallis(dfr1',num2str((1:(rnsize+1))'),'off');
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
					[dump,v(idx)] = vecsimilarity(tempfrmat',varargin{:});
				end
            end
            warning on MATLAB:divideByZero
			% compute the corrcoef
			% r = corrcoef(rmat2');
			% grab the correlation coefficients
			% r2 = triu(r,1);
			% r3 = r2(r2>0);
			% get the median and quartiles for the r values
            if(~Args.NoCorr)
                p3 = prctile(r3,prctv);
				mr3 = p3(1);
                lqr3 = p3(2);
                uqr3 = p3(3);
            end
            if(Args.VAng)
				% get the median and quartiles for the v values
                p3 = prctile(v,prctv);
				mv = p3(1);
				lqv = p3(2);
                uqv = p3(3);
            end		
        end
    end % if(Args.CellPopvar) 
    
	% create object data
	% this is a valid object
	% these are fields that are useful for most objects
	data.numSets = 1;
	data.setNames{1} = cwd;
	
    data.ncells = dsize;
	data.mediancorr = mr3;
    data.lqcorr = lqr3;
    data.uqcorr = uqr3;
    data.medianvang = mv;
    data.lqvang = lqv;
    data.uqvang = uqv;
    data.medianjep = mjep;
    data.lqjep = lqjep;
    data.uqjep = uqjep;
    data.pszscore = pszscore;
    data.kwp = kwp;
    data.mfr = [mfr' sfr'];
	
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
data.mediancorr = [];
data.lqcorr = [];
data.uqcorr = [];
data.medianvang = [];
data.lqvang = [];
data.uqvang = [];
data.medianjep = [];
data.lqjep = [];
data.uqjep = [];
data.pszscore = [];
data.kwp = [];
data.mfr = [];

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
