function [data,surrmean,surrstd,surrhist,zscore] = computeSurrData(sdatafile)

% set constants
DATA_COL = 1;
SHIFT_COL = 2;
SURR_COL = 3;
% the kstest2 is not acurate if n1*n2/(n1+n2) < 4. Since
% the number of surrogate spikes is always going to be 1000
% times the number of real spikes, this simplifies to
% 1000n^2/1001n < 4, which further simplies to 
% n < 4001/1000, which is pretty close to 4
KSNSPIKES = 4;

% load the surrogate data file
load(sdatafile);
% load the data file
load(cellpairdatafile);

% get list of surrogate files
surrlist = nptDir([surrprefix '*.mat']);
% check to make sure we have the complete set of files
surrnum = length(surrlist);
if(surrnum<numwindows)
	error('Incomplete number of surrogate files!');
else
	% allocate memory for output arguments
	data = zeros(numwindows,1);
	shiftp = data;
	sdata = zeros(numSurrogates,numwindows);
	% create vector for summing appropriate bins in xcorr vector
	% create vector using windowlength
	sumvec = zeros(1,lagslength);
    % allocate memory for storing mean and std of surrogate xcorr
    smstmp = sumvec';
	% find the center of lagslength
	llcenter = ceil(lagslength/2);
	centerbins = (NumCentralBins-1)/2;
	% create ones in places that correspond to NumCentralBins
	sumvec((llcenter-centerbins):(llcenter+centerbins)) = 1;
	% generate indices for the surrogates so we don't have to keep
	% generating it inside the for loop
	numsurrvec = 1:numSurrogates;
	srind = numsurrvec + 2;
	srind2 = srind - 1;
	% create memory for saving xcorr data
	dxcorr = zeros(lagslength,numwindows);
	dxcz = repmat(nan,lagslength,numwindows);
	dxcmean = dxcz;
	dxcstd = dxcz;
	% create matrix to save result of testing psth distributions
	hpval = repmat(nan,2,numwindows);
	% generate vector of lagslength
	lagsvec = 1:lagslength;
	% generate matrix for result of jbtest
	jbmat = repmat(nan,lagslength,numwindows);
	% generate matrix for result of lillietest
	limat = jbmat;
	
    % get vector that picks out just the data without the shuffle edges
    sadd = shuffle + 1;
    % get the time of the subbins
    bintimes1 = subbins;
    % get length of bintimes1
    bt1l = length(bintimes1);
    % get binsize so we can add the necessary bins for the shuffles that 
    % go off the edge
    bsize = bintimes1(2)-bintimes1(1);
    % get extra bins
    extrabins = (1:shuffle)'*bsize;
    bintimes = [-flipud(extrabins); bintimes1; bintimes1(bt1l) + extrabins];
    dvec = sadd:(windowlength+shuffle);
    pbins = windowlength+(2*shuffle);
    psthbins = repmat(nan,pbins,numwindows);
    % set second index first so we don't have to change memory size
    psths{2} = psthbins;
    psths{1} = psthbins;
    
	ns1 = numSurrogates + 1;
	ns2 = numSurrogates + 2;
	ws1 = psthwin - 1;
	% get half window size
	halfwinsize = ws1 / 2;
	% get length of central window
	lcentralwin = windowlength;
	% compute length of psthbins and save for use later
	lpsthbins = lcentralwin + ws1;
    % columns representing real data and surrogates for first cell
	ps1ind = [numsurrvec ns1];
    % columns representing real data and surrogates for second cell
	ps2ind = ps1ind + ns1;
    % columns representing real data and surrogates for first cell in the
    % preceding window
	ps1ind2 = ps2ind + ns1;
    % columns representing real data and surrogates for first cell in the
    % subsequent window
	ps1ind3 = ps1ind2 + ns1;
    % columns representing real data and surrogates for second cell in the
    % preceding window
	ps2ind2 = ps1ind3 + ns1;
    % columns representing real data and surrogates for second cell in the
    % subsequent window
	ps2ind3 = ps2ind2 + ns1;
	% get indices representing window
	cwinidx = 1:lcentralwin;
	cwinidx2 = cwinidx+lcentralwin;
	% compute the relevant indices
	didx = cwinidx + ws1;
	surridx1 = 2:ns1;
	surridx2 = surridx1 + ns1;
	smoothmtx = convmtx(ones(psthwin,1),lpsthbins);
    defaultspt = repmat(nan,1,ns1);
	dspsthz = zeros(numwindows,2*lcentralwin);
    
    % pdf for truncated normal distribution
    % taken from example from Mathworks:
    % http://www.mathworks.com/products/statistics/demos.html?file=/products/demos/shipping/stats/customdist2demo.html
    % but since we have values of 0, we normalize by the cummulative
    % pdf for the value of -1
    % pdf_truncnorm = @(x,mu,sigma) normpdf(x,mu,sigma) ./ (1-normcdf(-1,mu,sigma));

    % turn off warning
	warning off MATLAB:divideByZero
	for sidx = 1:numwindows
		l = load([surrprefix num2str(sidx,'%04d')]);
	    % if there were no spikes or if there were no reps in which there
	    % were spikes in both cells, skip this analysis since there would
	    % no xcorr values to look at.
	    if(~isempty(l.sptrains))
	    	% *** Compute Synchronous Spikes ***
			% sum counts within maxlags for data and surrogates by creating 
			% appropriate vector
			nsyncspikes = sumvec * l.xc;
			% separate the surrogates from the data
			data(sidx) = nsyncspikes(DATA_COL);
			shiftp(sidx) = nsyncspikes(SHIFT_COL);
			sdata(:,sidx) = nsyncspikes(srind)';
			
			% *** Compute XCORR Z-Scores and test for Normality ***
			% save the xcorr data so we don't have to compute it in the plot 
			% function
			xcdata = l.xc(:,DATA_COL);
			% get surrogate data
			xcsdata = l.xc(:,srind)';
            % initialize smean and sstd to zeros
            smean = smstmp;
            sstd = smstmp;
			% compute the jbtest to test for presence of normal distributions
			for lagi = lagsvec
				testdata = full(xcsdata(:,lagi));
				[dummy,jbmat(lagi,sidx)] = jbtest(testdata);
				[dummy,limat(lagi,sidx)] = lillietest(testdata);
%                 if(sum(testdata)>0)
%                     start = [mean(testdata) std(testdata)];
%                     try
%                         [phat,pci] = mle(testdata,'pdf',pdf_truncnorm,'start',start,'lower',[-Inf 0]);
%                         % fprintf('%d mean: %f std: %f\n',lagi,phat(1),phat(2));
%                         smean(lagi) = phat(1);
%                         sstd(lagi) = phat(2);
%                     catch
%                         fprintf('%d %s\n',lagi,lasterr);
%                         % phat = nan;
%                     end
%                 %else
%                 %    fprintf('%d Skipped\n',lagi);
%                 end
			end
			% compute the mean and std of the surrogates
			smean = mean(xcsdata)';
			sstd = std(xcsdata)';
			% convert data to z-score using the mean and std of the surrogates
			dxcz(:,sidx) = (xcdata-smean)./sstd;
			dxcorr(:,sidx) = xcdata;
			dxcmean(:,sidx) = smean;
			dxcstd(:,sidx) = sstd;	 
			
			% *** Test smoothed PSTH for differences between data and surrogates ***
			% load frame before and frame after, which is 2 windows before 
			% and after since we are stepping in half frame increments
			% convert cell arrays to matrices
			% store in cell array since there is a for-loop later which
			% will make use of that data
            if(~isempty(l.sptrains{2}))
    			spt{2} = full(cell2mat(l.sptrains{2}));
            else
                spt{2} = defaultspt;
            end
            if(~isempty(l.sptrains{1}))
    			spt{1} = full(cell2mat(l.sptrains{1}));
            else
                spt{1} = defaultspt;
            end
			l0spt1 = zeros(1,ns1);
			l0spt2 = l0spt1;
			l1spt1 = l0spt1;
			l1spt2 = l0spt1;
			if(sidx>2)
				l0 = load([surrprefix num2str(sidx-2,'%04d')]);
                if(~isempty(l0.sptrains))
                    if(~isempty(l0.sptrains{1}))
						l0spt1 = full(cell2mat(l0.sptrains{1}));
                        % l0sp1 = l0spt1(:,1);
                    end
                    if(~isempty(l0.sptrains{2}))
						l0spt2 = full(cell2mat(l0.sptrains{2}));
                        % l0sp2 = l0spt2(:,1);
                    end
                end
			end
			if(sidx<(numwindows-1))
				l1 = load([surrprefix num2str(sidx+2,'%04d')]);
                if(~isempty(l1.sptrains))
                    if(~isempty(l1.sptrains{1}))
						l1spt1 = full(cell2mat(l1.sptrains{1}));
                        % l1sp1 = l1spt1(:,1);
                    end
                    if(~isempty(l1.sptrains{2}))
						l1spt2 = full(cell2mat(l1.sptrains{2}));
                        % l1sp2 = l1spt2(:,1);		
                    end
                end
			end
			% concatenate data together in order to perform histogram just
			% once
			smat = concat(spt{1},spt{2},l0spt1,l1spt1, ...
				l0spt2,l1spt2,'Columnwise');
			% shift by 0.5 at the beginning and the end so we can use
			% histcie instead of hist as hist includes all values outside
			% the range in the first and last bins
            spsthbins = (startbins(sidx)-halfwinsize):(endbins(sidx)+halfwinsize);
			psthhistcbins = [spsthbins-0.5 spsthbins(lpsthbins)+0.5];
			psth = histcie(smat,psthhistcbins,'DropLast','DataCols');
			% add data from frame before and frame after before computing
			% smoothed PSTH
			psth1 = psth(:,ps1ind) + psth(:,ps1ind2) + psth(:,ps1ind3);
			psth2 = psth(:,ps2ind) + psth(:,ps2ind2) + psth(:,ps2ind3);
			% now compute the sliding window smoothed psth
			spsth = smoothmtx * [psth1 psth2];
            % grab the data
            dspsth1 = spsth(didx,1);
            dspsth2 = spsth(didx,ns2);
            % grab the surrogates
            sspsthmat = [spsth(didx,surridx1)' spsth(didx,surridx2)'];
            % compute mean and std of the surrogates
            sspm = mean(sspsthmat);
            sspstd = std(sspsthmat);
            % compute z-score for the data
            dspsthz(sidx,:) = ([dspsth1' dspsth2'] - sspm) ./ sspstd;

			loopbins = ( (startbins(sidx)-shuffle):(endbins(sidx)+shuffle) )';
			for cidx = cellvec
				% get the spike times of the first cell
				% call the full function to make sure any sparse matrices 
				% gets converted to full matrices otherwise kstest2 will
				% not work. Sparse matrices are possible now since binidx 
				% is saved as a sparse matrix in shufflesync and the code
				% in shufflesyncsurr will index into binidx and save it 
				% directly to l.sptrains.
                sptrains = spt{cidx};
				if(~isempty(sptrains))
					% get the spike times of the data
					sptimes = sptrains(:,DATA_COL);
					% the kstest2 is not acurate if n1*n2/(n1+n2) < 4. Since
					% the number of surrogate spikes is always going to be 1000
					% times the number of real spikes, this simplifies to
					% 1000n^2/1001n < 4, which further simplies to 
					% n < 4001/1000, which is pretty close to 4
					if(length(sptimes)>KSNSPIKES)
						% get the spike times of the surrogates
						surrtimes = sptrains(:,srind2);
						srtimes = surrtimes(:);
						% do test for difference
						[dummy,hpval(cidx,sidx)] = kstest2(sptimes,srtimes);
					end
	
					% check to see if we need to pad sptrains so the histogram is
					% taken in the right dimension
					if(size(sptrains,1)==1)
						sptrains = concatenate(sptrains,nan);
					end
					loopcounts = hist(sptrains,loopbins);
					psths{cidx}(:,sidx) = mean(loopcounts(:,surrvec1),2);
				end
			end % for cidx = cellvec
			psthbins(:,sidx) = loopbins;
		end % if(~isempty(l.sptrains))
	end % for sidx = 1:numwindows
	
	% compute histogram of surrogates
	surrmean = mean(sdata)';
	surrstd = std(sdata)';
	% set bins so it covers the max of the data as well as of the surrogates
	maxn = max([sdata(:); data(:); shiftp(:)]);
	bins = SurrHistMin:SurrHistStep:(ceil(maxn/SurrHistStep) ...
		* SurrHistStep);
	% do hist since values are integers
	surrhist = hist(sdata,bins);
	% find bins with std not equal to 0
	si = find(surrstd);
	% create z-score vector consisting of nans
	dzscore = repmat(nan,numwindows,1);
	% fill bins that don't have std of 0
	dzscore(si) = (data(si) - surrmean(si)) ./ surrstd(si);
	% convert matrices with no nan's to sparse matrices to save space
	data = sparse(data);
	shiftp = sparse(shiftp);
	surrmean = sparse(surrmean);
	surrstd = sparse(surrstd);
	surrhist = sparse(surrhist);
	dxcorr = sparse(dxcorr);
	
	% calculate the synchronous spikes for the raster display
	repvec = 1:reps;
	sprepvec = circshift(repvec,[1 -1]);
	repvec2 = repvec + reps;
	sprepvec2 = sprepvec + reps;
	nspikes = binidxsize(1);
	mat1 = repmat(-1,nspikes,2);
	mat2 = ones(2,nspikes);
	markmat1 = zeros(nspikes,reps);
	markmat2 = markmat1;
	spmarkmat1 = markmat1;
	spmarkmat2 = markmat1;
	% replace 0's in binidx with nan's so we can distinguish
	% 0's due to synchronous spikes
	binidx2 = binidx;
	binidx2(binidx2==0) = nan;
	% find the spikes that are within NumCentralBins. Add 1 so that we
	% can just use the < operator instead of the <= operator 
	stdiff = (NumCentralBins+1)/2;
	for repi = repvec
		% get spikes from each cell
		mat1(:,1) = binidx2(:,repi);
		mat2(2,:) = binidx2(:,repvec2(repi))';
		% compute xcorr and take the absolute value so we can more easily
		% search for differences that are within NumCentralBins
		xc = abs(mat1 * mat2);
		% find the spikes that are smaller than stdiff
		[m1,m2] = find(xc<stdiff);
		markmat1(m1,repi) = 1;
		markmat2(m2,repi) = 1;
		% get shift predictor spikes
		mat2(2,:) = binidx2(:,sprepvec2(repi))';
		xc = abs(mat1 * mat2);
		% find the spikes that are smaller than stdiff
		[m1,m2] = find(xc<stdiff);
		spmarkmat1(m1,repi) = 1;
		spmarkmat2(m2,sprepvec(repi)) = 1;
	end
	% convert to sparse matrices to save space
	markmat1 = sparse(markmat1);
	spmarkmat1 = sparse(spmarkmat1);
	markmat2 = sparse(markmat2);
	spmarkmat2 = sparse(spmarkmat2);
	warning on MATLAB:divideByZero

	% save data
	save(resultfile,'data','shiftp','surrmean','surrstd','surrhist', ...
		'dzscore','dxcorr','dxcmean','dxcstd','dxcz','hpval', ...
		'markmat1','spmarkmat1','markmat2','spmarkmat2','psthbins', ...
		'psths','jbmat','limat','dspsthz');
end % if(surrnum<numwindows)
