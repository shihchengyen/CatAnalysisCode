function asxcSurrSynchrony(startwindow,endwindow,datafile,surrprefix)

% needed to compile standalone executable
if(ischar(startwindow))
	startwindow = str2num(startwindow);
end
if(ischar(endwindow))
	endwindow = str2num(endwindow);
end

% load data file
mat = load(datafile);
% get histogram indices
idx = mat.binidx;
% get reps, e.g. 100
reps = mat.reps;
% get 2*reps which is used pretty often so we don't have to keep calculating it
% e.g. 200
tworeps = mat.tworeps;
% get total number of bins in analysis windows, e.g. 36
wlength = mat.windowlength;
% get number of bins in a subwindow, e.g. 12
swinbins = mat.SubWindowBins;
% get number of subwindows
subwindows = mat.subwindows;
% get start and end bins for each window
sbins = mat.startbins;
ebins = mat.endbins;
% get bins used to compute histogram of the xcorr values, e.g. -35:1:35
lagbins = mat.lagbins;
% get length of lagbins, e.g. 2*35+1
lagslength = mat.lagslength;
% get vector of incrementing subwindow bins, e.g. 1:12
swbvector = mat.swbvector;
% matrix of nan's used to store window bin numbers with spikes, 
% e.g. size is 36x200
winspikes = mat.windowspikes;
wssize = mat.wssize;
% matrix of nan's used to store subwindow bin numbers with spikes, 
% e.g. size is 12x200x3
swinspikes = mat.subwindowspikes;
swssize = mat.swssize;
% matrix of subwindow indices, e.g. size is [1:12]'x200x3
subwindowindices = mat.subwindowindices;
% matrix 1 used to compute xcorr, e.g. size is 36x2
m1 = mat.m1;
% matrix 2 used to compute xcorr, e.g. size is 2x36
m2 = mat.m2;
% matrix to store xcorr for data and surrogates, e.g. size is (35*2+1)x1001
xctemp = mat.xctemp;
% size of xc
dnssize = mat.dnssize;
% vector to store histogram of one particular xcorr
nxctemp = mat.nxctemp;
% vector contain indices from 2 to SurrSets+1 used to index into xc
setnvec = mat.setnvec;
% vector containing repetition indices for second cells, e.g. 101:200
repnvec2 = mat.repnvec2;

for windown = startwindow:endwindow
	% get bins corresponding to window n
	fmin = sbins(windown);
	fmax = ebins(windown);
	% find indices corresponding to window n
	% subtract and add 1 to fmin and fmax so we don't have to use <= and >=
	[trialspiken,repspike] = find(idx>(fmin-1) & idx<(fmax+1));
	% check to see if no indices were found, which means that there were no
	% spikes in either cell in this frame, or if there were no repetitions where
	% there were spikes in both cells in this frame
	% we first do unique to identify repetitions with spikes
	% then we subtract 1 and take the mod using the number of repetitions, e.g. if
	% there were 100 reps, columns 1 to 100 become 0 to 99 and columns 101 to 200 
	% become 0 to 99
	% then we sort the output so that if there were at least one spike in both cells
	% we will have a repeated repetition number, which will show up as a 0 in the
	% diff
	% find the unique cols (i.e. repetition)
	[urepspike,urepspikea,urepspikeb] = unique(repspike);
	spikereps = sort(mod(urepspike-1,reps));
	spikerepsi = find(diff(spikereps)==0);
	% save rep numbers so we can compute synchrony only on these reps
	repnvec = vecr(spikereps(spikerepsi)+1);
	if(isempty(trialspiken) || isempty(repnvec))
		% no pairs of spikes found in any rep so just save matrix of all zeros
		% don't use matrix of nan's since a matrix of zeros is more accurate
		xc = xctemp;
	else % is(~isempty(trialspiken))
		% grab the actual idx values. Can't use the values returned by find 
		% since it will be all 1's since the matrix passed to find is binary		
		trialspikebin = idx(sub2ind(size(idx),trialspiken,repspike));
		% find the first row number for each col. Can't just take the 2nd output 
        % of unique to get the first row for each column since it instead returns 
        % the last row for each column and we can't just add 1 to it since 
        % there might entries where there is only 1 unique value, i.e. last and
        % first row for each unique value is the same
		firstrowperrep = [1; 1+find(diff(repspike))];
		% create corresponding vector which is the first row of each column
		trialspiken1 = trialspiken(firstrowperrep(urepspikeb));
		% subtract the first row from each entry in trialspiken and then add 1 
		% to make sure it is an index
		winspiken = trialspiken - trialspiken1 + 1;
		% we should add 1 to windowspikebin to go from 0-indexed to 1-indexed
		% but since we are taking the subtraction between the values here, and
		% because we are dividing by swinbins and taking the floor to convert 
		% the bin numbers in the window to bin numbers in the subwindow, we are
		% going to leave out the addition of 1 so we don't have to subtract 1
		% again later
		windowspikebin = trialspikebin - fmin;
		winspikes(sub2ind(wssize,winspiken,repspike)) = windowspikebin;
		% compute xcorr for real data
		nxc = nxctemp;
		for repn = repnvec
			% replace data in matrices with spike times from each
			% repetition
			m1(:,1) = winspikes(:,repn);
			m2(2,:) = winspikes(:,repnvec2(repn))';
			mxc = m1 * m2;
			nxc = nxc + hist(mxc(:),lagbins)';
		end
		% first column in xc is for data
		xc(:,1) = nxc;
		% want to create 3-D matrix with repetition number in columns, bins in
		% sub-window in rows and sub-windows in pages. This allows us to do a .*
		% operation with another matrix which is essentially randperm in each
		% column of every page. This will produce surrogate spike trains that
		% observe an absolute refractory period without having to explicitly check
		% for violations.
		% The tricky part in all of this is to create increment indices so that
		% multiple spikes that appear in one sub-window are incremented properly,
		% e.g.   spike #   repetition #   sub-window #   increment #
		%          1            1             1              1
		%          2            1             1              2
		%          3            1             2              1
		%          4            1             2              2
		%          5            1             3              1
		%          6            1             3              2
		%          1            2             1              1
		%          2            2             1              2
		%          3            2             1              3
		%          4            2             2              1
		% get sub-window number (0-indexed). Don't have to subtract 1 from
		% windowspikebin since we didn't add 1 earlier so the division will
		% return the right value
		swinspikebin0 = floor(windowspikebin/swinbins);
		% change to 1-indexed
		swinspikebin = swinspikebin0+1;
		% create unique sub-window number
		repspikewlength = (repspike-1)*wlength;
		repswinspikebin = repspikewlength + swinspikebin;
		% ga -> unique values in sub-window number, repswinspikebin
		% gb -> indices before the sub-window number changes
		% gc -> indices that correspond to values in ga
		[urepswinspikebin,urepswinspikebina,urepswinspikebinb] ...
			= unique(repswinspikebin);
		repswinspikebinl = length(repswinspikebin);
		% produce incremental numbers 
		repswinspikebini = (1:repswinspikebinl)';
		% find the indices corresponding to the second unique value on since the
		% incremental values associated with the first unique value are fine and
		% don't have to be changed.
		repswinspikebini2 = (urepswinspikebina(1)+1):repswinspikebinl;
		% find the incremental values corresponding to the transition to an unique
		% value
		transincrvalue = urepswinspikebina(urepswinspikebinb(repswinspikebini2)-1);
		% subtract the appropriate incremental values from the actual incremental
		% values to get the correct incremental values for spikes within a
		% sub-window
		swinspikeind = [repswinspikebini(1:urepswinspikebina(1)); ...
				repswinspikebini(repswinspikebini2)-transincrvalue];
		% create two matrices of the appropriate size, one to multiply by a
		% randperm matrix to get bin numbers within a sub-window, which will then
		% be added to another matrix to produce bin numbers across sub-windows,
		% i.e. bin 2 in second sub-window will become bin 12+2=14
		swinspikes1 = swinspikes;
		swinspikevals = swinspikes;
		swinspikeind = sub2ind(swssize,swinspikeind,repspike,swinspikebin);
		% change the values in the appropriate indices to 1
		swinspikes1(swinspikeind) = 1;
		% change the values in swinspikevals to correspond to first bin number 
		% in the appropriate sub-window, e.g. 0, 12, 24
		swinspikevals(swinspikeind) = (swinspikebin-1)*swinbins;

		% initialize random number generator which is capable of generating 
		% 2^1492 values so a typical run uses 12x200x3x1000=7200000 random
		% numbers so we won't have to reset the state between calls to randn
		rand('state',sum(windown*clock));
		
		for setn = setnvec
			% create matrix of random numbers
			rmat = randn(swinbins,tworeps,subwindows);
			% use the indices from sorting the random numbers to reorder the randperm
			% matrix
			[csort,csi] = sort(rmat);
			swrand = subwindowindices(csi);
			% multiply the 2 matrices to produce appropriate bin numbers for each
			% sub-window
			swinspikebins = swinspikes1 .* swrand;
			% now add the 2 matrices togther
			surrspikes = swinspikevals + swinspikebins;
			% reshape back to 2-D matrix with bins containing spikes in rows and each
			% column is a repetition
			ss1 = sort(reshape(permute(surrspikes,[1 3 2]),wssize));
			nxc = nxctemp;
			for repn = repnvec
				% replace data in matrices with spike times from each
				% repetition
				m1(:,1) = ss1(:,repn);
				m2(2,:) = ss1(:,repnvec2(repn))';
				mxc = m1 * m2;
				nxc = nxc + hist(mxc(:),lagbins)';
			end
			xc(:,setn) = nxc;
		end % for setn = setnvec
	end % if(isempty(trialspiken) || isempty(repnvec))
	% convert xc to sparse matrix to save space
	xc = sparse(xc);
	% save data to file
	save([surrprefix num2str(windown,'%04d')],'xc');
end % for windown = startwindow:endwindow
