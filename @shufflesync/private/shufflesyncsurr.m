function shufflesyncsurr(sdatafile,varargin)

% define some constants
% the values for the data are stored in the first column of the xcorr 
% matrix as well as the spike times matrix
DATACOL = 1;
% data for the shift predictor is stored in the second column of the xcorr 
% matrix
SHIFTCOL = 2;
% column numbers for startlimits and endlimits
EXT1COL = 1;
EXT2COL = 2;
CEN1COL = 3;
CEN2COL = 4;

% load surrogate data file
load(sdatafile);
% load data file
load(cellpairdatafile);

% needed to compile standalone executable
if(nargin>1)
    startwindow = varargin{1};
    if(ischar(startwindow))
        startwindow = str2num(startwindow);
    end
    endwindow = varargin{2};
    if(ischar(endwindow))
        endwindow = str2num(endwindow);
    end
else
    startwindow = 1;
    endwindow = numwindows;
end

for windown = startwindow:endwindow
    % get bins corresponding to window n
    % use startbins2 and endbins2 so we can eliminate boundary effects
    fmin = startbins2(windown) - 1;
    fmax = endbins2(windown) + 1;
    % find indices corresponding to window n
    % subtract and add 1 to fmin and fmax so we don't have to use <= and >=
    [trialspiken,repspike] = find(binidx>fmin & binidx<fmax);
    % check for empty trialspiken so we don't try to use sub2ind which will
    % fail with empty matrices
    emptytsn = isempty(trialspiken);
    if(~emptytsn)
        % grab the actual binidx values. Can't use the values returned by find 
        % since it will be all 1's since the matrix passed to find is binary		
        trialspikebin = binidx(sub2ind(binidxsize,trialspiken,repspike));
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
        [eurepspike,eurepspikea,eurepspikeb] = unique(repspike);    
        % find the reps belonging to the second cell
        ecell2i = eurepspike>reps;
        ecell1i = ~ecell2i;
        erepnums1 = eurepspike(ecell1i);
        erepnums2 = eurepspike(ecell2i) - reps;
        
        % find spikes in central window
		fmin2 = startbins(windown) - 1;
		fmax2 = endbins(windown) + 1;
		ctsbi = (trialspikebin>fmin2 & trialspikebin<fmax2);
        % grab the actual binidx values. Can't use the values returned by find 
        % since it will be all 1's since the matrix passed to find is binary		
        ctrialspikebin = trialspikebin(ctsbi);
        crepspike = repspike(ctsbi);
		% use the above indices to get bin numbers
		% can't use the values returned by find since the array passed to
		% find is a logical array
		[urepspike,urepspikea,urepspikeb] = unique(crepspike);
		cell2i = urepspike>reps;
		cell1i = ~cell2i;
		repnums1 = urepspike(cell1i);
		repnums2 = urepspike(cell2i) - reps;
	
        % find the reps that have spikes in the central window of at least 1 cell
        % can't use union since it will include reps with spikes in the central
        % window for 1 cell and no spikes in the extended window for the other
        [intrepnums1,sria1,srib2] = intersect(erepnums1,repnums2);
        [intrepnums2,srib1,sria2] = intersect(erepnums2,repnums1);
        repsboth = unique([intrepnums1; intrepnums2]);
        % save rep numbers so we can compute synchrony only on these reps
        repnvec = vecr(repsboth);
        repnl = length(repnvec);
    end
    if(emptytsn || (repnl==0))
        % no pairs of spikes found in any rep so just save matrix of all zeros
        % don't use matrix of nan's since a matrix of zeros is more accurate
        xc = xctemp;
        sptrains = {};
        repswspike1 = [];
        repswspike2 = [];
    else % if(~isempty(trialspiken))
        xcreps = cell(repnl,1);
        % get number of repetitions with spikes in 1st cell
        lrepn1 = length(repnums1);
        % get number of repetitions with spikes in 2nd cell
        lrepn2 = length(repnums2);
        % allocate memory for the required number of repetitions
        sptrains = spiketrains;
        sptrains{1} = cell(lrepn1,1);
        sptrains{2} = cell(lrepn2,1);
        % get number of repetitions with spikes in the extended window of 
        % the 1st cell
        elrepn1 = length(erepnums1);
        % get number of repetitions with spikes in 2nd cell
        elrepn2 = length(erepnums2);
        % allocate memory for the required number of repetitions
        esptrains = spiketrains;
        esptrains{1} = cell(elrepn1,1);
        esptrains{2} = cell(elrepn2,1);
        
		% grab the limits for each group from the output of unique
		ecell1limits = eurepspikea(ecell1i);
		ecell2limits = eurepspikea(ecell2i);
		% add 1 to cell1limits and cell2limits so we don't have to add 1
		% inside the for-loop
		ec1l1 = [0; ecell1limits(1:(elrepn1-1))] + 1;
		ec2l1 = [ecell1limits(elrepn1); ecell2limits(1:(elrepn2-1))] + 1;
		% grab the limits for each group from the output of unique
		cell1limits = urepspikea(cell1i);
		cell2limits = urepspikea(cell2i);
		% add 1 to cell1limits and cell2limits so we don't have to add 1
		% inside the for-loop
        blr1 = lrepn1 > 0;
        if(blr1 && lrepn2>0)
			c1l1 = [0; cell1limits(1:(lrepn1-1))] + 1;
    		c2l1 = [cell1limits(lrepn1); cell2limits(1:(lrepn2-1))] + 1;
        elseif(blr1)
			c1l1 = [0; cell1limits(1:(lrepn1-1))] + 1;
    		c2l1 = [];
        else
			c1l1 = [];
    		c2l1 = [0; cell2limits(1:(lrepn2-1))] + 1;
        end
		
		% set the limits for the indices in trialspikebin
		% this array will keep track of which windows are empty so we don't
		% have to add if statements inside the loop to make sure things work
        startlimits = tmpstartlimits;
        endlimits = tmpendlimits;
		% set start and end limits for extended window in cell 1
		startlimits(erepnums1,EXT1COL) = ec1l1;
		endlimits(erepnums1,EXT1COL) = ecell1limits;
		% set start and end limits for extended window in cell 2
		startlimits(erepnums2,EXT2COL) = ec2l1;
		endlimits(erepnums2,EXT2COL) = ecell2limits;
		% set start and end limits for central window in cell 1
		startlimits(repnums1,CEN1COL) = c1l1;
		endlimits(repnums1,CEN1COL) = cell1limits;
		% set start and end limits for central window in cell 2
		startlimits(repnums2,CEN2COL) = c2l1;
		endlimits(repnums2,CEN2COL) = cell2limits;

		% initialize random number generator which is capable of generating 
		% 2^1492 values so a typical run uses 12x200x3x1000=7200000 random
		% numbers so we won't have to reset the state between calls to randn
		rand('state',sum(windown*clock));
		
        for repi = 1:repnl
            % get rep number
            repn = repsboth(repi);
            % get the spike times for the central window
            c1lindices = startlimits(repn,CEN1COL):endlimits(repn,CEN1COL);
            c2lindices = startlimits(repn,CEN2COL):endlimits(repn,CEN2COL);
            sptimes{2} = ctrialspikebin(c2lindices);
            sptimes{1} = ctrialspikebin(c1lindices);
            nsc1 = length(c1lindices);
            nsc2 = length(c2lindices);
            nspikescells = [nsc1 nsc2];
            
            % get the spike times for the extended window
            ec1lindices = startlimits(repn,EXT1COL):endlimits(repn,EXT1COL);
            ec2lindices = startlimits(repn,EXT2COL):endlimits(repn,EXT2COL);
            espt1 = trialspikebin(ec1lindices);
            espt2 = trialspikebin(ec2lindices);
			esptimes{2} = espt2;
			esptimes{1} = espt1;
			ensc1 = length(espt1);
			ensc2 = length(espt2);
			enspikescells = [ensc1 ensc2];

            % get the spike times for the extended window after removing the
            % spikes in the central window. This makes computing the xcorr 
            % easier since we don't have to subtract out the double-counting 
            % of the central spikes
            excspt1 = setdiff(esptimes{1},sptimes{1});
            excspt2 = setdiff(esptimes{2},sptimes{2});
            excsptimes{2} = excspt2;
            excsptimes{1} = excspt1;
            excnsc1 = length(excspt1);
            excnsc2 = length(excspt2);
            excnspikescells = [excnsc1 excnsc2];
            
            % check spike times
%             concat([full([trialspikebin(ec1lindices) trialspiken(ec1lindices) repspike(ec1lindices)])], ...
%             		[full([trialspikebin(ec2lindices) trialspiken(ec2lindices) repspike(ec2lindices)])], ...
%             		[full([ctrialspikebin(c1lindices) crepspike(c1lindices)])], ...
%             		[full([ctrialspikebin(c2lindices) crepspike(c2lindices)])],'Columnwise')
            
            % create matrix to compute cross-correlation between central 
            % windows of cells 1 and 2
            mat1 = [sptimes{1}(:) repmat(-1,nsc1,1)];
            mat2 = [ones(1,nsc2); sptimes{2}(:)'];
            % compute xcorr for the data
            xcvals = mat1 * mat2;
            nxcvals = nsc1 * nsc2;
                        
            % compute xcorr between central window of cell 1 and extended
            % window of cell 2
            emat2 = [ones(1,excnsc2); excsptimes{2}(:)'];
            excvals1 = mat1 * emat2;
            nexcvals1 = nsc1 * excnsc2;
            
            % compute xcorr between central window of cell 2 and extended
            % window of cell 1
            emat1 = [excsptimes{1}(:) repmat(-1,excnsc1,1)];
            excvals2 = emat1 * mat2;
            nexcvals2 = excnsc1 * nsc2;

            % initialize with nan since 0's mean 0 time difference or
            % synchronous spikes
            xc1 = repmat(nan,nxcvals+nexcvals1+nexcvals2,numSurrogates2);
            xc1(:,DATACOL) = [xcvals(:); excvals1(:); excvals2(:)];
            % initialize cellvecr to cellvec
            cellvecr = cellvec;
            % create array for storing spike times
            if(nsc1==0)
	            % sptrains{1}{repi} = [];
	            cellvecr = 2;
                spiketrains{1} = [];
                noskip1 = 0;
            else
	            sptrains{1}{repi} = zeros(nsc1,numSurrogates1);
				% save the spike times of the data for the psth analysis
				sptrains{1}{repi}(:,DATACOL) = sptimes{1}(:);
                noskip1 = 1;
	        end
	        if(nsc2==0)
            	% sptrains{2}{repi} = [];
            	cellvecr = 1;
                spiketrains{2} = [];
                noskip2 = 0;
	        else
            	sptrains{2}{repi} = zeros(nsc2,numSurrogates1);
				% save the spike times of the data for the psth analysis
				sptrains{2}{repi}(:,DATACOL) = sptimes{2}(:);
                noskip2 = 1;
            end
            % set flags to see if we can skip the while loop
            % do this outside the for surri = surrvec loop so we don't do
            % it over and over again
            nfc = needrefcheck;
            for celli = cellvecr
                % only do while loop if there is more than 1 spike
                % and the intervals between the spikes is such that 
                % there will be overlaps
                bspikes = enspikescells(celli) > 1;
                bisi = ~isempty(find(diff(esptimes{celli})<shufflevecl));
                % check if there are isi's smaller than 1 in the extended
                % window
                nref = length(find(diff(excsptimes{celli})<1));
                if( bspikes && bisi)
					nfc(celli) = 1;
                end
            end

            % loop over number of surrogates
            for surri = surrvec
            	surri1 = surrvec1(surri);
                surri2 = surrvec2(surri);
                for celli = cellvecr
                    % check if we need to check for refractory period 
                    % violations
                    if(nfc(celli))
                        % make sure we go into the while loop
                        nisi0 = 1;
                        % check for refractory period violations
                        while(nisi0)
                            % get random numbers for cell i
                            randnums = floor(rand(nspikescells(celli),1) ...
                                * shufflevecl) - shuffle;
                            % add random numbers to spike times for cell i and
                            % sort the spike times
                            spt1 = sort(sptimes{celli}+randnums);
                            % take the diff and find 0's and then subtract
                            % the pre-existing 0's to see if they were 
                            % introduced by the surrogates
                            % no need to worry about negative numbers since
                            % the first number is always going to be larger
                            % or equal to nref
                            nisi0 = length(find(diff(spt1)==0)) - nref;
                        end % while(nisi0)
                        spiketrains{celli} = spt1;
                    else % if(nfc(celli))
                        % don't have to worry about overlap between the 
                        % shuffled spikes so just add the shuffle to the
                        % spike times
                        spiketrains{celli} = sort(sptimes{celli} + ...
                            floor(rand(nspikescells(celli),1) ...
                            * shufflevecl) - shuffle);
                    end % if(nfc(celli))
                end % for celli = cellvec
                % now set the first cell's spiketrains
                if(noskip1)
                    mat1(:,1) = spiketrains{1}(:);
                end
                % now set the second cell's spiketrains
                if(noskip2)
                    mat2(2,:) = spiketrains{2}(:)';
                end
                % now compute the cross-correlation between the first and
                % second cell
                xcvals = mat1 * mat2;                
                % compute xcorr between central window of cell 1 and extended
                % window of cell 2
                excvals1 = mat1 * emat2;
                % compute xcorr between central window of cell 2 and extended
                % window of cell 1
                excvals2 = emat1 * mat2;
                
                % save the differences so we can do a histogram at the end
                % instead of after every surrogate
                xc1(:,surri2) = [xcvals(:); excvals1(:); excvals2(:)];
                % store the spike times for psth analysis
                if(noskip1)
	                sptrains{1}{repi}(:,surri1) = spiketrains{1}(:);
	            end
	            if(noskip2)
	                sptrains{2}{repi}(:,surri1) = spiketrains{2}(:);
	            end
            end % loop over number of surrogates
            xcreps{repi} = xc1;
        end % for repi = 1:repnl
        % convert cell array xcreps to matrix
        xcmat = cell2mat(xcreps);
        % check if xcmat is a single row (i.e. only 1 rep had exactly 1
        % spike in each cell)
        if(size(xcmat,1)==1)
            % add a second row of nan so the histogram will be taken in
            % columns. Otherwise, hist will collapse the values across the
            % entire set of surrogates and return xc which is 1 x lagslength
            % instead of lagslength x (numsurrogates+1)
            xcmat = concatenate(xcmat,nan);
        end
        % do histogram of the values in xcreps
        xc = histcie(xcmat,lagbins2,'DropLast');        
        
        % now compute the shift predictor
        % subtract rep number for the second cell so that the 1st rep of
        % cell 1 is compared to 2nd rep of cell 2
        sperepnums2 = mod(erepnums2 - 2,reps) + 1;
        % subtract the rep number for the central window for the second cell
        sprepnums2 = mod(repnums2 - 2,reps) + 1;
		% find the reps that have spikes in the central window of at least 1 cell
		% can't use union since it will include reps with spikes in the central
		% window for 1 cell and no spikes in the extended window for the other
		[spintrepnums1,spsria1,spsrib2] = intersect(erepnums1,sprepnums2);
		[spintrepnums2,spsrib1,spsria2] = intersect(sperepnums2,repnums1);
		sprepsboth = unique([spintrepnums1; spintrepnums2]);
		% save rep numbers so we can compute synchrony only on these reps
		sprepnvec = vecr(sprepsboth);
		sprepnl = length(sprepnvec);
        
        if(sprepnl>0)
            xcreps = cell(sprepnl,1);
			% set the limits for the indices in trialspikebin
			startlimits(:,[EXT2COL CEN2COL]) = circshift(startlimits(:,[EXT2COL CEN2COL]),-1);
			endlimits(:,[EXT2COL CEN2COL]) = circshift(endlimits(:,[EXT2COL CEN2COL]),-1);

            for repi = 1:sprepnl
				% get rep number
				repn = sprepsboth(repi);
				% get spikes times for the central window
				c1lindices = startlimits(repn,CEN1COL):endlimits(repn,CEN1COL);
				c2lindices = startlimits(repn,CEN2COL):endlimits(repn,CEN2COL);
				sptimes1 = ctrialspikebin(c1lindices);
				sptimes2 = ctrialspikebin(c2lindices);
				nsc1 = length(c1lindices);
				nsc2 = length(c2lindices);
				
				% get the spike times for the extended window
				ec1lindices = startlimits(repn,EXT1COL):endlimits(repn,EXT1COL);
				ec2lindices = startlimits(repn,EXT2COL):endlimits(repn,EXT2COL);
				espt1 = trialspikebin(ec1lindices);
				espt2 = trialspikebin(ec2lindices);
	
				% get the spike times for the extended window after removing the
				% spikes in the central window. This makes computing the xcorr 
				% easier since we don't have to subtract out the double-counting 
				% of the central spikes
				excsptimes1 = setdiff(espt1,sptimes1);
				excsptimes2 = setdiff(espt2,sptimes2);
				excnsc1 = length(excsptimes1);
				excnsc2 = length(excsptimes2);
				
				% check spike times
% 				concat([full([trialspikebin(ec1lindices) trialspiken(ec1lindices) repspike(ec1lindices)])], ...
% 						[full([trialspikebin(ec2lindices) trialspiken(ec2lindices) repspike(ec2lindices)])], ...
% 						[full([ctrialspikebin(c1lindices) crepspike(c1lindices)])], ...
% 						[full([ctrialspikebin(c2lindices) crepspike(c2lindices)])],'Columnwise')
				
				% create matrix to compute cross-correlation between central 
				% windows of cells 1 and 2
				mat1 = [sptimes1(:) repmat(-1,nsc1,1)];
				mat2 = [ones(1,nsc2); sptimes2(:)'];
				% compute xcorr for the data
				xcvals = mat1 * mat2;
							
				% compute xcorr between central window of cell 1 and extended
				% window of cell 2
				emat2 = [ones(1,excnsc2); excsptimes2(:)'];
				excvals1 = mat1 * emat2;
				
				% compute xcorr between central window of cell 2 and extended
				% window of cell 1
				emat1 = [excsptimes1(:) repmat(-1,excnsc1,1)];
				excvals2 = emat1 * mat2;
				xcreps{repi} = [xcvals(:); excvals1(:); excvals2(:)];
            end % for repi = 1:repnl
            % do histogram of the values in xcreps
            xc(:,SHIFTCOL) = vecc(histcie(full(cell2mat(xcreps)),lagbins2,'DropLast'));
        end % if(~isempty(trialspiken) && (repnl>0))
        
        % now shuffle the spike times for the remaining reps that are
        % not in repnvec
		% create cell array to store cell1limits and cell2limits and its
		% corresponding c1l1 and c2l1
		clcell = {cell1limits,cell2limits};
		cl1cell = {c1l1,c2l1};
		eclcell = {ecell1limits,ecell2limits};
		ecl1cell = {ec1l1,ec2l1};
		repnumscell = {repnums1,repnums2};
		erepnumscell = {erepnums1,erepnums2};
        for celli = cellvec
			% determine the repetitions that are not covered by repnvec,
			% i.e. reps that do not have spikes in both cells
			[repsnopair{celli},srinopair] = setdiff(repnumscell{celli},repsboth);
			% sri1nopair = srinopair + 1;
			% find the indices for the extended windows
			[discard,esrinopair] = intersect(erepnumscell{celli},repsnopair{celli});
			% esri1nopair = esrinopair + 1;
			a = repsnopair{celli};
			for repi = 1:length(repsnopair{celli})
				repi2 = repnl+repi;
				% get the spike times for the central window
                srepi = srinopair(repi);
				clindices = cl1cell{celli}(srepi) : clcell{celli}(srepi);
				sptime = ctrialspikebin(clindices);
                sptrains{celli}{repi2}(:,DATACOL) = sptime;
				nspikes = length(sptime);
				% get the spike times for the extended window
                esrepi = esrinopair(repi);
				eclindices = ecl1cell{celli}(esrepi) : eclcell{celli}(esrepi);
				esptime = trialspikebin(eclindices);
				enspikes = length(esptime);
				% get the spike times for the extended window without the
				% spikes from the central window
				excsptime = setdiff(esptime,sptime);
				
                % check the spike times
                % concat(full([sptime crepspike(clindices)]),full([esptime repspike(eclindices)]),'Columnwise')
                
                % only do while loop if there is more than 1 spike
                % and the intervals between the spikes is such that 
                % there will be overlaps
                bspikes = enspikes > 1;
                bisi = ~isempty(find(diff(esptime)<shufflevecl));
                % check if there are isi's smaller than 1 in the extended
                % window
                nref = length(find(diff(excsptime)<1));
                if( bspikes && bisi)
					nfc = 1;
				else
					nfc = 0;
                end

				% loop over number of surrogates
				for surri = surrvec
					surri1 = surrvec1(surri);
					% check if we need to check for refractory period 
					% violations
					if(nfc)
						% make sure we go into the while loop
						nisi0 = 1;
						% check for refractory period violations
						while(nisi0)
							% get random numbers for cell i
							randnums = floor(rand(nspikes,1) ...
								* shufflevecl) - shuffle;
							% add random numbers to spike times for cell i and
							% sort the spike times
							spt1 = sort(sptime+randnums);
							% take the diff and find 0's
                            % no need to worry about negative numbers since
                            % the first number is always going to be larger
                            % or equal to nref
							nisi0 = length(find(diff(spt1)==0)) - nref;
						end % while(nisi0)
						sptrains{celli}{repi2}(:,surri1) = spt1;
					else % if(nfc)
						% don't have to worry about overlap between the 
						% shuffled spikes so just add the shuffle to the
						% spike times
						sptrains{celli}{repi2}(:,surri1) = sort(sptime + ...
							floor(rand(nspikes,1) ...
							* shufflevecl) - shuffle);
					end % if(nfc)
				end % for surri = surrvec
			end % for repi = 1:length(repsnopair{celli})
		end % for celli = [1 2]

        repswspike1 = [repsboth; repsnopair{1}(:)];
		repswspike2 = [repsboth; repsnopair{2}(:)];
    end % if(~isempty(trialspiken) && ~isempty(repnvec))
    % convert xc to sparse matrix to save space
    xc = sparse(xc);
    % save data to file
    save([surrprefix num2str(windown,'%04d')],'xc','sptrains', ...
    	'repswspike1','repswspike2');
end % for windown = startwindow:endwindow
