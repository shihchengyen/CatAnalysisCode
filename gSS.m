function gSS(n,sets,filename,cstep,qscale)
%   gSS(n,sets,filename,cstep,qscale)

% make sure arguments are numeric
if(ischar(n))
	n = str2num(n);
end
if(ischar(sets))
	sets = str2num(sets);
end
if(ischar(cstep))
	cstep = str2num(cstep);
end
if(ischar(qscale))
	qscale = str2num(qscale);
end

% print out arguments
fprintf('Cell# %i filename %s cstep %i qscale %f\n',n,filename,cstep,qscale);
load('rfd');
wt = rfd.wt{n};
% use length of wt to step through qt
lwt = length(wt);

% initialize rand seed
rand('state',sum(100*clock)); 

% get free firing rate, q(t)
rqt = rfd.qt(:,n) * qscale;
% get number of repetitions
reps = rfd.repetitions(n);
% get recovery function, w(t) and pad with 1's
wt = [wt; ones(cstep-lwt,1)];

% get number of bins in qt
nbins = rfd.duration/rfd.qtbinsize;
% wrap q(t) function around but since we are only moving in windows of
% cstep size, we just need to pad with cstep values
qt = [rqt; rqt(1:cstep)];
% calculate running sum for cstep points
matSum = tril(ones(cstep,cstep));

% pre-allocate memory
sptrain = cell(1,sets);

% loop on sets
for setn = 1:sets	
	% create sptrain so we don't have to keep changing memory size when
	% looping within one repetition
	sptrain{1,setn} = cell(reps,1);
	% assume you can't have more than 1 spike in 1 ms so make sptrain
	% equal to the duration of a repetition in ms
	rep = 1;
	sptrain{1,setn}{rep} = zeros(rfd.duration,1);
	spti = 0;
	% initialize cend for loop 
	cend = 0;
	eValue = 0;
	cSi = [];
	bstop = 0;
	
	% get random number uniformly distributed in the range {0,1}
	% get first spike by using w(t) = 1
	r = rand;
	rln = -log(r);
	while(1)
		while( isempty(cSi) )
			% restart the indicies if they exceed nbins
			if(cend>nbins)
				% finalize old repetition by removing extraneous zeros
				sptrain{1,setn}{rep} = sptrain{1,setn}{rep}(1:spti);
				fprintf('Finished set %i rep %i\n',setn,rep);
				% start a new repetition
				rep = rep + 1;
				if(rep>reps)
					% set flag so outer while loop will break as well
					bstop = 1;
					break;
				end
				sptrain{1,setn}{rep} = zeros(rfd.duration,1);
				spti = 0;
				cend = cend - nbins;
			end
			% get the next start and end values
			cstart = cend + 1;
			cend = cend + cstep;
			% compute the next cstep cummulative sums
			cSum = eValue + matSum * qt(cstart:cend);
			% find if there is an index greater than rln
			cSi = find(cSum>rln);
			% get the last value from the cummulative sum for next
			% calculation
			eValue = cSum(cstep);
		end
		if(bstop)
			% inner loop exceeded reps so stop
			break;
		end
		% value found so figure out where to put the spike
		% get index that first exceeds rln
		% if it was index 22, cstart will be 21, cSi(1) will be 2 so in 
		% order to get back 22, we subtract 1 from cstart + cSi(1)
		xrlni = cstart+cSi(1)-1;
		% check to see if we need to go to the next repetition
		if(xrlni>nbins)
			% finalize old repetition by removing extraneous zeros
			sptrain{1,setn}{rep} = sptrain{1,setn}{rep}(1:spti);
			fprintf('Finished set %i rep %i\n',setn,rep);
			% create new repetition
			rep = rep + 1;
			if(rep>reps)
				% break out of while loop
				break;
			end
			sptrain{1,setn}{rep} = zeros(rfd.duration,1);
			spti = 0;
			% adjust xrlni to the new repetition's time scale
			xrlni = xrlni - nbins;
		end
		% convert index to spike time
		spt = rfd.rtEdges(xrlni);
		% increment spike count
		spti = spti + 1;
		% add value to spike train
		sptrain{1,setn}{rep}(spti) = spt;
		% get new random number
		r = rand;
		rln =  -log(r);	
		% reset cstart, cend, eValue and cSi
		cstart = xrlni;
		% if cstep is 10 and cstart is 22 then cend = 22 + 10 - 1
		% will be 10 values
		cend = cstart + cstep - 1;
		% find cummulative sum again taking into account the
		% relative refractory period
		cSum = matSum * (qt(cstart:cend) .* wt);
		% find if there is an index greater than rln
		cSi = find(cSum>rln);
		% get the last value from the cummulative sum for next
		% calculation
		eValue = cSum(cstep);
	end	% end of while(1) loop
end % end of loop on sets

%save(filename,'sptrain')
fid = fopen([filename '.bin'],'w','ieee-le');
fwrite(fid,[sets reps],'int32');
for i = 1:sets
	for j = 1:reps
		st = sptrain{i}{j};
		% write number of spikes for set i, repetition j
		fwrite(fid,size(st,1),'int32');
		% write spike times for set i, repetition j
		% multiply by 10 so we can convert tenth of a ms to integer
		fwrite(fid,st*10,'int32');
	end
end
fclose(fid);
