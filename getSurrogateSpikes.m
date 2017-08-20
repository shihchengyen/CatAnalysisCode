function sptrain = getSurrogateSpikes(rfd,n,varargin)
%getSurrogateSpikes Generate surrogate spikes using q(t)
%   SPTRAIN = getSurrogateSpikes(RFD,CELLN,VARARGIN) generates
%   surrogate spikes for CELLN in the RFD structure. The following
%   fields in RFD are required:
%      rfd.qt(n) - the free firing rate, q(t) for cell n.
%      rfd.repetitions(n) - the number of repetitions for cell n.
%      rfd.wt{n} - the recovery function for cell n.
%      rfd.duration - the duration of q(t).
%      rfd.qtbinsize - the bin size used to compute q(t) in ms.
%      rfd.rtEdges - the time vector for the q(t) function.
%
%   These are the optional input arguments:
%      dbflag - prints out debug information while running.
%      cstep - number of data points to step through at a time
%              when creating spikes (default is the length of 
%              rfd.wt{n}).
%      displayflag - plots the intermediate calculations that were
%                    used to create surrogate spikes.
%
%   sptrain = getSurrogateSpikes(rfd,celln,'dbflag','cstep', ...
%                length(rfd.wt{n}),'displayflag')

wt = rfd.wt{n};
% use length of wt to step through qt
lwt = length(wt);
Args = struct('dbflag',0,'cstep',lwt,'displayflag',0);

Args = getOptArgs(varargin,Args,'flags',{'dbflag','displayflag'});
cstep = Args.cstep;

% initialize rand seed
rand('state',sum(100*clock)); 

% get free firing rate, q(t)
rqt = rfd.qt{n} * rfd.qtbinsize(n) / 1000;
% get number of repetitions
reps = rfd.repetitions(n);
% get recovery function, w(t) and pad with 1's
wt = [wt; ones(cstep-lwt,1)];

% get number of bins in qt
nbins = length(rqt);
% wrap q(t) function around but since we are only moving in windows of
% cstep size, we just need to pad with cstep values
qt = [rqt; rqt(1:cstep)];
% calculate running sum for cstep points
matSum = tril(ones(cstep,cstep));

if(Args.dbflag)
	fprintf('Reps: %i nbins: %i cstep: %i\n',reps,nbins,cstep);
end

% get first spike by using w(t) = 1
% create sptrain so we don't have to keep changing memory size
% assume you can't have more than 1 spike in 1 ms so make sptrain
% equal to the duration of a repetition in ms
rep = 1;
sptrain{rep} = zeros(rfd.duration,1);
spti = 0;
% initialize cend for loop 
cend = 0;
eValue = 0;
cSi = [];
bstop = 0;

if(Args.displayflag)
	% clear the figure
	cla
	% set hold to on
	hold on
end

% get random number uniformly distributed in the range {0,1}
r = rand;
rln = -log(r);
while(1)
	while( isempty(cSi) )
		% restart the indicies if they exceed nbins
		if(cend>nbins)
			% finalize old repetition by removing extraneous zeros
			sptrain{rep} = sptrain{rep}(1:spti);
			fprintf('Finished rep %i\n',rep);
			% start a new repetition
			rep = rep + 1;
			if(rep>reps)
				% set flag so outer while loop will break as well
				bstop = 1;
				break;
			end
			sptrain{rep} = zeros(rfd.duration,1);
			spti = 0;
			cend = cend - nbins;
		end
		% get the next start and end values
		cstart = cend + 1;
		cend = cend + cstep;
		% compute the next cstep cummulative sums
		cSum = eValue + matSum * qt(cstart:cend);
		if(Args.dbflag)
			fprintf('cstart: %i cend: %i eValue: %i rln: %f\n',cstart,cend,eValue,rln);
			% take the transpose since fprintf grabs values by columns
			fprintf('qt: %f cSum: %f\n',[qt(cstart:cend)';cSum']);
		end
		% find if there is an index greater than rln
		cSi = find(cSum>rln);
		if(Args.displayflag & isempty(cSi))
			pxvals = rfd.rtEdges{n}(cstart:cend);
			% plot qt from cstart to cend
			plot(pxvals,qt(cstart:cend),'.-')
			% plot cSum from cstart to cend
			plot(pxvals,cSum,'r.-')
			% plot rln from cstart to cend
			plot(pxvals,rln,'g.-')
		end
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
	if(Args.dbflag)
		fprintf('xrlni: %i\n',xrlni);
	end
	if(Args.displayflag)
		pxvals = rfd.rtEdges{n}(cstart:xrlni);
		plot(pxvals,qt(cstart:xrlni),'.-')
		plot(pxvals,cSum(1:cSi(1)),'r.-')
		plot(pxvals,rln,'g.-')
	end
	% check to see if we need to go to the next repetition
	if(xrlni>nbins)
		% finalize old repetition by removing extraneous zeros
		sptrain{rep} = sptrain{rep}(1:spti);
		fprintf('Finished rep %i\n',rep);
		% create new repetition
		rep = rep + 1;
		if(rep>reps)
			% break out of while loop
			break;
		end
		sptrain{rep} = zeros(rfd.duration,1);
		spti = 0;
		% adjust xrlni to the new repetition's time scale
		xrlni = xrlni - nbins;
	end
	% convert index to spike time
	spt = rfd.rtEdges{n}(xrlni);
	% check if spt is the same as the last spike time
	if( (spti>0) && (spt==sptrain{rep}(spti)) )
		fprintf('Same spike time!\n');
	end
	% increment spike count
	spti = spti + 1;
	% add value to spike train
	sptrain{rep}(spti) = spt;
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
	if(Args.dbflag)
		fprintf('xrlni: %i spti: %i spt: %f\n',xrlni,spti,spt);
		fprintf('cstart: %i cend: %i rln: %f\n',cstart,cend,rln);
		fprintf('qt: %f wt: %f csum: %f\n',[qt(cstart:cend)'; wt'; cSum']);
	end
	% find if there is an index greater than rln
	cSi = find(cSum>rln);
	if(Args.displayflag & isempty(cSi))
		pxvals = rfd.rtEdges{n}(cstart:cend);
		% plot qt from cstart to cend
		plot(pxvals,qt(cstart:cend),'.-')
		% plot cSum from cstart to cend
		plot(pxvals,cSum,'r.-')
		% plot rln from cstart to cend
		plot(pxvals,rln,'g.-')
	end
	% get the last value from the cummulative sum for next
	% calculation
	eValue = cSum(cstep);
end	% end of while(1) loop
