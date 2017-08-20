function [qt,Wt,rt,rtEdges,reps] = getFreeFiringRate(rasters,rec,varargin)
%getFreeFiringRate Get free firing rate for a set of repetitions
%
%   [QT,WT,RT,RTEDGES,REPS] = getFreeFiringRate(RASTERS,REC,VARARGIN) 
%   computes the free firing rate using RASTERS, a cell array of spike 
%   times across repetitions, and REC, a vector containing the values
%   of the recovery function, w(t). These are the output arguments:
%      QT - the free firing rate, which is RT divided by WT.
%      WT - the spike probability function.
%      RT - the PSTH. The last bin returned by histcie is dropped so
%           that RT is the same length as QT and WT.
%      RTEDGES - the edges used to compute the PSTH (using histcie).
%      REPS - the number of repetitions in RASTERS.
%
%   These are the optional input arguments:
%      qtbinsize - size of time bins used to compute the various 
%                 functions (default is 0.2 ms).
%      duration - the duration of a repetition (default is 30000 ms).
%
%   [QT,WT,RT,RTEDGES,REPS] = getFreeFiringRate(RASTERS,REC,VARARGIN)

% default values for optional arguments
Args = struct('qtbinsize',0.2,'duration',30000);

% get optional arguments
Args = getOptArgs(varargin,Args);

% get number of repetitions
reps = length(rasters);

% get length of recovery function
nrec = length(rec);
% get vector of 1:nrec so we don't have to keep generating it inside 
% the loop
recvec = 1:nrec;

% compute length of Wt
lWt = round(Args.duration/Args.qtbinsize);
% create Wt
Wt = zeros(lWt,1);
% create initial Wjt
Wt1 = ones(lWt,1);
% create row vector of 1's
h1 = Wt1(1:nrec)';

% get bins for PSTH
rtEdges = (0:Args.qtbinsize:Args.duration)';
% number of bins is going to be 1 more than lWt since we are going from
% 0 to Args.duration and histcie retains the last element.
lrt =  lWt + 1;
% create rt
rt = zeros(lrt,1);

% to get the spike probability function Wjt first convert spike times
% to indices into the Wt vector. Assume the binsize is 2 ms, and we have
% a spike at time 1 ms, then when we do floor(spt/Args.qtbinsize) we will
% get a index of 0, when it should be index 1. However, when we create 
% matrix 1 which is [spi1 1; spi2 1; ...; spin 1] and matrix2 which is
% [1 1 .... 1; 1 2 ... nrec] and take the matrix multiplication, we get
% [spi1+1 spi1+2 ... spi1+nrec; spi2+1 spi2+2 ... spi2+nrec; ...] which
% correctly adds 1 to the indices. We take the transpose before reshaping
% the matrix into a column vector since reshape grabs column by column.

% loop over reps
for r = 1:reps
	% get spike times for this repetition
	spt = vecc(rasters{r});
	% get number of spikes
	nsp = length(spt);
	% skip if there are no spikes
	if (nsp>0)
		% convert spike times to indices at the specified time resolution
		spi = floor(spt/Args.qtbinsize);
		% create matrix1
		m1 = [spi Wt1(1:nsp)];
		
		% create matrix2
		m2 = [h1; recvec];
		% do matrix multiplication
		m = m1 * m2;
		% reshape into appropriate column vector
		mi = reshape(m',[],1);
		% replicate recovery function nsp times
		mv = repmat(rec,nsp,1);	
		
		% initialize Wjt with all ones
		Wjt = Wt1;
		% insert the spike probability function
		% make sure indices in mi don't exceed lWt
		miX = find(mi>lWt);
		if isempty(miX)
			Wjt(mi) = mv;
		else
			miend = miX - 1;
			Wjt(mi(1:miend)) = mv(1:miend);
		end
	
		% add Wjt to Wt
		Wt = Wt + Wjt;
		
		% compute PSTH
		rjt = histcie(spt,rtEdges);
		% make sure rjt is a column vector if spt contains only 1 spike
		% get size of rjt
		rs = size(rjt);
		if (rs(1)<rs(2))
			rjt = rjt';
		end
		% add rjt to rt
		rt = rt + rjt;
	end
end

% divide by reps even if there were repetitions without spikes since 
% the absence of spikes also tells us something
Wt = Wt/reps;
% divide by qtbinsize to yield rate instead of count
% drop the last bin value which should be 0 because we used histcie
rt = rt(1:lWt)/(Args.qtbinsize*reps);
qt = rt ./ Wt;
