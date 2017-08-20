function [qt,Wt,rt,rtEdges,reps,rtbinsize,rec,nbins] = getFreeFiringRate2(cellstruct,varargin)
%getFreeFiringRate2 Get free firing rate for a neuron
%   [QT,WT,RT,RTEDGES,REPS,RTBINSIZE,REC,QTFRAMEBINS] = getFreeFiringRate2(CELLSTRUCT,
%   VARARGIN) computes the free firing rate using CELLSTRUCT, a 
%   structure containing the following fields:
%      CELLSTRUCT.stimulus_info.adjusted_frames_vector
%      CELLSTRUCT.stimulus_info.start_frame
%      CELLSTRUCT.stimulus_info.end_frame
%      CELLSTRUCT.cell_info.adjusted_spiketrain
%   These are the output arguments:
%      QT - the free firing rate, which is RT divided by WT.
%      WT - the spike probability function.
%      RT - the PSTH. The last bin returned by histcie is dropped so
%           that RT is the same length as QT and WT.
%      RTEDGES - the edges for displaying QT, WT and RT (using histcie).
%      REPS - the number of repetitions in the spike train.
%      RTBINSIZE - actual bin size used to scale RT.
%      REC - the values of the recovery function.
%
%   These are the optional input arguments:
%      (Numeric) - the first numeric argument will be used as the 
%                  recovery function, wt. 
%      qtbinsize - size of time bins used to compute the various 
%                 functions in ms (default is 0.2).
%      duration - the duration of a repetition in ms (default is 30000).
%      curvefitobj - curve fit object used to compute recovery function.
%                    The presence of this argument overides the numeric
%                    argument.
%      FrameSplit - integer specifying how the free firing rate within a
%                   a frame should be split before averaging (e.g.
%                   'FrameSplit',2 would create firing rates averaged by
%                   approximately 35/2 = 17 ms) (default is 0).
%
%   [qt,Wt,rt,rtEdges,reps,rtbinsize] = getFreeFiringRate2(cellstruct,[], ...
%      'qtbinsize',0.2,'duration',30000,'curvefitobj',[],'frameavg');

% default values for optional arguments
Args = struct('ISIBinSize',0.2,'duration',30000,'curvefitobj',[], ...
    'FrameSplit',0);

% get optional arguments
Args = getOptArgs(varargin,Args);

% default is no curve fit object
curvefit = 0;
% if there is a curve fit object use that instead
if(~isempty(Args.curvefitobj))
	curvefit = 1;
end
if(~isempty(Args.NumericArguments))
	% if there is an numeric argument assume the first one is rec
	rec = Args.NumericArguments{1};
else
	rec = [];
end

% get rtEdges by getting difference between adjusted_frame_vectors
% and then creating bins as close to those requested
afv = vecc(cellstruct.stimulus_info.adjusted_frames_vector);
% get length of afv 
afvl = length(afv);
dafv = diff(afv);
% get number of bins if we used specified binsize
% make sure data is in column then take mean to get just 1 number
% hopefully there is no reason the nbins will be different for any of
% the frames
nbins = round(mean(dafv/Args.ISIBinSize));
% get size of bins
binsizes = dafv / nbins;
% create matrix of bin limits
mat1 = tril(ones(nbins));
mat2 = [afv(1:(afvl-1))'; repmat(binsizes',(nbins-1),1)];
blimits = mat1 * mat2;
% append time of end of last frame
histEdges = [reshape(blimits,[],1); afv(afvl)];

% compute histogram
[binSpikeCounts,spikeind] = histcie(cellstruct.cell_info.adjusted_spiketrain,histEdges,'DropLast');
% remove any spikes that might be outside our bins of interest, e.g. any
% spikes in the blank period in the old data
spikeind = spikeind(spikeind~=0);
% compute number of frames presented
nframes = cellstruct.stimulus_info.end_frame - cellstruct.stimulus_info.start_frame + 1;
% compute number of bins in each repetition
repbins = nframes * nbins;
% reshape binSpikeCounts into repetitions and then take transpose so we
% can take the mean of each column and get the mean spike count in each bin
bSC1 = reshape(binSpikeCounts,repbins,[])';
% compute reps since we need to return that value
reps = size(bSC1,1);
% get mean spike counts
bSC1mean = mean(bSC1);
% get rt by dividing by the average bin size, which hopefully should all
% be pretty similar. We are dividing by average bin size to avoid small
% differences in the rate due to unequal bin sizes
rtbinsize = mean(binsizes);
rt = (bSC1mean / rtbinsize * 1000)';

% compute length of binSpikeCounts which should be nframes*nbins*reps
lWt = length(binSpikeCounts);
% create initial Wjt which is all ones and the same length as binSpikeCounts
Wjt = ones(lWt,1);
% get number of spikes
nSpikes = length(spikeind);
% if there is a curve fit object generate recovery function using rtbinsize
if(curvefit)
	fv = getFunctionValues(Args.curvefitobj,'binsize',Args.ISIBinSize);
	if(isempty(rec))
		rec = fv;
	else
        % alter fv to observe the absolute refractory period
        % find first non-zero value
        fnz = find(rec);
        fv(1:(fnz(1)-1)) = 0;
		rec = fv;
	end
end
% make sure recovery function is column vector
rec = vecc(rec);
% get length of recovery function
nrec = length(rec);
% get vector of 1:nrec so we don't have to keep generating it inside 
% the loop
recvec = 0:(nrec-1);
% create row vector of 1's
h1 = ones(1,nrec);
% create matrix1
m1 = [spikeind ones(nSpikes,1)];
% create matrix2
m2 = [h1; recvec];
% do matrix multiplication
m = m1 * m2;
% reshape into appropriate column vector
mi = reshape(m',[],1);
% replicate recovery function nsp times
mv = repmat(rec,nSpikes,1);	
% insert the spike probability function
% make sure indices in mi don't exceed lWt
miX = find(mi>lWt);
if isempty(miX)
	Wjt(mi) = mv;
else
	miend = miX - 1;
	Wjt(mi(1:miend)) = mv(1:miend);
end
% reshape Wjt into repetition and then take transpose so we can take the 
% mean of each column and get the mean Wt
Wjt1 = reshape(Wjt,repbins,[])';
% get mean Wjt, i.e. Wt
Wt = (mean(Wjt1))';
% get free firing rate
qt = rt ./ Wt;
if(~Args.FrameSplit)
	% if frameavg was specified, average qt within a frame
	qt1 = reshape(qt,nbins,[]);
	qt2 = mean(qt1);
	% expand qt2 to get constant qt at time resolution of rt
	qt3 = repmat(qt2,nbins,1);
	% reshape back to vector
	qt = reshape(qt3,[],1);
elseif(Args.FrameSplit>1)
% 	% if FrameSplit was specified, split the frame into n parts and then
%     % average qt within each part
%     splitbins = nbins/Args.FrameSplit;
% 	qt1 = reshape(qt,splitbins,[]);
% 	qt2 = mean(qt1);
% 	% expand qt2 to get constant qt at time resolution of rt
% 	qt3 = repmat(qt2,splitbins,1);
% 	% reshape back to vector
% 	qt = reshape(qt3,[],1);

    splitbins = roundup(nbins,Args.FrameSplit)/Args.FrameSplit;

    %create first multiplier
    mat1 = ones(splitbins, 2);
    mat1(:, 2)= 0:(splitbins-1);

    %create second multiplier
    a = ones(1, Args.FrameSplit*nframes);
    b = 0:(Args.FrameSplit*nframes-1);
    c1 = [0:(nframes-1)]';
    c2 = repmat(c1, 1, Args.FrameSplit);
    c = reshape(c2', [], 1);
    difference = Args.FrameSplit*splitbins - nbins;
    d = a+splitbins*b-difference*(c');

    mat2 = ones(2, Args.FrameSplit*nframes);
    mat2(1,:) = d;

    %create index matrix
    idx = mat1*mat2;

    %include NaN
    nanmat1 = zeros(splitbins, Args.FrameSplit);
    if(nbins<(splitbins*Args.FrameSplit))
        nanmat1((nbins+1):(splitbins*Args.FrameSplit)) = NaN;
    end
    nanmat2=repmat(nanmat1, 1, nframes);
    idx = idx+ nanmat2;

    qt1 = nanindex(qt, idx);
    qt2 = nanmean(qt1);
    qt3 = repmat(qt2, splitbins, 1);
    %put back NaN
    qt4 = qt3 + nanmat2;
    qt = qt4(~isnan(qt4));
end

% return the first repbins+1 points of blimits as rtEdges (since the bins
% used for the old style histcie is expected so we need to add one bin to
% the end). Moreover since the first bin might not be 0, we subtract the 
% value of the first bin
rtEdges = blimits(1:(repbins+1)) - blimits(1,1);

%round up an integer number to the nearest integer that is divisible by 'base'

function a = roundup(b, base)
c=b/base;
d=ceil(c);
a=d*base;


