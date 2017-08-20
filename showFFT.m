function [MX,f] = nptSpikeTimesFFT(sptimes,varargin)
%nptSpikeTimesFFT Computes FFT magnitude given spike times
%   [MX,F] = nptSpikeTimesFFT(SPTIMES,VARARGIN) plots the magnitude of 
%   the FFT of the spike train created from the spike times contained 
%   in the SPTIMES vector. If SPTIMES is a matrix, the FFT is taken 
%   columnwise. This function computes the FFT in segments to increase 
%   the bin size in the frequency domain. The Welch window is used to 
%   window data in the time domain.
%
%   The optional input arguments are:
%      'BinSize' - followed by number specifying size in ms to use to 
%                  create spike train (default 4).
%      'Duration' - followed by number specifying duration of repetition
%                   in ms (default: 30000).
%      'Segment' - followed by number specifying duration of a segment
%                  in ms (default: 4000).
%      'Overlap' - followed by number specifying fraction of segment 
%                  windows that overlap between consecutive segments 
%                  (default: 0.5).

Args = struct('BinSize',4,'Duration',30000,'Segment',4000, ...
			'Overlap',0.5);
Args = getOptArgs(varargin,Args);

% get number of columns in sptimes
ntrains = size(sptimes,2);
% create bins for histcie
bins = 0:Args.BinSize:Args.Duration;
% get spike counts
spcounts = histcie(sptimes,bins,'DropLast');
% get number of rows in spcounts
spcrows = size(spcounts,1);

% figure out how many segments we can fit into Args.Duration
% first subtract 1 segment length and then see how many non-overlapped 
% windows we can fit. We will add 1 if we need total number of segments.
noOverlapLength = Args.Segment*(1-Args.Overlap);
n1segments = floor( (Args.Duration-Args.Segment) / noOverlapLength );
nsegments = n1segments + 1;

% create matrices that will be multiplied together to get indices for
% spcounts to pass to PlotFFT
% first determine number of bins in each segment
segmentBins = Args.Segment/Args.BinSize;
mat1 = tril(ones(segmentBins));
% create starting indicies for each segment
startind = ((0:n1segments) * (noOverlapLength/Args.BinSize)) + 1;
b = repmat(0:(ntrains-1),nsegments,1);
c = reshape(b,[],1);
d = repmat(startind',ntrains,1);
e = [c d] * [spcrows; 1];
% get total number of segments
totalnsegments = nsegments * ntrains;
% mat 2 will contain the starting indices for each segment in the first
% row
mat2 = [e'; ones(segmentBins-1,totalnsegments)];
spind = mat1 * mat2;

% window data with Welch window before computing FFT
[MX,f] = nptFFTMag(spcounts(spind).*repmat(welch(segmentBins),1,totalnsegments), ...
	1000/Args.BinSize);
