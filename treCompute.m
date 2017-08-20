function treCompute(framebins,filename,sgname)
%treCompute Compute Temporal Reliability Entropy
%   treCompute(FRAMEBINS,FILENAME) computes the Temporal Reliability
%   Entropy for the surrogate data. FRAMEBINS specifies the number of
%   bins each frame is broken up into. The resulting entropy measurements
%   are saved in FILENAME.mat.

if(ischar(framebins))
	framebins = str2num(framebins);
end

% load data
mat = load('refractory.mat');
rf = struct(mat.rf);

surrFiles = 10;
surrSetsPerFile = 100;
surrTotal = surrFiles * surrSetsPerFile;

% get frame boundaries
rtlength = length(rf.data.rtEdges{1});
% take out last value in rtEdges so there will be an even number of points
rt1 = reshape(rf.data.rtEdges{1}(1:(rtlength-1)),rf.data.qtframebins,[]);
% get the last value which will be used a couple of times
rtlast = rf.data.rtEdges{1}(rtlength);
% number of columns in rt1 is number of frames
numframes = size(rt1,2);
% grab the first row which should be the frame limits 
frameLimits = [rt1(1,:) rtlast];
% get difference between frame limits
dfL = diff(frameLimits);
% get bin sizes per frame
binsize = dfL/framebins;
% create matrix of bin limits
mat1 = tril(ones(framebins));
mat2 = [frameLimits(1:numframes); repmat(binsize,framebins-1,1)];
binLimits = mat1 * mat2;
% append time of last frame
histEdges = [reshape(binLimits,[],1); rtlast];

% preallocate memory
entropy = zeros(surrTotal,numframes);

% initialize loop variables
index = 1;

% loop over surrogates
for i = 1:surrFiles
	% load surrogates
	sptrain = readSurrogateBin([sgname num2str(i) '.bin']);
	% loop over sets in each file
	for j = 1:surrSetsPerFile
		fprintf('Computing surrogate %d\n',index);
		% binSC dimensions should be (numframes * framebins) x repetitions
		binSC = histcie(cell2array(sptrain{j}),histEdges,'DropLast');
		% reshape into framebins x (numframes * repetition) so we can sort
		% bin ranks within each frame
		frameSC = reshape(binSC,framebins,[]);
		% code adapted from Jean-Philippe's TRESMeasure.m
		% first sort gives indices
		[fy,frank1] = sort(frameSC);
		% second sort uses indices to get ranks
		[fy,frank1] = sort(frank1);
		% flip updown to control for equal ranks
		[fy,frank2] = sort(flipud(frameSC));
		[fy,frank2] = sort(frank2);
		frank2 = flipud(frank2);
		% get final ranks
		m_rank = (frank1 + frank2)/2;
		% reshape into (framebins * numframes) x repetitions and then transpose
		% to get repetitions x (framebins * numframes)
		m_rank = reshape(m_rank,framebins*numframes,[])';
		
		% now compute frequencies of ranks for each subbin
		% this is a little trick to compute faster	
		npossranks = (2 * framebins) - 1;
		m_rank=repmat(m_rank,[1 1 npossranks]);
		v_possiblerank=zeros([1 1 npossranks]);
		v_possiblerank(1,1,:)=(2:2*framebins)/2;
		m_possiblerank=repmat(v_possiblerank,[rf.data.repetitions (framebins * numframes) 1]);	
		m_rank=(m_rank==m_possiblerank);	
		% m_cases becomes 1 x (framebins * numframes) x npossranks
		m_cases=mean(m_rank,1);
		% m_cases now becomes (framebins * numframes) x npossranks
		m_cases=squeeze(m_cases);
		% we replace all the zero values by non zero values	so we can take log
		% zero values will be zero'ed out since we are doing element-wise 
		% multiplication below
		m_ca=m_cases+(m_cases==0);
		% do element-wise multiplication and take mean across rows first 
		% to retain frame order
		m_ca2 = mean(m_cases.*log(m_ca),2);
		% reshape so we can get the final mean for each frame 
		entropy(index,:)=-mean(reshape(m_ca2,framebins,[]));
		index = index + 1;
	end
end

% save data
save(filename,'entropy')
