function entropy = getTREntropy(frameSC,numframes);

[framebins,frameSCcols] = size(frameSC);
repetitions = frameSCcols/numframes;
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
m_possiblerank=repmat(v_possiblerank,[repetitions (framebins * numframes) 1]);	
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
entropy=-mean(reshape(m_ca2,framebins,[]));
