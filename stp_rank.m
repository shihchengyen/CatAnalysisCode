function [m_cases, s_entropy, s_score, s_std]=stp_rank(m_matspike,s_nbbins,s_nbsurrogate,s_graph)

%I have this new measure of reliability
%it's quite fast
%The proceed the following way:
%1. I divide a frame into a couple of subbins (let's say 3, non-overlapping)
%2. For each subbin I count the number of spikes, and then I rank the bin (the first bin has the 3-rd largest number of spikes, and so on ...)
%3. The assumption is that if there is some temporal organization, then it should more or less always be the same subbin that has the largest number of spikes, it should always be the same one that has the less spikes, and so on ...
%4. In fact, I can count for each bin the number of times it has been ranked number 1, number 2, number 3, and so on ...
%since I may have ex-aequos, I also end up with fractions . like 2.5 , so if I have n subbins I have 2n-1 possible rankings
%any way, I end up with a matrix that is n x (2n-1)
%
%If the temporal precision is perfect then in this matrix, only n boxes have non-zero values (60 = the number of trials)
%so, I can probably compute an entropy measure here, that goes from the minimum, where any subbin can receive any ranking, to the extreme case above.
%
%I should compute a small function that does exactly that :
%that is : you give a matrix of spikes., that is p reps with q time steps where you have either a 1 or a 0.
%then you give a number of subbins that divides up this number q, you specify how many surrogates you want to compute
%and it gives you the nx(2n-1) matrix + an entropy, plus the proportion of surrogates that have a worse entropy (i.e less order)
%
%that would be stp_rank (for spike temporal precision - rank)
%
% m_matspike must be trial x sample
%
% [m_cases, s_entropy, s_score, s_std]=stp_rank(m_matspike,s_nbbins,s_nbsurrogate,s_graph)

[s_trial,s_time]=size(m_matspike);
if (rem(s_time,s_nbbins)>0)
   error('Oups, slight mistake ... the number of bins is not a divider of the number of samples in your time interval');
end;

s_binsize=s_time/s_nbbins;
m_data=zeros(s_trial,s_nbbins);
for s_b=1:s_nbbins
   m_sub=m_matspike(:,1+(s_b-1)*s_binsize:s_b*s_binsize);
   v_nbspike=sum(m_sub');
   m_data(:,s_b)=v_nbspike';
end;

[m_y,m_ranka]=sort(m_data,2);
[m_y,m_ranka]=sort(m_ranka,2);

[m_y,m_rankb]=sort(fliplr(m_data),2);
[m_y,m_rankb]=sort(m_rankb,2);

m_rankb=fliplr(m_rankb);

m_rank=(m_ranka+m_rankb)/2;
m_ranko=m_rank;
m_rank=repmat(m_rank,[1 1 2*s_nbbins-1]); % this is a little trick to compute faster

v_possiblerank=zeros([1 1 2*s_nbbins-1]);
v_possiblerank(1,1,:)=(2:2*s_nbbins)/2;

m_possiblerank=repmat(v_possiblerank,[s_trial s_nbbins 1]);

m_rank=(m_rank==m_possiblerank);
m_cases=mean(m_rank,1);
m_cases=squeeze(m_cases); % this is the matrix with the frequencies of ranking for each subbins

if (s_graph)
   figure;
	bar3(m_cases);
	colorbar;
   figure;
   imagesc(m_matspike);
end;
m_ca=m_cases+(m_cases==0); % we replace all the zero values by non zero values
s_entropy=-mean(mean(m_cases.*log(m_ca)));



% and now the same thing for the surrogates

%m_data=rand(s_trial*s_nbsurrogate,s_nbbins);

%[m_y,m_ranka]=sort(m_data,2);
%[m_y,m_ranka]=sort(m_ranka,2);

%[m_y,m_rankb]=sort(fliplr(m_data),2);
%[m_y,m_rankb]=sort(m_rankb,2);

%m_rankb=fliplr(m_rankb);

%m_rank=(m_ranka+m_rankb)/2;

m_h1=repmat(m_ranko,[1 1 s_nbsurrogate]);
m_h2=shiftdim(m_h1,2);
m_h3=reshape(m_h2,[s_nbsurrogate*s_trial s_nbbins]);

[m_v,m_h4]=sort(rand(s_nbsurrogate*s_trial,s_nbbins),2);
m_h4=(m_h4-1)*(s_nbsurrogate*s_trial);
m_h4=m_h4+repmat((1:s_nbsurrogate*s_trial)',1,s_nbbins);
m_rank=m_h3(m_h4);

m_rank=repmat(m_rank,[1 1 2*s_nbbins-1]); % this is a little trick to compute faster

m_possiblerank=repmat(v_possiblerank,[s_trial*s_nbsurrogate s_nbbins 1]);

m_rank=(m_rank==m_possiblerank);

m_rank=reshape(m_rank,[s_nbsurrogate s_trial s_nbbins 2*s_nbbins-1]);
m_rank=mean(m_rank,2);

m_cases=m_rank+(m_rank==0).*ones(size(m_rank)); % this is to avoid taking log(0)

m_prob=-m_cases.*log(m_cases); % this is still of dim 4 - s_nbsurrogate x 1 x s_nbbins x 2*...
m_prob=mean(m_prob,4);
m_prob=mean(m_prob,3);

m_prob=squeeze(m_prob);
s_score=length(find(m_prob>s_entropy));

s_std=-(s_entropy-mean(m_prob))/std(m_prob);


