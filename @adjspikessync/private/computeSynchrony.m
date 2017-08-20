function xc = computeSynchrony(sp1,sp2,nsbins,reps,maxlag)

spprod = sp1 .* sp2;
% now take the sum within a frame
sppsum = sum(spprod);
% now reshape so we can average over repetitions
spps1 = reshape(sppsum,[],reps);
xc = sum(spps1,2);
% now do the same for different shifts
for(idx=1:maxlag)                  
	% take product of both spike trains after shift
	nsb1 = 1:(nsbins-idx);
	nsb2 = (1+idx):nsbins;
	spprod = sp1(nsb1,:) .* sp2(nsb2,:);
	sppsum = sum(spprod);
	spps1 = reshape(sppsum,[],reps);
	% add sum of products to xc
	xc = xc + sum(spps1,2);
	% take product of both spike trains after opposite shift
	spprod = sp1(nsb2,:) .* sp2(nsb1,:);
	sppsum = sum(spprod);
	spps1 = reshape(sppsum,[],reps);
	% add sum of products to xc
	xc = xc + sum(spps1,2);
end % for(idx=1:maxlag)
