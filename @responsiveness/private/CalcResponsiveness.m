function [R,r,p,r2,p2,r3,p3] = CalcResponsiveness(obj,Args)
%%% This function takes the adjusted spiketrain and the adjusted frame
%%% points vectors and calculates the responsiveness of the cell. The
%%% responsiveness value is calculated by doing a repetition by repetition
%%% spike count correlation for all possible repetition comparisons. The
%%% responsiveness value is based on the percentage of significant
%%% correlations for all possible comparisons. The spike counts are
%%% calculated at the frame resolution.

Spiketrain = obj.data.adjSpiketrain';
FramePoints = obj.data.adjFramePoints';
num_repetitions = size(obj.data.raster,2);
numFrames = (length(FramePoints)-1)/num_repetitions;
%%%%% Spike Count Matrix for each Repetitions %%%%
frameCounts = histcie(Spiketrain,FramePoints,'DropLast');
spike_matrix = reshape(frameCounts,numFrames,[]);

% turn off warning in corrcoef function
warning off MATLAB:divideByZero
% compute correlation coefficient
[ccr,ccp] = corrcoef(spike_matrix);
% get unique coefficients
uidx = logical(tril(ones(size(ccr)),-1));
cvals = ccr(uidx);
% test to see if this distribution is normal
[r,p] = lillietest(cvals);
% get unique p values
pvals = ccp(uidx);
% get length of unique values
luvals = num_repetitions*(num_repetitions-1)/2;
% select only values that are not nan
pv2 = pvals(~isnan(pvals));
R = sum(pv2<Args.Significance)/length(pv2);

% randomize to create control distribution 
% create random matrix
rsets = Args.RandSets;
rcoefs = zeros(luvals,rsets);
if(Args.RandTimeTrial)
    rnum = num_repetitions*numFrames;
    for idx = 1:rsets
        % randomize across time and repetitions
		rn = rand(rnum,1);
		[rns,rni] = sort(rn);
		rni2 = reshape(rni,numFrames,[]);
		sm2 = spike_matrix(rni2);
		% compute correlation coefficients
		ccr2 = corrcoef(sm2);
		rcoefs(:,idx) = ccr2(uidx);
    end
else
    rsize = [numFrames num_repetitions];
    idxmat = repmat((1:numFrames)',1,num_repetitions);
    for idx = 1:rsets
        % randomize across repetitions
        rn = rand(numFrames,num_repetitions);
        [rns,rni] = sort(rn,2);
        sm2 = spike_matrix(sub2ind(rsize, ...
            idxmat,rni));
		% compute correlation coefficients
		ccr2 = corrcoef(sm2);
		rcoefs(:,idx) = ccr2(uidx);
    end
end
warning on MATLAB:divideByZero

cvals2 = rcoefs(:);
% test to see if this distribution is normal
[r2,p2] = lillietest(cvals2);

% now test if the two distributions are diffferent
% if r1 or r2 is 1, then we have to use kstest2
if(r==1 || r2==1)
    [r3,p3] = kstest2(cvals,cvals2);
else
    [r3,p3] = ttest2(cvals,cvals2);
end

if(Args.PlotDist)
    subplot(3,1,1)
    plot(obj,'Fast');
    subplot(3,1,2)
    hist(cvals)
    xv1 = xlim;
    subplot(3,1,3)
    hist(cvals2)
    xv2 = xlim;
    xv3 = [xv1; xv2];
    xv4 = [min(xv3(:,1)) max(xv3(:,2))];
    subplot(3,1,2)
    xlim(xv4)
    title(sprintf('Lillietest H: %d P: %f, Disttest H: %d P: %f',r,p,r3,p3))
    subplot(3,1,3)
    xlim(xv4)
    title(sprintf('Lillietest H: %d P: %f',r2,p2))
    pause
end
