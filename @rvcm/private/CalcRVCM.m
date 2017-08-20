function data = CalcRVCM(obj,Args)

pdir = pwd;
data.numSets = 1;
data.setNames{1} = pdir;

ST = obj.data.adjSpiketrain';
FP = obj.data.adjFramePoints';
numReps = size(obj.data.raster,2);
numFrames = (length(FP)-1)/numReps;
%%%%% Spike Count Matrix for each Repetitions %%%%
FC = histcie(ST,FP,'DropLast')';
SM = reshape(FC,numFrames,[]);

% allocate memory outside the for-loop
data.rvcmc = zeros(numReps,1);
data.rvcmp = data.rvcmc;
data.rvcmcn = data.rvcmc;

if(Args.Original)
    for a = 1:numReps
        if a == 1
            sp2=mean(SM(:,2:5),2);
        elseif a == 2
            sp2=mean(SM(:,[1 3:5]),2);
        elseif a == numReps-1
            sp2=mean(SM(:,[(numReps-4):(numReps-2) numReps]),2);
        elseif a == numReps
            sp2=mean(SM(:,[(numReps-4):(numReps-1)]),2);
        else
            sp2=mean(SM(:,[a-2:a-1 a+1:a+2]),2);
        end
        [ccr,ccp] = corrcoef(SM(:,a),sp2);
        data.rvcmc(a)=ccr(1,2);
        data.rvcmp(a)=ccp(1,2);    
    end
else    
    for a = 1:numReps
        if a == 1
            sp1=mean(SM(:,1:2),2);
            sp2=mean(SM(:,3:4),2);
        elseif a == 2
            sp1=mean([SM(:,1) SM(:,3)],2);
            sp2=mean(SM(:,4:5),2);
        elseif a == numReps-1
            sp1=mean(SM(:,numReps-4:numReps-3),2);
            sp2=mean([SM(:,numReps-2) SM(:,numReps)],2);
        elseif a == numReps
            sp1=mean(SM(:,numReps-4:numReps-3),2);
            sp2=mean(SM(:,numReps-2:numReps-1),2);
        else
            sp1=mean(SM(:,a-2:a-1),2);
            sp2=mean(SM(:,a+1:a+2),2);
        end
        [ccr,ccp] = corrcoef(sp1,sp2);
        data.rvcmc(a)=ccr(1,2);
        data.rvcmp(a)=ccp(1,2);    
    end
end
% subtract all correlation coefficients from the first 
data.rvcmcn = data.rvcmc-data.rvcmc(1);
% compute correlations between all repetitions
allccr = corrcoef(SM);
% save the unique correlations
ctrilind = logical(tril(ones(numReps,numReps),-1));
data.repcc = allccr(ctrilind);

return

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
