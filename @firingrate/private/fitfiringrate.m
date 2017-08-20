function Fit = fitfiringrate(data,Args,cwd)
% This private function will fit the firingrate data with specified
% function. It will also bootstrap the firingrate data based on the
% repetitions.

firingRate = data.firingRate;

% Toss out the zeros in the firingRate and replace with NaN's
if Args.TossZeros
    firingRate(find(firingRate==0)) = NaN;
end

%%%% Histogram the firingrate values %%%%
max_value = nanmax(firingRate);
min_value = nanmin(firingRate);
if Args.HistRate > 0
    bins = (min_value:Args.HistRate:max_value+Args.HistRate)';
else
    bins = (min_value:(max_value-min_value)/Args.numBins:max_value)';
end
[N] = histcie(firingRate,bins,'DropLast');bins(1)=[];

%%%% Normalize to a Probablility %%%%
if Args.Probability
    tot = sum(N);
    N = N./tot;
end

%%% Values to be filled in
Fit.N = N;
Fit.bins = bins;
Fit.Surrogate_Mean_R_Values=[];
Fit.Surrogate_H_Value=[];
Fit.Surrogate_P_Value=[];
Fit.BootStrap_R_Values=[];
Fit.BootStrap_Mean_R_Values=[];
Fit.BootStrap_H_Value=[];
Fit.BootStrap_P_Value=[];
Fit.Firing_Rate_R_Values=[];

if Args.Surrogates
    
    %%% Load the Surrogate Spike data
    cd(cwd); %Change to the data directory
    files = nptDir('frames*');
    if ~isempty(files)
        cd ../..; st = stiminfo('auto');
        cd(cwd); % Change to the data directory
        [R,H,S_Inv] = fitsurrogatefr(files,st,Args);
        Surrogate_Inv = S_Inv;
        R_Values = R;
        Fit.Surrogate_Mean_R_Values = mean(R);  
        [H,P] = ttest2(R(:,1),R(:,2),.01);
        Fit.Surrogate_H_Value = H;
        Fit.Surrogate_P_Value = P;        
    else
        fprintf('No Surrogates Found')
    end
    
elseif Args.BootStrap
    
    %%% Run the Bootstrap function to get R value statistics
    [R] = bootstrp(Args.NumBoots,'fitfr',data.spike_matrix,Args); % First Column values are exp and second are
    %% for the power fits for the r-squared values.
    Fit.BootStrap_R_Values = R;
    Fit.BootStrap_Mean_R_Values = mean(R);
    [H,P] = ttest2(R(:,1),R(:,2),.01);
    Fit.BootStrap_H_Value = H;
    Fit.BootStrap_P_Value = P;        
    
end

% Use this fit function to find the best fit
%[rmse,sse,adjrs,list]=jcurvefit(N,bins);

% Fit the Raw Data tossing the bins with zero spikes
non_zeros = find(N>0);
N = N(non_zeros); 
bins = bins(non_zeros);

% Matlab Function for Fitting
% Fit.FitName1 = 'exp1'; % Exp Fit
% Fit.FitName2 = 'power1'; % Power Fit fittype('a + b*log(x)'); % Power Fit

% Our Functions for Fitting
Fit.FitName1 = 'poly1'; % Exp Fit
Fit.FitName2 = fittype('a + b*log(x)'); % Power Fit
N=log(N);

[fitobj1,goodness1,output1] = fit(bins,N,Fit.FitName1);
if Args.AdjR
    Fit.Firing_Rate_R_Values(1) = goodness1.adjrsquare;
else
    Fit.Firing_Rate_R_Values(1) = goodness1.rsquare;
end
[fitobj2,goodness2,output2] = fit(bins,N,Fit.FitName2);
if Args.AdjR
    Fit.Firing_Rate_R_Values(2) = goodness2.adjrsquare;
else
    Fit.Firing_Rate_R_Values(2) = goodness2.rsquare;
end

return
figure
subplot(1,2,1)
plot(bins,N,'*-'); axis tight; hold on
plot(fitobj1); 
ylabel('Log(Counts)')
xlabel('Firing Rate')
title(['EXP: ',num2str(Fit.Firing_Rate_R_Values(1))])
subplot(1,2,2)
plot(bins,N,'*-'); axis tight; hold on
plot(fitobj2); 
ylabel('Log(Counts)')
xlabel('Firing Rate')
title(['POWER: ',num2str(Fit.Firing_Rate_R_Values(2))])