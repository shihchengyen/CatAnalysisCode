function [R] = FitFR(data,Args)
% Function used to get a distribution of R squared values for the
% fits of the FR data, called in the firingrate/plot function when the
% bootstrp function is called.

if Args.Repetitions
    data = (mean(data))';
    max_value = nanmax(data(:));
    min_value = nanmin(data(:));
    if Args.HistRate > 0
        bins = (min_value:Args.HistRate:max_value+Args.HistRate)';
    else
        bins = (min_value:(max_value-min_value)/Args.numBins:max_value)';
    end
    N = histcie(data,bins,'DropLast');bins(1)=[];
    non_zeros = find(N>0);
    N = N(non_zeros);
    bins = bins(non_zeros);
else
    data = reshape(data',size(data,1)*size(data,2),1);
    max_value = nanmax(data(:));
    min_value = nanmin(data(:));
    if Args.HistRate > 0
        bins = (min_value:Args.HistRate:max_value+Args.HistRate)';
    else
        bins = (min_value:(max_value-min_value)/Args.numBins:max_value)';
    end
    N = histcie(data,bins,'DropLast');bins(1)=[];
    non_zeros = find(N>0);
    N = N(non_zeros);
    bins = bins(non_zeros);
end

%% Fit types to be used
% FitName1 = 'exp1'; % Exp Fit
% FitName2 = 'power1'; % Power Fit
FitName1 = 'poly1'; % Exp Fit
FitName2 = fittype('a + b*log(x)'); % Power Fit
N=log(N);

%% Fit the data
if length(N)<2
    fprintf('Error in Length')
    R=[NaN NaN];    
else
    [fitobj,goodness,output] = fit(bins,N,FitName1); % Exp
    if Args.AdjR
        R(1) = goodness.adjrsquare;
    else
        R(1) = goodness.rsquare;
    end
    [fitobj,goodness,output] = fit(bins,N,FitName2); % Power
    if Args.AdjR
        R(2) = goodness.adjrsquare;
    else
        R(2) = goodness.rsquare;
    end
end