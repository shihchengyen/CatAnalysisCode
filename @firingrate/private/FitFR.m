function [R] = FitFR(data)
% Function used to get a distribution of R squared values for the
% fits of the FR data, called in the firingrate/plot function when the
% bootstrp function is called.

%% Create the edges for the histogram
binedges = (min(data):(max(data)-min(data))/100:max(data))';
%% Count the number of FR values
[N] = histcie(data,binedges,'DataCols','DropLast');
binedges = binedges(2:end);
non_zeros = find(N>0);
N = log(N(non_zeros)); % Take out the zero values for the fit and takes the log
binedges = binedges(non_zeros); % Take out the zero values for the fit
%% Fit types to be used
FitName1 = 'poly1';
FitName2 = fittype('a + b*log(x)');
%% Fit the data 
[fitobj,goodness,output] = fit(binedges,N,FitName1); % Exp
R(1) = goodness.adjrsquare;
[fitobj,goodness,output] = fit(binedges,N,FitName2); % Power
R(2) = goodness.adjrsquare;