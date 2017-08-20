function y = gaussian(x,xdata)
%GAUSSIAN Returns gaussian values
%   Y = GAUSSIAN(X,XDATA) where:
%      X(1) is the mean,
%      X(2) is the standard deviation,
%      X(3) is the scaling factor, and
%      X(4) is the pedestal.

y = x(4) + x(3) * exp(-((xdata-x(1)).^2)/(2*x(2)^2));
