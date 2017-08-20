function [vals,xvals] = getFunctionValues(cf,varargin)
%REFRACTORY/getFunctionValues - Return recovery function values
%   [VALS,XVALS] = getFunctionValues(CF,VARARGIN) returns the 
%   values of the recovery function in VALS and the corresponding 
%   time values in XVALS. CF is a curve fit object returned by the 
%   fit function in the Curve-Fitting Toolbox. 
%
%   The optional arguments are:
%      binsize - time resolution at which the recovery function 
%                  is evaluated (default is 0.2 ms).
%      cfmax - final value to use to when looking for the point in 
%              the recovery function that approaches 1 (default is
%              0.99).
%      cfstep - time step in ms to take when looking for value in
%               the recovery function that exceeds cfmax (default
%               is 2 ms).
%
%   [vals,xvals] = getFunctionValues(cf,'binsize',0.2,'cfmax',0.99, ...
%      'cfstep',2);

% default values for optional arguments
Args = struct('binsize',0.2, ...
			  'cfmax',0.99, ... % make sure the final value of the fitted 
			  				...	% curve is larger than this value
			  'cfstep',2); 	% amount of time (ms) to step when looking for 
			  				% the value of the fitted curve that exceeds
			  				% cfmax

% get optional arguments
Args = getOptArgs(varargin,Args);

% check curve fit values
cfstart = 0;
cfend = Args.cfstep;
cfx = (cfstart:Args.binsize:cfend)';
cfy = feval(cf,cfx);
while(cfy(end)<Args.cfmax)
	cfstart = cfend + Args.binsize;
	cfend = cfstart + Args.cfstep;
	newcfx = (cfstart:Args.binsize:cfend)';
	cfy = [cfy; feval(cf,newcfx)];
	cfx = [cfx; newcfx];
end
% find first point that exceeds cfmax to limit redundancy in the recovery function
cfi = find(cfy>Args.cfmax);
vals = cfy(1:cfi(1)); 
xvals = cfx(1:cfi(1));
