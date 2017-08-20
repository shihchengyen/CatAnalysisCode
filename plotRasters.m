function plotRasters(sptrain)
%plotRasters Create raster plot for repeated stimulus presentation
%   plotRasters(SPTRAIN) plots the spikes times contained in SPTRAIN,
%   which is expected to be a cell array with spike times for each
%   repetition in each row of the cell array.

% get number of repetitions
reps = length(sptrain);
pheight = 1/reps;
pvec = 0:pheight:1;
pindex = 1;
hold on
for i = 1:reps
	% get spike times. Convert to row vector
	sp = vecr(sptrain{i});
	y1 = pvec(pindex);
	pindex = pindex + 1;
    plot(sp,y1,'b.')
	% line([sp; sp],[y1 pvec(pindex)],'Color','b');
end
hold off
