function sptrain = readSurrogateBin(filename)
%readSurrogateBin Reads binary file containing surrogate spike trains
%   SPTRAIN = readSurrogateBin(FILENAME) reads the binary file FILENAME
%   and returns a cell array, SPTRAIN, which has the following format:
%   sptrain{set}{repetition}(spike). Spike times are in units of ms.

fid = fopen(filename,'r','ieee-le');
% read number of sets and number of reps
[a,n] = fread(fid,2,'int32');
if(n~=2)
	error('Unexpected file format!')
end
sets = a(1);
reps = a(2);
% allocate memory for sets
sptrain = cell(1,sets);
for i = 1:sets
	% allocate memory for reps
	sptrain{i} = cell(reps,1);
	for j = 1:reps
		% read number of spikes for set i, repetition j
		nst = fread(fid,1,'int32');
		% read spike times for set i, repetition j
		st = fread(fid,nst,'int32');
		% scale integer spike times back to double
		sptrain{i}{j} = st*0.1;
	end
end
fclose(fid);
