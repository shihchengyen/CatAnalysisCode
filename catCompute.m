function catCompute(sets,filename)
%catCompute Computes various statistics for cat data
%   catCompute(sets,filename)

if(ischar(sets))
	sets = str2num(sets);
end

% load refractory data
mat = load('refractory.mat');
rf = struct(mat.rf);

% get frame boundaries
rtlength = length(rf.data.rtEdges{1});
% take out last value in rtEdges so there will be an even number of points
rt1 = reshape(rf.data.rtEdges{1}(1:(rtlength-1)),rf.data.qtframebins,[]);
% grab the first row which should be the frame limits and add the last
% point in rtEdges
fbins = [rt1(1,:)'; rf.data.rtEdges{1}(rtlength)];
% number of frames should be number of columns in rt1
numframes = size(rt1,2);
% get number of repetitions
reps = rf.data.repetitions;

% allocate memory
sptrain = cell(1,sets);
sg.scmean = zeros(numframes,sets);
sg.scstd = zeros(numframes,sets);

% generate surrogates
for i = 1:sets
	fprintf('Set %d\n',i);
	% generate surrogate spike trains, each repetition is in a column
	sptrain{i} = getSurrogateSpikes(rf.data,1);
	% get spike counts for each frame of sptrain, one column for each rep
	sc = histcie(cell2array(sptrain{i}),fbins,'DropLast');
	% reshape into repetitions and then transpose so we can take the mean
	% and std easily
	sc1 = reshape(sc,numframes,[])';
	% store mean and std in columns
	sg.scmean(:,i) = mean(sc1)';
	sg.scstd(:,i) = std(sc1)';
end

% save surrogates and statistics for surrogates
fid = fopen([filename '.bin'],'w','ieee-le');
fwrite(fid,[sets reps],'int32');
for i = 1:sets
	for j = 1:reps
		st = sptrain{i}{j};
		% write number of spikes for set i, repetition j
		fwrite(fid,size(st,1),'int32');
		% write spike times for set i, repetition j
		% multiply by 10 so we can convert tenth of a ms to integer
		fwrite(fid,st*10,'int32');
	end
end
fclose(fid);
save([filename 'FF'],'sg');
