function eventCompute(useThresholds,filename)

if(ischar(useThresholds))
	useThresholds = str2num(useThresholds);
end

nSurrogateFiles = 10;
nSurrogateSets = 100;

% load data file for cell
matfile = nptDir('*g*c*.mat');
load(matfile(1).name)
repetition_duration = data.stimulus_info.frame_duration*(data.stimulus_info.end_frame-data.stimulus_info.start_frame+1);
frame_duration = data.stimulus_info.frame_duration;
event_ind = 1;
bin_sizes = data.cell_info.events(event_ind).binsize;
perc_reps = data.cell_info.events(event_ind).percentage_reps;
thresholds = data.cell_info.events(event_ind).thresholds;

% initialize output
event.eventIntervals = cell(nSurrogateFiles * nSurrogateSets,1);
event.substimreps = event.eventIntervals;
index = 1;

% loop over surrogate files
for i = 1:nSurrogateFiles
	% load surrogates
	sptrain = readSurrogateBin(['framesg' num2str(i) '.bin']);
	% loop over surrogate sets
	for j = 1:nSurrogateSets
		fprintf('Computing surrogate %d\n',index);
        if(useThresholds)
            ev = EventsAnalysis(sptrain{j},bin_sizes(1),repetition_duration, ...
                frame_duration,perc_reps,thresholds);
        else
            ev = EventsAnalysis(sptrain{j},bin_sizes(1),repetition_duration, ...
                frame_duration,perc_reps);
        end
		event.eventIntervals{index} = ev.event_interval_durations;
        event.substimreps{index} = ev.substimreps;
		index = index + 1;
	end
end

save(filename,'event');
