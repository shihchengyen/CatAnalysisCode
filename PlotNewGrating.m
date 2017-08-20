function PlotNewGrating(experiment_name,site,session,group_number,cell_number,colorcode)

% This function will load the data structure for the gratings using the
% current data hierarchy of day/site/session/group/cell
% After the data is loaded the grating response is reordered to display

%Change directory to load files
% eval(['cd c:/data/' experiment_name '/site' site '/session' session '/group' group_number '/cell' cell_number]);
cmdstr = ['cd /Users/syen/Documents/ShihCheng/Data/Neural/Cat/catdata/' experiment_name '/site' site '/session' session '/group' group_number '/cell' cell_number];
eval(cmdstr);
load data

%%% Get the information to plot %%%%
raster_order = data.stimulus_info.raster_order;
frame_duration = data.stimulus_info.frame_duration;
raster = data.cell_info.raster;
num_repetitions = data.stimulus_info.num_blocks;
repetition_duration = frame_duration*data.stimulus_info.num_frames;
if session == '31' %% Hack for this session of a431
    drifting_order = data.stimulus_info.drifting_order;
    num_cycles = max(drifting_order{1});
    cycle_number = drifting_order{1}(1);
    repetition_vector = 0:frame_duration:frame_duration*length(drifting_order{2});
end

if session == '31' %% Hack for this session of a431
    %%%%%% Adjust the spiketimes %%%%%%%
    for r = 1:length(raster)    
        spiketrain = raster{r};
        index = find(drifting_order{r} == cycle_number); index = index(1);
        adj_spikes = spiketrain(find(spiketrain >= repetition_vector(index) & spiketrain < repetition_vector(end)));
        if ~isempty(adj_spikes)
            adj_spikes = adj_spikes - repetition_vector(index);
            raster{r} = adj_spikes;
        else
            raster{r} = [];
        end    
    end
end

%%%%% reorder the grating responses %%%%%%%%
for g = 1:length(raster_order)
    grating_raster{g} = raster{raster_order(g)};
end

clf
if exist('colorcode')
    %% Plot a Rastergram of the spikes %%%%%
    for r = 1:length(grating_raster)
        if ~isempty(grating_raster{r})
            hold on
            plot(grating_raster{r},r,['.' colorcode])
        end
    end
    axis ij
    ylim([0 length(grating_raster)+1])
    xlim([0 repetition_duration+1])
    xlabel('Milliseconds')
    ylabel('Repetitions')
end
hold off
