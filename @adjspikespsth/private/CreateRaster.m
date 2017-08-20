function raster = CreateRaster(adjusted_frames_vector,adjusted_spiketrain,num_repetitions,num_frames)

%%% Creation of the raster %%%%%%
for ii = 1:num_repetitions     
    %%%% Frames for the repetition %%%%
    frames = adjusted_frames_vector(ii*num_frames-num_frames+1:(ii*num_frames)+1);
    %%%% Spiketimes for the given repetition %%%%%
    rep_spikes = adjusted_spiketrain(find(adjusted_spiketrain > frames(1) & adjusted_spiketrain <= frames(end)));    
    if ii > 1    
        raster.spikes = concatenate(raster.spikes,rep_spikes-time_adjust);
        time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
    else
        if frames(1)~= 1
            raster.spikes(ii,:) = rep_spikes-frames(1);
            time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
        else
            raster.spikes(ii,:) = rep_spikes;
            time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
        end
    end    
end