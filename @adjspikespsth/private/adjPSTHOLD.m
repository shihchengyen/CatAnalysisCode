function [p_psth,plot_vector] = adjPSTH(obj,stimInfo,binsize,slide_duration)
% PlotPSTH(binsize,slide_duration)
% will load the adjspikes and stiminfo objects to create and then display
% the post stimulus time histogram (PSTH) 


%%%%%%%%% Load information from the AdjRaster and StimInfo Objects %%%%%%%%%%%%%%%%%%%%%
spiketrain = obj.data.adjSpiketrain;
frames_vector = obj.data.adjFramePoints;

frame_duration = stimInfo.data.catInfo.frame_duration;
num_frames = stimInfo.data.catInfo.num_frames;
repetition_duration = stimInfo.data.catInfo.repetition_duration;
num_repetitions = stimInfo.data.catInfo.num_repetitions;

%%% Create empty space in memory %%%%%
if ~isnumeric(binsize)
    binsize = frame_duration;
    if isnumeric(slide_duration)
        spike_matrix = zeros(num_repetitions,length(0:1:repetition_duration));
        plot_vector = 0:1:repetition_duration;
    else
        spike_matrix = zeros(num_repetitions,length(0:frame_duration:repetition_duration));
        plot_vector = 0:frame_duration:repetition_duration;
    end
else
    if isnumeric(slide_duration)
        spike_matrix = zeros(num_repetitions,length(0:1:repetition_duration));
        plot_vector = 0:1:repetition_duration;
    else
        spike_matrix = zeros(num_repetitions,length(0:binsize:repetition_duration));
        plot_vector = 0:binsize:repetition_duration;
    end
end

%%%%% Create a Spike Count Matrix %%%%%%%%
for r = 1:num_repetitions    
    psth_vector = frames_vector(1:num_frames+1);
    frames_vector(1:num_frames)=[];
    if isempty(binsize)
        binsize = frame_duration;
        if isnumeric(slide_duration)
            psth_vector = psth_vector(1):1:psth_vector(end);
        end
        spike_matrix(r,1:length(psth_vector)) = histc(spiketrain,psth_vector);    
    else
        if isnumeric(slide_duration)
            psth_vector = psth_vector(1):1:psth_vector(end);
        else
            psth_vector = psth_vector(1):binsize:psth_vector(end);
        end
        spike_matrix(r,1:length(psth_vector)) = histc(spiketrain,psth_vector);
    end        
end

%%%%%%%%%%%%% PSTH based on the average spike counts per bin in time %%%%%%%
if isnumeric(slide_duration)    
    if ~isempty(slide_duration)        
        slide_duration = round((binsize-1)/2); % This is the window of time for the sliding window.                
        for s = 1:length(psth_vector)            
            m_psth = [];
            if s < slide_duration+1
                m_psth = spike_matrix(:,s:s+slide_duration);
                dur=size(m_psth,2);
                m_psth = sum(m_psth,2);
                p_psth(s) = mean(m_psth)/(dur/1000);
            elseif length(psth_vector) - s  < slide_duration
                m_psth = spike_matrix(:,s-slide_duration:s);
                dur=size(m_psth,2);
                m_psth = sum(m_psth,2);
                p_psth(s) = mean(m_psth)/(dur/1000); 
            else
                m_psth = spike_matrix(:,s-slide_duration:s+slide_duration);
                dur=size(m_psth,2);
                m_psth = sum(m_psth,2);
                p_psth(s) = mean(m_psth)/(dur/1000);    
            end
        end        
    else
        p_psth = mean(spike_matrix)/(binsize/1000);
    end
else
    p_psth = mean(spike_matrix)/(binsize/1000);
end

if length(plot_vector) ~= length(p_psth)
    p_psth = p_psth(1:length(plot_vector));
end

%%%% Plot the PSTH %%%%%
%plot(plot_vector,p_psth)
