function [JointEvents] = JointEventsAnalysis(JointEventsInfo,stimInfo,Args)

%%% This is a private function used in the construction of the jointevents
%%% object. It also calls a couple other functions to calculate jointevents
%%% durations, the probabilities and the spike count correlations within
%%% the overlap time.

overlap = [];
for ii = 1:length(JointEventsInfo)-1    
    if isempty(overlap)
        eventline1=[];
        for e = 1:length(JointEventsInfo(ii).StartEvent)
            eventline1 = [eventline1 JointEventsInfo(ii).StartEvent(e):JointEventsInfo(ii).EndEvent(e)];
        end
        eventline2=[];
        for e = 1:length(JointEventsInfo(ii+1).StartEvent)
            eventline2 = [eventline2 JointEventsInfo(ii+1).StartEvent(e):JointEventsInfo(ii+1).EndEvent(e)];
        end 
        overlap = intersect(eventline1,eventline2);
    else
        eventline=[];
        for e = 1:length(JointEventsInfo(ii+1).StartEvent)
            eventline = [eventline JointEventsInfo(ii+1).StartEvent(e):JointEventsInfo(ii+1).EndEvent(e)];
        end 
        overlap = intersect(eventline,overlap);
    end    
end

%%%%%% If there is no overlap between the two cells then skip the lengthy
%%%%%% calculations.
if isempty(overlap)    
    JointEvents.Joint_Overlap_Probability = 0;
    expected_prob=1;
    for ii = 1:length(JointEventsInfo)
        if isempty(JointEventsInfo(ii).EventProbability)
            expected_prob=0;
            theoretical_min_prob(ii)=0;
        else
            expected_prob = [expected_prob * JointEventsInfo(ii).EventProbability];
            theoretical_min_prob(ii) = JointEventsInfo(ii).EventProbability;
        end
    end
    JointEvents.Expected_Overlap = expected_prob;
    %%% Similarity Index is the actual joint probablily / by the
    %%% expected_prob %%%
    if JointEvents.Joint_Overlap_Probability >= JointEvents.Expected_Overlap
        JointEvents.Similarity_Index = (JointEvents.Joint_Overlap_Probability - JointEvents.Expected_Overlap) / (min(theoretical_min_prob) - JointEvents.Expected_Overlap);
    else
        JointEvents.Similarity_Index = (JointEvents.Joint_Overlap_Probability - JointEvents.Expected_Overlap) / JointEvents.Expected_Overlap;
    end   
    JointEvents.Joint_Overlap_Start=[];
    JointEvents.Joint_Overlap_End=[];
    JointEvents.Joint_Overlap_Durations=[];
    JointEvents.Spike_Count_Correlations=[];
    JointEvents.Spike_Count_Significances=[];
    JointEvents.Shift_Spike_Count_Correlations=[];
    JointEvents.Shift_Spike_Count_Significances=[];
    R=[];
    for i_r = 1:Args.numRand
        ind1 = (randperm(size(JointEventsInfo(1).FramePSTH,1)))';
        ind2 = (randperm(size(JointEventsInfo(2).FramePSTH,1)))';
        [r,p] = corrcoef(JointEventsInfo(1).FramePSTH(ind1),JointEventsInfo(2).FramePSTH(ind2));
        R = [R r(1,2)];
    end
    JointEvents.RandomCC_Values = R';    
    [r,p] = corrcoef(JointEventsInfo(1).FramePSTH,JointEventsInfo(2).FramePSTH);
    if r(1,2) > prctile(R,95) | r(1,2) < prctile(R,5)
        JointEvents.RandomCC_Values_Significance_Value = 0;
    else
        JointEvents.RandomCC_Values_Significance_Value = 1;
    end    
    JointEvents.Correlation_Coefficent_Value = r(1,2);
    JointEvents.Correlation_Significance_Value = p(1,2);
    
else
    
    %%% Joint overlap Probability %%%%
    JointEvents.Joint_Overlap_Probability = length(overlap)/ceil(stimInfo.data.catInfo.repetition_duration);
    %%% Expected Joint Probablility via the product of the the individual
    %%% probablilites.
    expected_prob=1;
    theoretical_min_prob=[];
    for ii = 1:length(JointEventsInfo)
        expected_prob = [expected_prob * JointEventsInfo(ii).EventProbability];
        theoretical_min_prob(ii) = JointEventsInfo(ii).EventProbability;
    end
    JointEvents.Expected_Overlap = expected_prob;
    %%% Similarity Index is the actual joint probablily / by the
    %%% expected_prob %%%
    if JointEvents.Joint_Overlap_Probability >= JointEvents.Expected_Overlap
        JointEvents.Similarity_Index = (JointEvents.Joint_Overlap_Probability - JointEvents.Expected_Overlap) / (min(theoretical_min_prob) - JointEvents.Expected_Overlap);
    else
        JointEvents.Similarity_Index = (JointEvents.Joint_Overlap_Probability - JointEvents.Expected_Overlap) / JointEvents.Expected_Overlap;
    end
    %%% Joint overlap Durations %%%%    
    overlap_index = find(diff(overlap)>1);    
    overlap = sort([overlap(1) overlap(overlap_index) overlap(overlap_index+1) overlap(end)]);    
    overlap_start = overlap(1:2:end);
    overlap_end = overlap(2:2:end);    
    %%% Store the Start and End Points for the Joint Event Overlap %%%
    JointEvents.Joint_Overlap_Start = overlap_start;
    JointEvents.Joint_Overlap_End = overlap_end;
    JointEvents.Joint_Overlap_Durations = diff([JointEvents.Joint_Overlap_Start ; JointEvents.Joint_Overlap_End]);
    
    %%% Creation of the rasters for the spike counts %%%%%%
    for ii = 1:length(JointEventsInfo)    
        for r = 1:stimInfo.data.catInfo.num_repetitions     
            %%%% Frames for the repetition %%%%
            frames = JointEventsInfo(ii).framepoints(r*stimInfo.data.catInfo.num_frames-stimInfo.data.catInfo.num_frames+1:(r*stimInfo.data.catInfo.num_frames)+1);
            %%%% Spiketimes for the given repetition %%%%%
            rep_spikes = JointEventsInfo(ii).spiketrain(find(JointEventsInfo(ii).spiketrain > frames(1) & JointEventsInfo(ii).spiketrain <= frames(end)));    
            if r > 1    
                rasters(ii).raster{r} = rep_spikes-time_adjust;
                time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
            else
                if frames(1)~= 1
                    rasters(ii).raster{r} = rep_spikes-frames(1);
                    time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
                else
                    rasters(ii).raster{r} = rep_spikes;
                    time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
                end
            end  
        end
    end   
    
    %%%% Spike Count Correlations via the spike counts across the
    %%%% repetitions %%%
    for ii = 1:length(JointEvents.Joint_Overlap_Start)
        for r = 1:stimInfo.data.catInfo.num_repetitions
            spike_counts_1(r) = length(find(rasters(1).raster{r}>=JointEvents.Joint_Overlap_Start(ii) & rasters(1).raster{r}<JointEvents.Joint_Overlap_End(ii)));
            spike_counts_2(r) = length(find(rasters(2).raster{r}>=JointEvents.Joint_Overlap_Start(ii) & rasters(2).raster{r}<JointEvents.Joint_Overlap_End(ii)));
        end
        warning off MATLAB:divideByZero
        [r,p] = corrcoef(spike_counts_1,spike_counts_2);
        JointEvents.Spike_Count_Correlations(ii) = r(1,2);
        JointEvents.Spike_Count_Significances(ii) = p(1,2);
        spike_counts_1 = circshift(spike_counts_1',1)';
        warning off MATLAB:divideByZero
        [r,p] = corrcoef(spike_counts_1,spike_counts_2);
        JointEvents.Shift_Spike_Count_Correlations(ii) = r(1,2);
        JointEvents.Shift_Spike_Count_Significances(ii) = p(1,2);
        spike_counts_1 = [];
        spike_counts_2 = [];
    end 
    
    %%%% Correlation Coefficents for the two PSTH's at the Frame Resolution %%%%%
    R=[];
    for i_r = 1:Args.numRand
        ind1 = (randperm(size(JointEventsInfo(1).FramePSTH,1)))';
        ind2 = (randperm(size(JointEventsInfo(2).FramePSTH,1)))';
        [r,p] = corrcoef(JointEventsInfo(1).FramePSTH(ind1),JointEventsInfo(2).FramePSTH(ind2));
        R = [R r(1,2)];
    end
    JointEvents.RandomCC_Values = R';    
    [r,p] = corrcoef(JointEventsInfo(1).FramePSTH,JointEventsInfo(2).FramePSTH);
    if r(1,2) > prctile(R,95) | r(1,2) < prctile(R,5)
        JointEvents.RandomCC_Values_Significance_Value = 0;
    else
        JointEvents.RandomCC_Values_Significance_Value = 1;
    end    
    JointEvents.Correlation_Coefficent_Value = r(1,2);
    JointEvents.Correlation_Significance_Value = p(1,2);
    
end % if isempty(overlap)

%%% Creation of the raster %%%%%%
rasters=[];
for r = 1:length(JointEventsInfo)    
    for rr = (1:stimInfo.data.catInfo.num_repetitions)*r
        frames = JointEventsInfo(r).framepoints((rr/r)*stimInfo.data.catInfo.num_frames-stimInfo.data.catInfo.num_frames+1:((rr/r)*stimInfo.data.catInfo.num_frames)+1);
        %%%% Spiketimes for the given repetition %%%%%
        rep_spikes = JointEventsInfo(r).spiketrain(find(JointEventsInfo(r).spiketrain > frames(1) & JointEventsInfo(r).spiketrain <= frames(end)));    
        if (rr/r) > 1    
            rast = rep_spikes-time_adjust;
            time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
        else
            if frames(1)~= 1
                rast = rep_spikes-frames(1);
                time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
            else
                rast = rep_spikes;
                time_adjust = frames(end); % This adjusts the spiketimes for each repetition of the stimulus
            end
        end         
        if rr>1
            rasters = concatenate(rasters,rast);
        else
            rasters = rast;
        end        
    end
end   

JointEvents.Rasters = [1]; %rasters;
JointEvents.NumReps = size(rasters,1);
JointEvents.PSTH = [1]; %[JointEventsInfo(1).PSTH';JointEventsInfo(2).PSTH'];
JointEvents.UpperThresholds = [JointEventsInfo(1).Thresholds.meanupperthreshold JointEventsInfo(2).Thresholds.meanupperthreshold];
JointEvents.LowerThresholds = [JointEventsInfo(1).Thresholds.meanlowerthreshold JointEventsInfo(2).Thresholds.meanlowerthreshold];
JointEvents.StartEvent = concatenate(JointEventsInfo(1).StartEvent,JointEventsInfo(2).StartEvent,'ColumnWise');
JointEvents.EndEvent = concatenate(JointEventsInfo(1).EndEvent,JointEventsInfo(2).EndEvent,'ColumnWise');