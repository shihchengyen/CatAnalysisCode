function [events] = EventsAnalysis(obj,stimInfo,Args)
% This function is the main private function used to create the information
% for the events class of objects.

events.PSTH=[];
events.PSTHBins=[];
events.StartEvent=[];
events.EndEvent=[];
events.EventProbability=[];
events.EventDurations=[];
events.Thresholds=[];
events.EventMeanSpikeCounts=[];
events.EventFanoFactor=[];
events.EventMeanSpikeTimes=[];
events.MaxPSTHValue=[];
events.MeanPSTHValue=[];
events.EventStartingPoints=[];
events.EventIntervalDurations=[];
events.EventSubStim=[];
events.EventSubStimIndex=[];
events.SubStimReps=[];

%%% This value is the threshold of repetitons that need to have spikes,
%%% this is also a way of controlling which events are to be considered
%%% based on their reliability.
events.PercentageReps = Args.PercentageReps;
events.SlideDuration = Args.SlideDuration;
events.FrameDuration = stimInfo.data.catInfo.frame_duration;
events.RepetitionDuration = stimInfo.data.catInfo.repetition_duration;
events.NumFrames = stimInfo.data.catInfo.num_frames;
events.NumRepetitions = stimInfo.data.catInfo.num_repetitions;
events.ThresholdType = Args.ThresholdType;
events.Upper = Args.Upper;
events.Lower = Args.Lower;

if isnumeric(Args.Binsize)
    events.Binsize = Args.Binsize;
else
    events.Binsize = round(events.FrameDuration);
end
convraster = obj.data.raster';
Spiketrain = obj.data.adjSpiketrain;
FramePoints = obj.data.adjFramePoints;
%%% Creation of the raster %%%%%%
for r = 1:events.NumRepetitions
    %%%% Frames for the repetition %%%%
    Frames = FramePoints(r*events.NumFrames-events.NumFrames+1:(r*events.NumFrames)+1);
    %%%% Spiketimes for the given repetition %%%%%
    RepSpikes = Spiketrain(find(Spiketrain > Frames(1) & Spiketrain <= Frames(end)));    
    if r > 1    
        raster{r} = RepSpikes-TimeAdjust;
        TimeAdjust = Frames(end); % This adjusts the spiketimes for each repetition of the stimulus
    else
        raster{r} = RepSpikes-Frames(1);
        TimeAdjust = Frames(end); % This adjusts the spiketimes for each repetition of the stimulus
    end    
end

% Obtain the initial event boundries from the raster and psth via a created threshold 
[events] = EventsLimits(raster,events,convraster);
if ~isempty(events.StartEvent)       
    % Check the event boundaries and make sure the valleys are below the
    % mean firing rate and the boundaries are wide enough.
    [events] = EventsLimitsConfirm(events);
    if ~isempty(events.StartEvent)        
        % Obtain the spikes within the event boundries and calculate the
        % event intervals the last input variable is the percentage of
        % repetitions with a single spike. This will toss many events if
        % the percentage level is too high. 
        [events] = EventsSpikes(raster,events);
        if ~isempty(events.StartEvent)
            % Gives the relative weight for each event based on the area under the psth for the event
            % It will also depend on the event duration.
            %[events] = EventsWeights(events);
            % Calculate the durations for each event and the amount of time the cell
            % fires within an event, event probability.
            [events] = EventsDurations(events);
            %% Get the mean spike time for the first bin with a value of 75
            %% percent of the maximum firing rate bin.
            [events] = EventsStartingPoints(raster,events);
        else
            events.StartEvent=[];
            events.EndEvent=[];
            events.EventProbability=[];
            events.EventDurations=[];
            events.EventMeanSpikeCounts=[];
            events.EventFanoFactor=[];
            events.EventMeanSpikeTimes=[];
            events.MaxPSTHValue=[];
            events.MeanPSTHValue=[];
            events.EventStartingPoints=[];
            events.EventIntervalDurations=[];
            events.EventSubStim=[];
            events.EventSubStimIndex=[];
            events.SubStimReps=[];
        end
    else
        events.StartEvent=[];
        events.EndEvent=[];
        events.EventProbability=[];
        events.EventDurations=[];
        events.EventMeanSpikeCounts=[];
        events.EventFanoFactor=[];
        events.EventMeanSpikeTimes=[];
        events.MaxPSTHValue=[];
        events.MeanPSTHValue=[];
        events.EventStartingPoints=[];
        events.EventIntervalDurations=[];
        events.EventSubStim=[];
        events.EventSubStimIndex=[];
        events.SubStimReps=[];
    end
end 