function [y,t] = adjPSTH(obj,stimInfo,width,slide)
% adjPSTH(obj,stimInfo,binsize,slide_duration)
% input obj is the adjspikes and stiminfo objects to create and then
% calculate using the slidinghist program the post stimulus time histogram (PSTH) 

frame_duration = stimInfo.data.catInfo.frame_duration;
repetition_duration = stimInfo.data.catInfo.repetition_duration;
num_repetitions = stimInfo.data.catInfo.num_repetitions;
% Check for Spontaneous Activity, only found in cells from P1 and on.
if isfield(stimInfo.data.iniInfo, 'spontaneous_activity_duration_in_milliseconds')
    repetition_duration = repetition_duration+stimInfo.data.iniInfo.spontaneous_activity_duration_in_milliseconds;
end

if ~isnumeric(width)
    width = frame_duration;
end
if ~isnumeric(slide)
    slide = width;
end

[y,t] = slidingHist(obj.data.raster',slide,width,repetition_duration);

%change to spikes per second
y=(y/num_repetitions)/(width/1000);