function [mseq_frames,spike_counts] = getmseqframes(stimInfo,adjspikes,num_frames_back,Args)
% function to get the mseq frames for the stcov class
mseq_size = num2str(stimInfo.data.iniInfo.m_seq_size(1));
if strcmp(num2str(mseq_size),'16')
    numRows = 16;
    numCols = 16;
    order = 16;
    load mseq16grid
elseif strcmp(num2str(mseq_size),'64')
    numRows = 64;
    numCols = 64;
    order = 16;
    load mseq64grid % Load the mseq grid file.
end
%%%%%%%%%% Get the Spikes %%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(stimInfo.data.iniInfo,'last_frames')
    f_numbers = stimInfo.data.iniInfo.first_frames:stimInfo.data.iniInfo.last_frames;
else
    f_numbers = stimInfo.data.iniInfo.first_frames:stimInfo.data.iniInfo.frames_displayed;
    f_numbers = f_numbers(1:stimInfo.data.iniInfo.frames_displayed);
end
spikes = adjspikes.data.adjSpiketrain;
frames_vector = adjspikes.data.adjFramePoints;
if frames_vector(1) < 1 %% This is due to the conversion from 1 datapoint to ms
    frames_vector(1) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spike_counts = histcie(spikes,frames_vector,'DropLast')';
frame_numbers = f_numbers(find(spike_counts > Args.NumSpikes));
spike_counts = spike_counts(find(spike_counts > Args.NumSpikes));
frame_numbers = frame_numbers-num_frames_back;

clear stimInfo adjspikes
mseq_frames=zeros(numRows,numCols,length(frame_numbers));
for f = 1:length(frame_numbers)
    mseq_frames(:,:,f) = MSeqframe(mseqgrid,numRows,numCols,order,frame_numbers(f))/255; 
    %fprintf(['Frames Left: ',num2str(length(frame_numbers)-f) '\n'])
end