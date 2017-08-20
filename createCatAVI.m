function createCatAVI(filename,frameNums,repNums)
%createCatAVI(filename,clusterlist,frameNums,repNums)
%
%This function will create two files, a sound file and a movie file.
%These can be merged in Adobe Premiere to create an AVI.
%
%inputs:
%frameNums - vector of start and end frame
%repNums - vector of start and end repetitions
%
%session dir must be pwd.
%the output files will be named as the clusterlist file.

%%%% Load the stimulus Information
si = stiminfo('auto','save');
clusterlist = textread('cell_list.txt','%s');

%%%% Get the start and end times %%%%
startTime = si.data.catInfo.frame_duration*(frameNums(1)-1);  %start time of clip
finishTime = si.data.catInfo.frame_duration*(frameNums(2));  %finish time of clip
x_off = (si.data.iniInfo.ScreenWidth - si.data.iniInfo.frame_cols)/2;
y_off = (si.data.iniInfo.ScreenHeight - si.data.iniInfo.frame_rows)/2;
RF_Points = round(mean(si.data.iniInfo.ReceptiveField.Points,2));
RF_Points(1:2:end) = RF_Points(1:2:end)-x_off;
RF_Points(2:2:end) = RF_Points(1:2:end)-y_off;
Center_Points = [round(mean(si.data.iniInfo.ReceptiveField.CenterX))-x_off round(mean(si.data.iniInfo.ReceptiveField.CenterY))-y_off];

frameNums = (si.data.iniInfo.start_frame + frameNums)-1;

switch si.data.iniInfo.stimulus_root
    case('Biglebowski')
        movienum = 1;
    case('Everest')
        movienum = 2;
    case('Animals')
        movienum = 3;
    case('Cats')
        movienum = 4;
end

%%%% Create the image array and then write out an avi file %%%%%
count=1;
for ii = frameNums(1):frameNums(2)   
    image = readpglframe(ii,movienum);
%     for c = 1:4
%         points = RF_Points(1:2);
%         image(points(1)-5:points(1)+5,points(2)-5:points(2)+5) = 255;
%         RF_Points(1:2)=[];
%     end   
    image(Center_Points(2)-2:Center_Points(2)+2,Center_Points(1)-2:Center_Points(1)+2) = 255;
    frames(:,:,1,count) = image;
    count = count+1;
end
map = (0:(1/255):1)';
map = [map map map]; %%% Colormap for immovie function 
mov = immovie(frames,map); %% Converts the frames matrix into a matlab movie format
fps = si.data.catInfo.video_refresh/si.data.iniInfo.refreshes_per_frame;
movie2avi(mov,filename,'compression','None','quality',100,'fps',fps)
fprintf('Finished Movie')

%%%% Get the cell info %%%%%%
asObj = processSession(adjspikes,'Cells',clusterlist,'save','redo'); %cummulative adjspikes object
numClusters = size(clusterlist,1);
tt = [0:.0001:.1]; %each spike chirp will last 100 milliseconds
for ss = 1:length(repNums)
    for ii = 1:numClusters
        freq = ii*250; %freq of chirp at 200 HZ intervals
        chirp = cos(2*pi*freq*tt);
        spiketimes = round(asObj.data.raster(ii).spikes(repNums(ss),find(asObj.data.raster(ii).spikes(repNums(ss),:) >= startTime  &  asObj.data.raster(ii).spikes(repNums(ss),:) < finishTime))-startTime);   
        wave_file = zeros(1,length((startTime:1:finishTime)-startTime)*1000);
        for a = 1:length(spiketimes)
            wave_file(spiketimes(a)*1000:(spiketimes(a)*1000)+999) = chirp(1:1000);
        end
        clusters_file(ii,:) = wave_file;
        fprintf(['Finished with Cluster: ',num2str(ii) '\n'])
    end
    wave_file = sum(clusters_file,1); clear clusters_file
    wavwrite(wave_file,1000000,[filename num2str(ss)])      
    fprintf(['Finished with Repetition: ',num2str(ss) '\n'])
end