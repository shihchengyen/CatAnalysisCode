function displayeventspsth(cell_num,bin)

load jpdisplayevents %jpdisplaynormalmaxrand0

for x = 1:length(displayevents)
    
    cell(x) = displayevents(x).cell_num;
    
end

cell_ind = find(cell == cell_num);
rez = [3,5:5:15];
rez_ind = find(rez == bin);

psth = displayevents(cell_ind).normalpsths{rez_ind};
threshold = displayevents(cell_ind).thresholds(rez_ind);
raster = displayevents(cell_ind).raster;
rep_duration = displayevents(cell_ind).rep_duration;
frame_duration = displayevents(cell_ind).frame_duration;
samplingrate = displayevents(cell_ind).samplingrate;
num_repetition = displayevents(cell_ind).num_repetition;
events = displayevents(cell_ind).mean_event_times{rez_ind};

frames = fix(frame_duration/(samplingrate/1000));
frames = (1:frames:rep_duration);

figure

%%%%%%%PLOT THE RASTER%%%%%%%%%%%%%

for i_r = 1:num_repetition
    
    if isempty(raster{i_r})
        raster{i_r} = NaN;
    end
    hold on
    plot(raster{i_r},i_r+threshold+5,'.k');
    
end

%%%%%%%%%%%%%% PLOT THE FRAME MARKERS %%%%%%%%%%%%%%%%%%%
    
        
    for f = 1:length(frames)
    
        hold on
        line([frames(f) frames(f)],[num_repetition+1 num_repetition+2],'Color','g',...
                                                     'LineWidth',4)
    end
    

    
 %%%%%%%%%%%%%% PLOT THE EVENTS %%%%%%%%%%%%%%%%%%%%%%
    
 for f = 1:length(events)
    
        hold on
        line([events(f) events(f)],[num_repetition+3 num_repetition+4],'Color','r',...
                                                     'LineWidth',4)
 end
        
    

%%%%%%%%%%%%% GET THE PSTH INFO AND THRESHOLD %%%%%%%%%%%%%%%%

frame = frame_duration/(samplingrate/1000);
psthline = 0:bin:rep_duration;

if length(psthline) > length(psth)
    psthline = psthline(1:length(psth));
else
    psth = psth(1:length(psthline));
end

%%%%%%%%%%%%% PLOT THE PSTH AND THRESHOLD %%%%%%%%%%%%%%%%

hold on
stairs(psthline,psth,'b')
plot([0 rep_duration],[threshold threshold],'--r')
title(['Cell Number ',num2str(cell_num) ' Binsize ',num2str(bin)])
axis tight

iei_ind = find(diff(events)<frame);
iei = iei_ind+1;

if isempty(iei_ind)
    fprintf('NO SUB STIMULUS EVENTS')
    return;
end

for e = 1:length(iei_ind);
    
    key = input('RETURN - Next Trial; Q - Quit: ','s');
    n = str2num(key);
	
    if strcmp(key,'q')
        return;
        
    else
        ylim([threshold-30 threshold+5+num_repetition])
        xlim([events(iei(e))-2*frame events(iei(e))+2*frame])
    end
        
end

fprintf('END OF EVENTS')
    
    