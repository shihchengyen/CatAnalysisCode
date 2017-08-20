function displaysubstimevents(cell_num,binsize)

load jonsdisplayevents 

for x = 1:length(displayevents)
    
    cell(x) = displayevents(x).cell_num;
    
end

cell_ind = find(cell == cell_num);
rez = [3 5];
rez_ind = find(rez == binsize);

psth = displayevents(cell_ind).normalpsths{rez_ind};
threshold = displayevents(cell_ind).thresholds(rez_ind);
threshold2 = displayevents(cell_ind).thresholds2(rez_ind);
raster = displayevents(cell_ind).raster;
rep_duration = displayevents(cell_ind).rep_duration;
frame_duration = displayevents(cell_ind).frame_duration;
samplingrate = displayevents(cell_ind).samplingrate;
num_repetition = displayevents(cell_ind).num_repetition;
events = displayevents(cell_ind).mean_event_times{rez_ind};
pos_lim = displayevents(cell_ind).pos_limit{rez_ind};
neg_lim = displayevents(cell_ind).neg_limit{rez_ind};
psthline = displayevents(cell_ind).bin{rez_ind};
ss_ind = displayevents(cell_ind).ss_event_ind{rez_ind};

frames = (0:frame_duration:rep_duration);

figure

%%%%%%%PLOT THE RASTER%%%%%%%%%%%%%

for i_r = 1:num_repetition
    
    if isempty(raster{i_r})
        raster{i_r} = NaN;
    end
    hold on
    plot(raster{i_r},i_r+threshold+10,'.k');
    
end

%%%%%%%%%%%%%% PLOT THE FRAME MARKERS %%%%%%%%%%%%%%%%%%%
    
        
    for f = 1:length(frames)
    
        hold on
        line([frames(f) frames(f)],[threshold+10 threshold+7],'Color','g',...
                                                     'LineWidth',3)
    end
    

    
 %%%%%%%%%%%%%% PLOT THE EVENTS %%%%%%%%%%%%%%%%%%%%%%
    
 for f = 1:length(events)
    
        hold on
        line([events(f) events(f)],[threshold+7 threshold+4],'Color','r',...
                                                     'LineWidth',3)
 end
        
 
 %%%%%%%%%%%%% PLOT THE EVENT LIMITS %%%%%%%%%%%%%%%%

 for f = 1:length(pos_lim)
    
        hold on
        line([pos_lim(f)*binsize pos_lim(f)*binsize],[threshold+7 threshold+4],'Color','m',...
                                                     'LineWidth',3)
        hold on
        line([(neg_lim(f)-1)*binsize (neg_lim(f)-1)*binsize],[threshold+7 threshold+4],'Color','c',...
                                                     'LineWidth',3)
                                                
 end

  
%%%%%%%%%%%%% PLOT THE PSTH AND THRESHOLD %%%%%%%%%%%%%%%%

hold on
stairs(psthline,psth,'b')
plot([0 rep_duration],[threshold threshold],'--r')
plot([0 rep_duration],[threshold2 threshold2],'--c')
title(['Cell Number ',num2str(cell_num) ' Binsize ',num2str(binsize)])
axis tight

if isempty(ss_ind)
    fprintf('NO SUB STIMULUS EVENTS')
    return;
end

for e = 1:length(ss_ind);
    
    key = input('RETURN - Next Trial; Q - Quit: ','s');
    n = str2num(key);
	
    if strcmp(key,'q')
        return;
        
    else
        ylim([0 threshold+5+num_repetition])
        xlim([events(ss_ind(e))-3*frame_duration events(ss_ind(e))+3*frame_duration])
    end
        
end

fprintf('END OF EVENTS')
 
return

while 1
	% get keyboard input to see what to do next
	key = input('RETURN - Next; p - Previous; q - Quit: ','s');
	if strcmp(key,'q')
		break;
	elseif strcmp(key,'p')
		xmin = xmin - steps;
		xmax = xmax - steps;
	else
		xmin = xmin + steps;
		xmax = xmax + steps;
	end
	axis([xmin xmax ax1(3) ax1(4)])
end
    
