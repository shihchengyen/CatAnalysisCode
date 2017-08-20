% function figure(cell_num,binsize,num_sub_bins)

% define some constants
cell_num = 6;
binsize = 3;
num_sub_bins = 10;
% show frame makers
showFMarkers = 1;

% Position values are [left bottom width height]
pos1 = [.075 .45 .875 .525];
pos2 = [.075 .15 .875 .29];
pos3 = [.075 .1 .875 .04];
pos4 = [.075 .05 .875 .04];

% Loads the 3 files below and calls no other functions, must include a 
% cell number, binsize for the psth and events, and the number of 
% sub_bins for the tres score.

load jonsdisplayevents
load jonscellsdata
load jbcells

cell_ind = find(cell_numbers == cell_num);
rez = [3 5];
rez_ind = find(rez == binsize);

psth = displayevents(cell_ind).normalpsths{rez_ind};
threshold = displayevents(cell_ind).thresholds(rez_ind);
threshold2 = displayevents(cell_ind).thresholds2(rez_ind);
events = displayevents(cell_ind).mean_event_times{rez_ind};
pos_lim = displayevents(cell_ind).pos_limit{rez_ind};
neg_lim = displayevents(cell_ind).neg_limit{rez_ind};
psthline = displayevents(cell_ind).bin{rez_ind};
ss_event_ind = displayevents(cell_ind).ss_event_ind{rez_ind};
raster = jonscells(cell_ind).raster;
sc_num = jonscells(cell_ind).scene_changes;
frames = jonscells(cell_ind).frame_times;
rep_duration = jonscells(cell_ind).rep_duration;
frame_duration = jonscells(cell_ind).frame_duration;
num_repetition = jonscells(cell_ind).num_repetition;

figure
% turn double buffering on so updates will be faster
set(gcf,'DoubleBuffer','on');
% set(gcf,'Color',[1 1 1])

%%%%%%%PLOT THE RASTER%%%%%%%%%%%%%
%subplot('Position',[left bottom width height])
subplot('Position',pos1)
set(gca,'XAxisLocation','top') 
for r = 1:num_repetition
    if isempty(raster{r})
        raster{r} = NaN;
    end
    hold on
%     line([raster{r} raster{r}],[r+.5 r-.5],'Color','b','Linewidth',2)
    plot(raster{r},r,'.k')
end
fram = frames(2:end);

%%%%%%%%%PLOT THE SCENE CHANGES %%%%%%%%%%%%%

for sc = 1:length(sc_num)
   plot(fram(sc_num(sc)),num_repetition+1.5,'k^','LineWidth',2)
end

%%%%%%%%%%%%%% PLOT THE EVENTS %%%%%%%%%%%%%%%%%%%%%%
    
for f = 1:length(events)
   line([events(f) events(f)],[num_repetition+3 num_repetition+4.5],'Color','k',...
                                                     'LineWidth',1.5)
end
 
%%%%%%%%%%%%%% PLOT THE FRAME MARKERS %%%%%%%%%%%%%%%%%%%
if showFMarkers
	for f = 1:length(frames)
	   line([frames(f) frames(f)],[num_repetition+5 num_repetition+6.5],'Color','k',...
														'LineWidth',1.5)
	end
end

%  %%%%%%%%%%%%%%%%%% PLOT THE EVENT BOUNDARIES %%%%%%%%%%%%%%%
%  
%  
%  for f = 1:length(pos_lim)
%     
%        line([pos_lim(f)*binsize pos_lim(f)*binsize],[num_repetition+4 num_repetition+5.5],'Color','c',...
%                                                      'LineWidth',4)
%                                                  
%        line([neg_lim(f)*binsize neg_lim(f)*binsize],[num_repetition+4 num_repetition+5.5],'Color','c',...
%                                                      'LineWidth',4)
%  end
        
% ylabel('Repetitons')
% axis tight
if showFMarkers
	liney = num_repetition + 6.5;
	ylim([0 num_repetition+6.5])
else
	liney = num_repetition + 5.5;
	ylim([0 num_repetition+5.5])
end
line([0 rep_duration],[liney liney],'Color','k') 
axis ij
set(gca,'box','off')
set(gca,'Tickdir','out')
set(gca,'XTickmode','manual'); set(gca,'XTick',[]);
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'YTick',[0:10:num_repetition]);


%%%%%%%%%%%%%% SUBPLOT 2 %%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%% PLOT THE PSTH AND THRESHOLD %%%%%%%%%%%%%%%%

subplot('Position',pos2)
stairs(psthline,psth,'k')
% ylabel('FR')
hold on
plot([0 rep_duration],[threshold threshold],'--k')
% hold on
plot([0 rep_duration],[threshold2 threshold2],'--k')
axis tight
set(gca,'box','off')
set(gca,'Tickdir','out')
% set(gca,'XTickmode','manual'); 
set(gca,'XTick',[]);
set(gca,'YTick',[50 100 150]);
% set(gca,'XAxisLocation','top')
set(gca,'XColor',[1 1 1])
set(gca,'YColor',[0 0 0])


%%%%%%%%%%%%%% SUBPLOT 3 %%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%% GET THE FANO FACTOR INFO %%%%%%%%%%%%%%%%

fano = jonscells(cell_ind).fanofactor;
msc = jonscells(cell_ind).meanspikecount;

msc_ind = find(msc<1);
fano(msc_ind) = NaN;

% fano_ind = find(fano>=5);
% fano(fano_ind) = 5;

%%%%%%%%%%%%% PLOT THE FANO FACTOR %%%%%%%%%%%%%%%%

subplot('Position',pos3)
hb = bar(frames(1:end-1)+(frame_duration/2),fano,'k');
% vert=get(hb,'vertices');
% vert(vert==0) = floor(min(fano)*100)/100;
% set(hb,'vertices',vert);
% set(gca,'yscale','log','ytick',[(floor(min(fano)*10)/10):.1:1 2:1:5]);
% set(gca,'yscale','log','ytick',[.3:.2:.9]);
ylim([0.4 1])
set(gca,'YTick',[0.6 1]);

hold on
line([0 rep_duration],[1 1],'Color','k')
% ylabel('FF')
% axis tight
set(gca,'box','off')
set(gca,'Tickdir','out')
% set(gca,'XTickmode','manual'); 
set(gca,'XTick',[]);
set(gca,'XColor',[1 1 1])
set(gca,'YColor',[0 0 0])


%%%%%%%%%%%%%% SUBPLOT 4 %%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%% PLOT THE TRES SCORE %%%%%%%%%%%%%%%%%%
 
sub_bins = [10 6 3 2];

nsb = find(sub_bins == num_sub_bins);

tres = jonscells(cell_ind).tresscores{nsb};

tres(msc_ind) = NaN;

subplot('Position',pos4)
bar(frames(1:end-1)+(frame_duration/2),tres,'k')
% ylabel('TRES')
% axis tight
set(gca,'box','off')
set(gca,'Tickdir','out')
% xlabel('Time (msec)')
ylim([50 100])
set(gca,'YTick',[50 100])
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])

% print -dtiff sciencefig1b.tiff



