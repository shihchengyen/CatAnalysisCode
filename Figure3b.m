if ispresent('fig3data.mat','file')
	load fig3data
else
	% for figure a & b
	load jonsevents

	% set a few constants
	binsize = 3;
	bin = [3 5];
	rez = find(bin == binsize);
	num_sub_bins = 10;
	sub_bins = [10 6 3 2];
	nsb = find(sub_bins == num_sub_bins);
	
	% get number of cells
	ncells = length(events);
	
	% initialize variables
	cevents = []; % all event intervals
	nevents = []; % number of event intervals for each cell
	nssevents = []; % number of sub-stimulus intervals for each cell
	eventcells = {}; % sub-stimulus intervals for each cell
	eventcellid = [];
	
	% for figure c
	load jonsdisplayevents
		
	% initialize index for eventcells
	ec = 1;
	for i = 1:ncells
		% check to see if there are sub-stimulus events in this cell
		if ~isempty(events(i).ss_event_times{rez})
			% add all events from this cell to total
			cevents = [cevents events(i).all_event_intervals{1}];
			nevents = [nevents; length(events(i).all_event_intervals{rez})];
			nssevents = [nssevents; length(events(i).ss_event_times{rez})];
			eventcellid = [eventcellid; events(i).cell];
			eventcells{ec} = events(i).ss_event_times{rez};
			ec = ec + 1;
		end
		
		raster = displayevents(i).raster;
		
		% find sub-stimulus events for cell e
		ss_event_ind = displayevents(i).ss_event_ind{rez};
				
		if isempty(ss_event_ind);
			% if there were no sub-stimulus events
			all_ss_events(i).ss_iei_mat = [];
		else
			% find limits of sub-stimulus events for cell e
			pos_lim = displayevents(i).first_pos_limit{rez};
			pos_lim1 = pos_lim(ss_event_ind);
			pos_lim2 = pos_lim(ss_event_ind+1);
			
			neg_lim = displayevents(i).first_neg_limit{rez};
			neg_lim1 = neg_lim(ss_event_ind);
			neg_lim2 = neg_lim(ss_event_ind+1);
			
			% get number of sub-stimulus events
			nss = length(pos_lim1);
			ss_events = zeros(1,nss);
			
			% loop over sub-stimulus events in cell e
			for x = 1:nss
				% get number of repetitions - can't move this outside for 
				% loop since some cells have different number of reps
				reps = length(raster);
				ss_event_mat = zeros(reps,2);
				
				% loop over repetitions
				for r = 1:reps
					% for repetition r check to see if there were spikes 
					% within event limits
					iei_ind_1 = find((raster{r}>=(neg_lim1(x)*binsize)-binsize)... 
									 & (raster{r}<pos_lim1(x)*binsize));
					if isempty(iei_ind_1)
						iei_ind_1 = 0;
					else
						iei_ind_1 = 1;
					end
										
					iei_ind_2 = find((raster{r}>=(neg_lim2(x)*binsize)-binsize)... 
									 & (raster{r}<pos_lim2(x)*binsize));
					if isempty(iei_ind_2)
						iei_ind_2 = 0;
					else
						iei_ind_2 = 1;
					end
					
					ss_event_mat(r,1) = iei_ind_1;
					ss_event_mat(r,2) = iei_ind_2;
				end % end loop over repetitions
				
				% get percentage of repetitions that have spikes in both 
				% sub-stimulus events
				% ss_event_mat = ss_event_mat(:,1).*ss_event_mat(:,2);
				% ss_events(x) = sum(ss_event_mat)/length(ss_event_mat)*100;
				ss_events(x) = (ss_event_mat(:,1)' * ss_event_mat(:,2)) / reps * 100;			
			end % end loop over sub-stimulus events in cell e
			
			all_ss_events(i).cell_num = displayevents(i).cell_num;
			all_ss_events(i).ss_iei_mat = ss_events;
			% ss_events=[];
		end %if empty ss events
	end %for num cells
	
	nev = [nevents nssevents eventcellid];
	[nevsort,nevsi] = sortrows(nev);

	cmax = max(cevents);
	bins = 0:5:(cmax+5);
	ncevents = histcie(cevents,bins);
	bins2 = 0:1:cmax;
	ncevents2 = histcie(cevents,bins2);
	
	% get overlap percentages for all cells
    ss_dist=[];
    for x = 1:ncells
        ss_dist = [ss_dist all_ss_events(x).ss_iei_mat];
    end
    
    % get cummulative sum for number of sub-stimulus events
    % add a zero in front to make it easier to use
    cnss = [0; tril(ones(10)) * nev(:,2)];    
    
    ss_overlap = [];
	bstart = 1;
	for i = 1:10
		% get number of sub-stimulus intervals
		nssint = nevsort(i,2);
		% get order of cell
		ci = nevsi(i);
		% get start and end indices
		nstart = cnss(ci) + 1;
		nend = cnss(ci+1);
		bend = bstart + nssint - 1;
		ssoverlap = sortrows([eventcells{nevsi(i)}' ss_dist(nstart:nend)']);
		ss_overlap = [ss_overlap; (bstart:bend)' ssoverlap(:,2)];
		nstart = nend + 1;
		bstart = bend + 2;
	end

    % for figures d, e, and f
	load jonscellsdata

	% for figure d	
	all_fano_scene_rev = [];
	% for figure e
	all_tres_scene_rev = [];
	% for figure f
	all_event_scene_rev = [];	
	
	for x = 1:ncells
		sc_num = jonscells(x).scene_changes;
		if sc_num(1) ~= 1
			sc_num = [1 sc_num];
		end
		
		fano = jonscells(x).fanofactor;
		msc = jonscells(x).meanspikecount;
		% find frames that have msc>=1 
		msc1 = find(msc >= 1);
		% find frames that have msc>=1 and FF<1
		fanos_below = intersect(find(fano < 1),msc1);
		fano_scene_rev = [];
		
		frames = jonscells(x).frame_times;
		frames = frames(1:end-1);
		tres = jonscells(x).tresscores{nsb};
		
		event_times = events(x).all_event_times{rez};
		ss_events = event_times(events(x).ss_event_ind{rez});
	
		for f = 1:length(fanos_below)
			sc_frame = fliplr(fanos_below(f)-1:-1:fanos_below(f)-5);
			sc_rev = [0 0 0 0 0];
			for s = 1:length(sc_num);
				sc_rev_ind = [sc_frame == sc_num(s)];
				sc_rev = sc_rev+sc_rev_ind;
			end
			
			fano_scene_rev = [fano_scene_rev;sc_rev];
		end
			
		all_fano_scene_rev = [all_fano_scene_rev;fano_scene_rev];

		% find frames that have msc>=1 and TRES > 99
		tres_above = intersect(find(tres > 99),msc1);
		tres_scene_rev = [];
		
		for f = 1:length(tres_above)
			sc_frame = fliplr(tres_above(f)-1:-1:tres_above(f)-5);
			sc_rev = [0 0 0 0 0];
			for s = 1:length(sc_num);
				sc_rev_ind = [sc_frame == sc_num(s)];
				sc_rev = sc_rev+sc_rev_ind;
			end
			
			tres_scene_rev = [tres_scene_rev;sc_rev];
		end
			
		all_tres_scene_rev = [all_tres_scene_rev;tres_scene_rev];

		event_scene_rev = [];
		if ~isempty(ss_events)
			ss_event_frames=[];
			for f = 1:length(frames)-1
				ss_ind = find(ss_events>frames(f) & ss_events<frames(f+1));
				if ~isempty(ss_ind)
					num_ss = length(ss_ind);
					ff = zeros(1,num_ss);
					ff = ff+f;
					ss_event_frames = [ss_event_frames ff];
				end
			end
				   
			for f = 1:length(ss_event_frames)
				sc_frame = fliplr(ss_event_frames(f)-1:-1:ss_event_frames(f)-5);
				sc_rev = [0 0 0 0 0];
				for s = 1:length(sc_num);
					sc_rev_ind = [sc_frame == sc_num(s)];
					sc_rev = sc_rev+sc_rev_ind;
				end
				event_scene_rev = [event_scene_rev;sc_rev];
			end
		end %if empty
		
		all_event_scene_rev = [all_event_scene_rev;event_scene_rev];
	end % end for loop over cells
end

figure
a1 = axes('Position',[0.13 0.11 0.2128 0.683]);
% h1 = subplot(1,3,1);
hb = bar(nevsort(:,1));
set(hb,'FaceColor',[0.5 0.5 0.5])
bv = 0.8;
vert = get(hb,'vertices');
vert(vert==0) = bv;
set(hb,'vertices',vert);
set(gca,'YScale','log','YTick',[1 2 5 10 20 55],'YMinorTick','off');
hold on
% get the handle so we can change the color
hb = bar(nevsort(:,2),'k');
vert = get(hb,'vertices');
vert(vert==0) = bv;
set(hb,'vertices',vert);
axis([0 11 0.8 55])
% axis square
hold off
xlabel('Cell number')
% ylabel('# of intervals')
ha = axes('Position',[0.13 0.8 0.2128 0.125]);
set(ha,'XAxisLocation','top','YTick',[],'TickDir','out','Box','on')
for i = 1:10
	% get x values
	xv = eventcells{nevsi(i)};
	% set y values
	y1 = i - 1;
	y2 = y1 + 0.5;
	line([xv;xv],[y1;y2],'Color','k');
end
axis([0 35 -0.5 10])

% axes('Position',[0.4111 0.11 0.2128 0.815])
subplot(1,3,2)
% get histogram of event intervals
hb = bar(bins,ncevents,0.8,'histc');
set(hb,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0]);
hold on
hb = bar(bins(1:7),ncevents(1:7),0.8,'histc');
set(hb,'FaceColor','k');
set(gca,'XTick',[0 35 70 105 140],'TickDir','out')
xlim([0 150])
sync1 = 1000/85.1;
sync2 = 2 * sync1;
sync3 = 3 * sync1;
xlabel('Inter-event interval (ms)')
% ylabel('# occurrences')
% draw inset with distribution of sub-stimulus event intervals
ha = axes('Position',[0.493 0.755 0.122 0.165]);
hb2 = bar(bins2,ncevents2,0.8,'histc');
set(hb2,'FaceColor','k','EdgeColor',[0 0 0]);
line([sync1 sync1],[0 6],'Color',[0.5 0.5 0.5],'LineStyle',':')
line([sync2 sync2],[0 6],'Color',[0.5 0.5 0.5],'LineStyle',':')
line([sync3 sync3],[0 6],'Color',[0.5 0.5 0.5],'LineStyle',':')
axis([0 35 0 6])
set(ha,'XTick',[0 12 23 35],'TickDir','out','TickLength',[0.03 0.025],'YTick',[2 4])
hold off

% axes('Position',[0.6922 0.11 0.2128 0.815])
subplot(1,3,3)
% hold on
% hb = bar(bstart:bend,ssoverlap(:,2));
hb = bar(ss_overlap(:,1),ss_overlap(:,2));
set(hb,'FaceColor',[0.5 0.5 0.5]);
set(gca,'XTick',[1.5 4 7 10 12 14.5 20 26.5 30 34.5],'XTickLabel',1:10,...
	'box','on');
xlim([0 38])
xlabel('Cell number')
% ylabel('Repetitions(%)')
hold off

% ssbins = 0:5:100;
% sshc = histcie(ss_dist,ssbins);
% hb = bar(ssbins,sshc,0.8,'histc');
% set(hb,'FaceColor',[0.5 0.5 0.5]);
% xlim([-5 105])
% xlabel('Repetitions (%)')
% ylabel('# of occurrences')

figure
subplot(1,3,1)
hb = bar(sum(all_fano_scene_rev)/size(all_fano_scene_rev,1)*100,'k');
set(hb,'FaceColor',[0.5 0.5 0.5]);
set(gca,'XTickLabel',5:-1:1)    
axis square
xlabel('Frames before FF < 1')
% ylabel('# of scene changes (%)')
ylim([0 8])

subplot(1,3,2)
hb = bar(sum(all_tres_scene_rev)/size(all_tres_scene_rev,1)*100,'k');
set(hb,'FaceColor',[0.5 0.5 0.5]);
set(gca,'XTickLabel',5:-1:1)    
axis square
xlabel('Frames before TRESScore=100')
% ylabel('# of scene changes (%)')

subplot(1,3,3)
hb = bar(sum(all_event_scene_rev)/size(all_event_scene_rev,1)*100,'k');
set(hb,'FaceColor',[0.5 0.5 0.5]);
set(gca,'XTickLabel',5:-1:1)    
axis square
xlabel('Frames before sub-stimulus interval')
% ylabel('# of scene changes (%)')
