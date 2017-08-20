% function displaysubstimulusevents(binsize)

% cd d:/data

binsize = 3;

load jonsevents

rez = [3 5];
rez_ind = find(rez == binsize);
ss_ind=[];
rank=[];
for c = 1:length(events)
    
    ss_events = events(c).ss_event_times{rez_ind};
    
    num_intervals = events(c).number_of_intervals(rez_ind);
    
    if ~isempty(ss_events)
        
        ss_ind = [ss_ind c];
        
        rank = [rank num_intervals];
        
    end
    
end

figure
all_ss_events=[];
for c = 1:length(ss_ind)
    
    rank_ind = find(rank == min(rank));
    rank_ind = rank_ind(1);
    
    rank(rank_ind) = [];
    
    ss_events = events(ss_ind(rank_ind)).ss_event_times{rez_ind};
    
    ss_ind(rank_ind) = [];
    
    all_ss_events = [all_ss_events ss_events];
    
    for x = 1:length(ss_events)
        line([c-.5 c+.5],[ss_events(x) ss_events(x)],'Color','k','LineWidth',2)
        hold on
    end
    
end
set(gca,'XTick',1:1:10)
xlabel('Cell Number')
ylabel('Event Interval (ms)')
ylim([0 35])
xlim([0 11])

figure

[N] = histc(all_ss_events,0:5:35);
barh(0:5:35,N,'b')
ylim([0 35])
xlabel('# Occurrences')
ylabel('Event Interval (ms)')

