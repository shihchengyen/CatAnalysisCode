function displayjonsevents()
% This function will load the event data from the desired structure and plot a summary figure. psthtype variable number of 1 is normal psth data, otherwise as sliding psth.

% cd d:/data

load jonsevents

rez = [3 5];

for f = 1:length(rez)
    
    event_intervals = [];
    ss_intervals = [];
    event_dist = zeros(length(events),2);
    
        for ff = 1:length(events)
            
            event_intervals = [event_intervals events(ff).all_event_intervals{f}];
            
            ss_intervals = [ss_intervals events(ff).ss_event_times{f}];
            
            if ~isempty(events(ff).ss_event_times{f})
                cell_below(ff) = 1;
                event_dist(ff,1) = 1;
                event_dist(ff,2) = length(events(ff).ss_event_times{f});
            else
                cell_below(ff) = 0;
            end
                        
        end
                      
        all_event_intervals{f} = event_intervals;
        num_intervals(f) = length(event_intervals);
        num_ss_intervals(f) = length(ss_intervals);
        num_ss_cells(f) = sum(cell_below);
        ss_event_dist{f} = event_dist;        
               
end
        
    
            % PLOT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

num = [1 3];
num2 = [2 4];

for f = 1:2
        
    subplot(4,1,num(f))
    h = histc(all_event_intervals{f},0:1:35);
    bar(h,1,'b')
    axis tight
    xlim([0 35])
    hold on
    if max(h) < 1
        h = 1;
    end
    perc = (round(num_ss_intervals(f)/num_intervals(f)*1000))*.1;
    text(.5,max(h)-(max(h)/3),[num2str(num_ss_intervals(f)) '/' num2str(num_intervals(f)) ' (' num2str(perc) '%) intervals'],'Color','w')
    perc = round(num_ss_cells(1)/length(events)*1000)*.1;
    text(.5,max(h)-1.5*(max(h)/3),[num2str(num_ss_cells(f)) '/' num2str(length(events)) ' (' num2str(perc) '%) cells'],'Color','w')
    text(.5,max(h)-2*(max(h)/3),[num2str(rez(f)) ' ms PSTH bins'],'Color','w')
    set(gca,'box','off')
    set(gca,'Tickdir','out')
    ylabel('Number of Intervals')
    xlabel('Inter Event Interval (msec)')
       
    event_dist = ss_event_dist{f};
    event = event_dist(:,2);
    event = event(find(event >0));
    [h,bins] = histc(event,0:1:max(event));
    
    subplot(4,1,num2(f))
    bar(0:1:max(event),h,'b')
    axis tight
    set(gca,'box','off')
    set(gca,'Tickdir','out')
    ylabel('Number of Cells')
    xlabel('Number of Sub-Stimulus Event Intervals')
        
end

figure

for f = 1:2
    
    subplot(2,1,f)
    [N,X] = hist(all_event_intervals{f},100000);
    plot(X,N,'b')
    set(gca,'box','off')
    set(gca,'Tickdir','out')
    set(gca,'XScale','log')
    xlabel('IEI (msec)')
    ylabel('Count')
    xlim([0 20000])
    
end


for f = 1:2
    for ff = 1:length(events)
        
            ss_intervals = [ss_intervals events(ff).ss_event_times{f}];
            
            if ~isempty(events(ff).ss_event_times{f})
                cell_below(ff) = 1;
            else
                cell_below(ff) = 0;
            end
            
    end
    cell_rez(f,:) = cell_below;
end

fprintf('\n')
fprintf('Cell # | 3ms | 5ms |')
cell = [cell_numbers' cell_rez']
fprintf('Cell # | 3ms | 5ms |')
