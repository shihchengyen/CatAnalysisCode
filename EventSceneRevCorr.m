function [all_event_scene_rev] = eventscenerevcorr(binsize,fig)

load jonscellsdata
load jonsevents

all_tres_scene_rev = [];
bins = [3 5];
eve_ind = find(bins == binsize);

for x = 1:length(jonscells)
        
        sc_num = jonscells(x).scene_changes;
        
        frames = jonscells(x).frame_times;
        
        frames = frames(1:end-1);
                       
        event_times = events(x).all_event_times{eve_ind};
        
        ss_events = event_times(events(x).ss_event_ind{eve_ind});
        
        if isempty(ss_events)
            
            event_scene_rev = [];
            
        else
            
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
                   
            if sc_num(1) == 1;
        
            else
                sc_num = [1 sc_num];
            end
    
               
            event_scene_rev = [];
    
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
                
end %for



if exist('fig')
    
    hold on
    bar(sum(all_event_scene_rev));
    set(gca,'XTickLabel',5:-1:1)    
    xlabel('Scene Change Frames Before Sub-Stimulus Event Intervals')
    ylabel('Number of Scene Changes')
    title(['Events using ',num2str(binsize) ' ms Bins'])

end
