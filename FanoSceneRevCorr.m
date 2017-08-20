function [all_fano_scene_rev] = fanoscenerevcorr(fig)

load jonscellsdata

all_fano_scene_rev = [];

for x = 1:length(jonscells)
        
        sc_num = jonscells(x).scene_changes;
        
        frames = jonscells(x).frame_times;
        
        fano = jonscells(x).fanofactor;
        
        if sc_num(1) == 1;
        
        else
            sc_num = [1 sc_num];
        end
    
    fanos_below = find(fano < 1);
    
    fano_scene_rev = [];
    
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
                
end %for



if exist('fig')
    
    bar(sum(all_fano_scene_rev));
    set(gca,'XTickLabel',5:-1:1)    
    xlabel('Scene Change Frames Before FF Frame')
    ylabel('Number of Scene Changes')

end
