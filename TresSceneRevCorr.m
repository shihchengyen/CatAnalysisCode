function [all_tres_scene_rev] = tresscenerevcorr(num_sub_bins,fig)

load jonscellsdata

all_tres_scene_rev = [];

for x = 1:length(jonscells)
        
        sc_num = jonscells(x).scene_changes;
        
        frames = jonscells(x).frame_times;
        
        frames = frames(1:end-1);
        
        sub_bins = [10 6 3 2];

        nsb = find(sub_bins == num_sub_bins);

        tres = jonscells(x).tresscores{nsb};
        
        if sc_num(1) == 1;
        
        else
            sc_num = [1 sc_num];
        end
    
    tres_above = find(tres > 99);
    
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
                
end %for



if exist('fig')
    
    hold on
    bar(sum(all_tres_scene_rev));
    set(gca,'XTickLabel',5:-1:1)    
    xlabel('Scene Change Frames Before TRES Score Frames > 99')
    ylabel('Number of Scene Changes')
    title(['TRES Scores using ',num2str(num_sub_bins) ' Sub Bins'])

end
