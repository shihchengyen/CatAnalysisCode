function fanotresdistribution(sub_bins)

load jonscellsdata
tres_index=[];
fano_index=[];
sub_ind = [10 6 3 2];
sub_ind = find(sub_ind == sub_bins);

for c = 1:length(jonscells)
    
    msc_ind = find(jonscells(c).meanspikecount>=1);
    
    if~isempty(msc_ind)
        
        fano_ind = find(jonscells(c).fanofactor(msc_ind)<1);
        
        if~isempty(fano_ind)
            
            fano_index = [fano_index c];
            
        end
        
        tres = jonscells(c).tresscores{sub_ind};
        tres_ind = find(tres(msc_ind)>99);
        
        if~isempty(tres_ind)
            
            tres_index = [tres_index c];
            
        end
        
    end
end


%%%%%%%%%%FANO%%%%%%%%%%%%%
figure

for c = 1:length(fano_index)
    
    msc_ind = find(jonscells(fano_index(c)).meanspikecount>=1);
    
    fanos = find(jonscells(fano_index(c)).fanofactor(msc_ind)<1);
        
    msc_dist(c) = length(msc_ind);
    fano_dist(c) = length(fanos);
    
    hold on
    plot(c,(length(fanos)/length(msc_ind))*100,'.b')
    text(c,((length(fanos)/length(msc_ind))*100)+1,num2str(jonscells(fano_index(c)).cell_num))
    title('Fano Factor Distribution')
    
end

figure
bar(msc_dist)
hold on
bar(fano_dist,'r')
legend MSC FANO
ylabel('Number of Frames')
xlabel('Cells')
title('FANO FACTOR')


figure

for c = 1:length(tres_index)
    
    msc_ind = find(jonscells(tres_index(c)).meanspikecount>=1);
    
    tres = jonscells(tres_index(c)).tresscores{sub_ind};
    tres = find(tres(msc_ind)>99);
        
    msc_dist(c) = length(msc_ind);
    tres_dist(c) = length(tres);
    
    hold on
    plot(c,(length(tres)/length(msc_ind))*100,'.b')
    text(c,((length(tres)/length(msc_ind))*100)+1,num2str(jonscells(tres_index(c)).cell_num))
    title('TRES Distribution')
    
end

figure
bar(msc_dist)
hold on
bar(tres_dist,'r')
legend MSC TRES
ylabel('Number of Frames')
xlabel('Cells')
title('TRES')
            
        
        
        
        
    
