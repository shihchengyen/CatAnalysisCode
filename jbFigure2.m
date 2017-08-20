function figure2(num_sub_bins)

% Must input the number of sub bins for the tres etc. calculations.
% set constants
% num_sub_bins = 10;

load jonscellsdata

all_fanos_msc=[];
all_tres_msc=[];
all_fanos=[];
all_tres=[];
all_tres_std=[];
all_entropy=[];
msc_cells=[];
msc_ff_cells=[];
msc_tres_cells=[];
tres_cells=[];
all_fanos_msc_1=[];
all_msc=[];

for f = 1:length(jonscells)
	
	msc_ind = find(jonscells(f).meanspikecount>=1);
	
	sub_bins = [10 6 3 2];
	
	sub_ind = find(sub_bins == num_sub_bins);
	
	tres_ind = find(jonscells(f).tresscores{sub_ind}>99);
	
	if ~isempty(tres_ind)
		
		tres_cells = [tres_cells f];
		
	end
	
	tres_std = jonscells(f).tres_std{sub_ind};
	
	entropy = jonscells(f).entropy{sub_ind};
		
	if ~isempty(msc_ind)
		
		msc_cells = [msc_cells f];
		
		fano = jonscells(f).fanofactor(msc_ind);
		
		fano_ind = find(fano<1);
		
		fano = fano(fano_ind);
		
		entropy = entropy(msc_ind);
		
		entropy = entropy(fano_ind);
		
		tres = jonscells(f).tresscores{sub_ind};
		
		tres = find(tres(msc_ind)>99);
			   
		tres_std = tres_std(msc_ind);
							 
		if ~isempty(fano)
			
			msc_ff_cells = [msc_ff_cells f];
			
		elseif ~isempty(tres)
			
			msc_tres_cells = [msc_tres_cells f];
			
		end
	   
		all_tres_std = [all_tres_std tres_std];
		
		all_entropy = [all_entropy entropy];
		
		all_fanos_msc_1 = [all_fanos_msc_1 fano];
		
	end
	
	all_fanos_msc = [all_fanos_msc jonscells(f).fanofactor(msc_ind)];
	
	all_msc = [all_msc jonscells(f).meanspikecount(msc_ind)];
	
	tres = jonscells(f).tresscores{sub_ind};
	
	all_tres_msc = [all_tres_msc tres(msc_ind)];
	
	all_fanos = [all_fanos jonscells(f).fanofactor];
	
	all_tres = [all_tres jonscells(f).tresscores{sub_ind}];
	
end

fano_ind = find(all_fanos == 100);

all_fanos(fano_ind) = NaN;

fano_ind = find(all_fanos_msc < 1);

Number_Cells = length(jonscells)

Number_Cells_Frame_MSC1 = length(msc_cells)

Number_Cells_Frame_MSC1_FF1 = length(msc_ff_cells)

Number_Cells_Frame_TRES100 = length(tres_cells)

Number_Cells_Frame_MSC1_TRES100 = length(msc_tres_cells)

Number_Frames = length(all_fanos)

Number_Frames_MSC1 = length(all_fanos_msc)

Number_Frames_MSC1_FF1 = length(find(all_fanos_msc < 1))

Number_Frames_MSC1_TRES100 = length(find(all_tres_msc > 99))

Number_Frames_MSC1_FF1_TRES100 = length(find(all_tres_msc(fano_ind) > 99))

Number_Frames_FF1 = length(find(all_fanos < 1))

%%%%%%%%%%%%%%%%%% PANEL A %%%%%%%%%%%%%%%%%%%%%%%%

figure

for f = 1:length(jonscells)
    
    hold on
    plot(jonscells(f).meanspikecount,jonscells(f).fanofactor,'k.')
    ylim([0 10])
    xlim([0 5])
    xlabel('Mean Count [spikes]')
    ylabel('FF')
    
end
hold on
line([1 1],[0 10],'Color','r','LineWidth',2)
line([0 5],[1 1],'Color','r','LineWidth',2)


%%%%%%%%%%%%%%%%%% PANEL B %%%%%%%%%%%%%%%%%%%%%%%%

figure

for f = 1:length(jonscells)
    
    hold on
    plot(jonscells(f).meanspikecount,jonscells(f).fanofactor,'k.')
    ylim([0 1])
    xlim([0 5])
    xlabel('Mean Count [spikes]')
    ylabel('FF')
    
end
hold on
line([1 1],[0 1],'Color','r','LineWidth',2)

%%%%%%%%%%%%%%%%%% PANEL C %%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:4
    
    bin_count = [];
    
    for c = 1:length(jonscells)

        msc_ind = find(jonscells(c).meanspikecount>=1);
        tres = jonscells(c).tresscores{s};
        bin_count = [bin_count length(find(tres(msc_ind)>99))];
        
    end
    
    frame_count_tres_sig(s) = sum(bin_count);
    
end

figure
bar((frame_count_tres_sig/Number_Frames_MSC1)*100,1,'b')
ylabel('% occurences')
set(gca,'XTickLabel',[16 11 7 3])

%%%%%%%%%%%%%%%%%% PANEL D %%%%%%%%%%%%%%%%%%%%%%%%

figure
[N,X] = histc(all_fanos_msc,0:.5:6);
bar(0:.5:6,N,1,'r')
hold on
[N,X] = histc(all_fanos_msc(find(all_tres_msc>99)),0:.5:6);
bar(0:.5:6,N,1,'b')
xlim([0 6])
xlabel('FF')
ylabel('# occurences')


%%%%%%%%%%%%%%%%%% PANEL E %%%%%%%%%%%%%%%%%%%%%%%%

figure

plot(all_tres_std,all_fanos_msc,'k.')
hold on
line([0 10],[1 1],'Color','r')
line([1 1],[0 10],'Color','r')
xlim([0 10])
ylim([0 5])
xlabel('TREStd(10)')
ylabel('FF')

%%%%%%%%%%%%%%%%%% PANEL F %%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(all_entropy,all_fanos_msc_1,'k.')
xlabel('ENT(10)')
ylabel('FF')


%%%%%%%%%%%%%%%%%% PANEL C %%%%%%%%%%%%%%%%%%%%%%%%

tres_ind = find(all_tres == 0);

all_tres(tres_ind) = NaN;

bin = 0:2:100;

figure
[N] = histc(all_tres,bin);
bar(bin,N,1,'r')
hold on
[N] = histc(all_tres_msc,bin);
bar(bin,N,1,'b')
xlabel('TRES Scores')
axis tight


figure

plot(all_msc,all_fanos_msc,'k.')
ylim([0 10])
xlim([0 5])
xlabel('Mean Count [spikes]')
ylabel('FF')
line([1 1],[0 10],'Color','r','LineWidth',2)
line([0 5],[1 1],'Color','r','LineWidth',2)



