function [obj, varargout] = plot(obj,varargin)
%@jointevents/plot Plot function for jointevents object.
%   OBJ = plot(OBJ) creates plot of the jointevents analysis result
%
Args = struct('showTitle',1,'GroupPlots',0,'GroupPlotIndex',0, ...
    'Color','b','xlabel',0,'stimInfoDir',['..' filesep '..' ],'linkedZoom',1, ...
    'showDurations',0,'showProbability',0,'showSimilarityIndex',0,'showSpikeCountCorrelations',0,...
    'showSpikeCountCorrelationsDurations',0,'showJointEventsPSTHCorrelations',0,...
    'showJointEventsPSTHCorrelationDistributions',0,'numBins',25,'Significance',0.05,...
    'showShiftedSpikeCountCorrelations',0,'showJointEvents',0,'showJointPSTHScatterPlot',0,'showJEDiffs',0);

Args = getOptArgs(varargin,Args,'flags',{'showTitle','xlabel','linkedZoom','showDurations',...
        'showProbability','showSimilarityIndex','showSpikeCountCorrelations',...
        'showSpikeCountCorrelationsDurations','showJointEventsPSTHCorrelations',...
        'showJointEventsPSTHCorrelationDistributions','showShiftedSpikeCountCorrelations',...
        'showJointEvents','showJointPSTHScatterPlot','showJEDiffs'});

% if there are numeric arguments, figure out who it is for
if(~isempty(Args.NumericArguments))
    Cum = 0;
    n = Args.NumericArguments{1};
    Joint_Overlap_Durations = obj.data.jointevents.Joint_Overlap_Durations(n,:);
    Joint_Overlap_Probability = obj.data.jointevents.Joint_Overlap_Probability(n);
    Joint_Overlap_Expected = obj.data.jointevents.Expected_Overlap(n);
    Similarity_Index = obj.data.jointevents.Similarity_Index(n);
    Spike_Count_Correlations = obj.data.jointevents.Spike_Count_Correlations(n,:);
    Spike_Count_Significances = obj.data.jointevents.Spike_Count_Significances(n,:);
    Shift_Spike_Count_Correlations = obj.data.jointevents.Shift_Spike_Count_Correlations(n,:);
    Shift_Spike_Count_Significances = obj.data.jointevents.Shift_Spike_Count_Significances(n,:);
    PSTH_Correlation_Coefficents_Dist = obj.data.jointevents.RandomCC_Values(:,n);
    PSTH_Correlation_Coefficents_Dist_Sig = obj.data.jointevents.RandomCC_Values_Significance_Value(n);
    PSTH_Correlation_Coefficents = obj.data.jointevents.Correlation_Coefficent_Value(n);
    PSTH_Correlation_Significances = obj.data.jointevents.Correlation_Significance_Value(n);
    Start_Overlap =  obj.data.jointevents.Joint_Overlap_Start(n,:);
    End_Overlap =  obj.data.jointevents.Joint_Overlap_End(n,:);
    if n==1
        Rasters = obj.data.jointevents.Rasters(1:obj.data.jointevents.NumReps(n),:);
        PSTH = obj.data.jointevents.PSTH(1:2,:);
        Start = obj.data.jointevents.StartEvent(:,1:2);
        End = obj.data.jointevents.EndEvent(:,1:2);
    else
        Rasters = obj.data.jointevents.Rasters(obj.data.jointevents.NumReps(n-1)+1:obj.data.jointevents.NumReps(n),:);
        PSTH = obj.data.jointevents.PSTH((n*2)-1:n*2,:);
        Start = obj.data.jointevents.StartEvent(:,(n*2)-1:n*2);
        End = obj.data.jointevents.EndEvent(:,(n*2)-1:n*2);
    end
    UpperThresholds = obj.data.jointevents.UpperThresholds(n,:);
    LowerThresholds = obj.data.jointevents.LowerThresholds(n,:);
else 
    Cum = 1;
    Joint_Overlap_Durations = reshape(obj.data.jointevents.Joint_Overlap_Durations,1,[]);
    Joint_Overlap_Probability = obj.data.jointevents.Joint_Overlap_Probability;
    Joint_Overlap_Expected = reshape(obj.data.jointevents.Expected_Overlap,1,[]);
    Similarity_Index = obj.data.jointevents.Similarity_Index;
    Spike_Count_Correlations = reshape(obj.data.jointevents.Spike_Count_Correlations,1,[]);
    Spike_Count_Significances = reshape(obj.data.jointevents.Spike_Count_Significances,1,[]);
    Shift_Spike_Count_Correlations = reshape(obj.data.jointevents.Shift_Spike_Count_Correlations,1,[]);
    Shift_Spike_Count_Significances = reshape(obj.data.jointevents.Shift_Spike_Count_Significances,1,[]);
    PSTH_Correlation_Coefficents_Dist = obj.data.jointevents.RandomCC_Values;
    PSTH_Correlation_Coefficents_Dist_Sig = obj.data.jointevents.RandomCC_Values_Significance_Value;
    PSTH_Correlation_Coefficents = reshape(obj.data.jointevents.Correlation_Coefficent_Value,1,[]);
    PSTH_Correlation_Significances = reshape(obj.data.jointevents.Correlation_Significance_Value,1,[]);
    PSTH = obj.data.jointevents.PSTH;
    Start = obj.data.jointevents.StartEvent;
    End = obj.data.jointevents.EndEvent;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showDurations
    %%%% Plot the Joint Event Durations %%%%
    bins = 0:nanmax(Joint_Overlap_Durations)/Args.numBins:nanmax(Joint_Overlap_Durations);
    median_dur =  nanmedian(Joint_Overlap_Durations);
    mean_dur =  nanmean(Joint_Overlap_Durations);
    [N] = histcie(Joint_Overlap_Durations,bins);
    bar(bins,N,'histc')
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',Args.Color,'EdgeColor','k')
    set(gca,'FontSize',16)
    xlabel('Overlap Durations (msec)')
    ylabel('Number of Occurrences')
    xlim([0 nanmax(Joint_Overlap_Durations)+10])
    if Cum
        title(['Overlap Durations for ',num2str(obj.data.numSets) ' Cell Pairs' '  Median: ',num2str(median_dur) '  Mean: ',num2str(mean_dur)])
    else
        title([obj.data.setNames{n} '  Median: ',num2str(median_dur) '  Mean: ',num2str(mean_dur)])
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showProbability
    %%%% Plot the Event Probablity %%%%
    bins = 0:nanmax(Joint_Overlap_Probability)/Args.numBins:nanmax(Joint_Overlap_Probability);
    median_prob =  nanmedian(Joint_Overlap_Probability);
    mean_prob =  nanmean(Joint_Overlap_Probability);
    [N1] = histcie(Joint_Overlap_Probability,bins);
    bar(bins,N1,'histc')
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',Args.Color,'EdgeColor','k')
    set(gca,'FontSize',16)
    xlabel('Joint Probability')
    ylabel('Number of Occurrences')
    if Cum
        title(['Population Probabilities for ',num2str(obj.data.numSets) ' Cell Pairs' '  Median: ',num2str(median_prob) '  Mean: ',num2str(mean_prob)])
        xlim([0 nanmax(Joint_Overlap_Probability)])
    else
        title([obj.data.setNames{n} '  Median: ',num2str(median_prob) '  Mean: ',num2str(mean_prob)])
        xlim([0 nanmax(Joint_Overlap_Probability)])
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showSimilarityIndex
    %%%% Plot the Similarity Index %%%%
    bins = -1:.05:1; 
    ind = find(Similarity_Index>0);
    median_sim =  nanmedian(Similarity_Index(ind));
    med_sim =  nanmedian(Similarity_Index);
    [N] = histcie(Similarity_Index,bins);
    bar(bins,N,'histc')
    set(gca,'FontSize',7)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',Args.Color,'EdgeColor','k')
    xlabel('Similarity Index')
    ylabel('Number')
    set(gca,'TickDir','out')
    %set(gca,'Box','off')
    xlim([-1.05 1.05])
    if Cum
        title(['Intragroup Cell Pairs (n=' ,num2str(obj.data.numSets) ')  Median(all): ',num2str(med_sim) '  Median(>0): ',num2str(median_sim)])
    else
        title([obj.data.setNames{n} '  Median: ',num2str(median_sim) '  Mean: ',num2str(mean_sim)])
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showSpikeCountCorrelations
    %%%% Plot the Spike Count Correlations %%%%
    index = isnan(Spike_Count_Correlations); % Take out the NaNs
    index = find(index == 1);
    Spike_Count_Correlations(index) = [];
    Spike_Count_Significances(index) = [];
    if ~isempty(Spike_Count_Correlations)
        index = find(Spike_Count_Significances<Args.Significance);
        Significant_Spike_Count_Correlations = Spike_Count_Correlations(index);
        Perc_Sign = (length(Significant_Spike_Count_Correlations)/length(Spike_Count_Correlations))*100;
        bins = -1:2/Args.numBins:1;
        index = find(Significant_Spike_Count_Correlations>=0);
        median_scc =  median(Significant_Spike_Count_Correlations(index));
        mean_scc =  mean(Significant_Spike_Count_Correlations(index));
        [N1] = histcie(Spike_Count_Correlations,bins);
        bar(bins,N1,1)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','none','EdgeColor','k')
        hold on
        [N2] = histcie(Significant_Spike_Count_Correlations,bins);
        bar(bins,N2,1)
        set(gca,'FontSize',16)
        %bar(bins,[N2 N1],'stack')
        xlabel('Correlation Coefficent')
        ylabel('Number of Occurrences')
        xlim([-1 1])
        hold off
        if Cum
            title(['Population ',num2str(obj.data.numSets) ' Cell Pairs' '  Median: ',num2str(median_scc) '  Mean: ',num2str(mean_scc) ' Perc. Sign.: ',num2str(Perc_Sign)])
        else
            title([obj.data.setNames{n} '  Median: ',num2str(median_scc) '  Mean: ',num2str(mean_scc)])
            Joint_Overlap_Indicies = find(Spike_Count_Correlations<0);
            Spike_Count_Correlations = Spike_Count_Correlations(Joint_Overlap_Indicies);
            if ~isempty(Joint_Overlap_Indicies)
                Joint_Overlap_Indicies
                Spike_Count_Correlations
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showSpikeCountCorrelationsDurations
    index = find(Spike_Count_Significances<Args.Significance);
    Significant_Spike_Count_Correlations = Spike_Count_Correlations(index);
    Significant_Durations = Joint_Overlap_Durations(index);
    index = find(Spike_Count_Significances>=Args.Significance);
    Non_Significant_Spike_Count_Correlations = Spike_Count_Correlations(index);
    Non_Significant_Durations = Joint_Overlap_Durations(index);
    plot(Non_Significant_Spike_Count_Correlations,Non_Significant_Durations,'k*','LineWidth',3)   
    hold on
    plot(Significant_Spike_Count_Correlations,Significant_Durations,[Args.Color '*'],'LineWidth',3)  
    hold off
    set(gca,'FontSize',16)
    xlabel('Correlation Coefficent')
    ylabel('Overlap Durations (msec)')
    set(gca,'FontSize',16)
    xlim([-1 1])
    if Cum
        title(['Population ',num2str(obj.data.numSets) ' Cell Pairs'])
    else
        title(obj.data.setNames{n})
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showShiftedSpikeCountCorrelations
    %%%% Plot the Spike Count Correlations %%%%
    bins = -1:2/Args.numBins:1;
    [N1] = histcie(Shift_Spike_Count_Correlations,bins);
    bar(bins,N1,'histc')
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor','none','EdgeColor','k')
    hold on
    [N2] = histcie(Spike_Count_Correlations,bins);
    bar(bins,N2,'histc')
    set(gca,'FontSize',16)
    xlabel('Correlation Coefficent')
    ylabel('Number of Occurrences')    
    xlim([-1 1])
    if Cum
        title(['Population ',num2str(obj.data.numSets) ' Cell Pairs  Mean: ',num2str(nanmean(Spike_Count_Correlations)) '  Mean Shifted: ',num2str(nanmean(Shift_Spike_Count_Correlations))])
    else
        title(obj.data.setNames{n})
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showJointEventsPSTHCorrelations
    index = find(PSTH_Correlation_Significances<Args.Significance);
    Significant_Spike_Count_Correlations = PSTH_Correlation_Coefficents(index);
    bins = -1:2/Args.numBins:1;
    index = find(Significant_Spike_Count_Correlations>=0);
    median_scc =  nanmedian(Significant_Spike_Count_Correlations(index));
    mean_scc =  nanmean(Significant_Spike_Count_Correlations(index));
    [N] = histcie(PSTH_Correlation_Coefficents,bins);
    bar(bins,N,'histc');hold on    
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','none','EdgeColor','k')
    [N] = histcie(Significant_Spike_Count_Correlations,bins);
    bar(bins,N,'histc');hold off
    set(gca,'FontSize',16)
    xlabel('Correlation Coefficent')
    ylabel('# of Occurrences')
    xlim([-1 1])
    if Cum
        title(['Population ',num2str(obj.data.numSets) ' Cell Pairs' '  Median: ',num2str(median_scc) '  Mean: ',num2str(mean_scc)])
    else
        title([obj.data.setNames{n} '  Median: ',num2str(median_scc) '  Mean: ',num2str(mean_scc)])
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showJointEventsPSTHCorrelationDistributions
    
    bins = -1:2/Args.numBins:1;
    [N] = histcie(PSTH_Correlation_Coefficents_Dist,bins);
    bar(bins,N,'histc'); hold on  
    plot(PSTH_Correlation_Coefficents,10,'*r','LineWidth',3)
    set(gca,'FontSize',16)
    xlabel('Correlation Coefficent')
    ylabel('# of Occurrences')
    xlim([min([PSTH_Correlation_Coefficents_Dist;PSTH_Correlation_Coefficents]) max([PSTH_Correlation_Coefficents_Dist;PSTH_Correlation_Coefficents])])
    hold off
    sig = (round(PSTH_Correlation_Significances*100))/100;
    title(['Significance: ' num2str(sig) '  Percentile Sig: ' num2str(PSTH_Correlation_Coefficents_Dist_Sig)])
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showJointEvents  
    NumReps = size(Rasters,1)/2;
    r = transpose([1:size(Rasters(1:NumReps,:),1)])*ones(1,size(Rasters,2));
    plot(Rasters(1:NumReps,:),r,'LineStyle','none','Marker','.','Color','b');  
    hold on
    plot(Rasters(NumReps+1:end,:),r,'LineStyle','none','Marker','.','Color','r'); 
    index = find(Start_Overlap>0);
    for e = 1:length(index);
        line([Start_Overlap(e) End_Overlap(e)],[-2 -2],'LineWidth',2,'Color','k')
        line([Start_Overlap(e) End_Overlap(e)],[NumReps+2 NumReps+2],'LineWidth',2,'Color','k')
        spc = round(Spike_Count_Correlations(e)*100)/100;
        text((Start_Overlap(e)+(End_Overlap(e)-Start_Overlap(e))/2),-12,num2str(spc),'LineWidth',2)
        text((Start_Overlap(e)+(End_Overlap(e)-Start_Overlap(e))/2),NumReps+6,num2str(spc),'LineWidth',2)
    end   
    for e = 1:size(Start,1)
        if Start(e,1) >= 0
            line([Start(e,1) End(e,1)],[-4 -4],'LineWidth',2,'Color','b')
        end
        if Start(e,2) >= 0
            line([Start(e,2) End(e,2)],[-6 -6],'LineWidth',2,'Color','r')
        end
    end    
    stairs(0:1:(size(PSTH,2)-1),PSTH(1,:)-max2(PSTH),'b')
    line([0 (size(PSTH,2)-1)],[UpperThresholds(1)-max2(PSTH) UpperThresholds(1)-max2(PSTH)],'Color','b')
    line([0 (size(PSTH,2)-1)],[LowerThresholds(1)-max2(PSTH) LowerThresholds(1)-max2(PSTH)],'Color','b','LineStyle','--')
    stairs(0:1:(size(PSTH,2)-1),PSTH(2,:)-max2(PSTH),'r')
    line([0 (size(PSTH,2)-1)],[UpperThresholds(2)-max2(PSTH) UpperThresholds(2)-max2(PSTH)],'Color','r')
    line([0 (size(PSTH,2)-1)],[LowerThresholds(2)-max2(PSTH) LowerThresholds(2)-max2(PSTH)],'Color','r','LineStyle','--')
    ylim([-max2(PSTH)-1 NumReps+12])
    ylabel('Repetitions')
    xlabel('Time (msec)')
    hold off
    title([obj.data.setNames{n} '   Similarity Index: ',num2str(Similarity_Index) '  PSTH Correlation Coefficient: ',num2str(PSTH_Correlation_Coefficents)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showJointPSTHScatterPlot
    plot(PSTH(1,:),PSTH(2,:),Args.Color)
    if ~isempty(Start_Overlap)
        index = find(Start_Overlap>0);
        if ~isempty(index)
            Start_Overlap = Start_Overlap(index);
            End_Overlap = End_Overlap(index);
            for e = 1:length(Start_Overlap)
                hold on
                plot(PSTH(1,Start_Overlap(e):End_Overlap(e)),PSTH(2,Start_Overlap(e):End_Overlap(e)),'r')
            end
        end
    end
    line([0 max2(PSTH)],[0 max2(PSTH)],'Color','k')
    xlabel('Firing Rate')
    ylabel('Firing Rate')
    title([obj.data.setNames{n} '/  Similarity Index: ',num2str(Similarity_Index) '/  PSTH Correlation Coefficient: ',num2str(PSTH_Correlation_Coefficents)])
    xlim([0 max2(PSTH)])
    ylim([0 max2(PSTH)])
    hold off
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showJEDiffs
    differences=[];
    for n = 1:get(obj,'Number')
        PSTH = obj.data.jointevents.PSTH((n*2)-1:n*2,:);
        Start = obj.data.jointevents.StartEvent(:,(n*2)-1:n*2);
        End = obj.data.jointevents.EndEvent(:,(n*2)-1:n*2);
        index = length(find(Start(1,:)>0));
        if index==2
            for e = 1:length(find(Start(:,1)>0))
                differences = [differences abs(diff([PSTH(1,Start(e,1):End(e,1));PSTH(2,Start(e,1):End(e,1))]))];
            end
            for e = 1:length(find(Start(:,2)>0))
                differences = [differences abs(diff([PSTH(1,Start(e,2):End(e,2));PSTH(2,Start(e,2):End(e,2))]))];
            end  
        end
    end
    bins = 0:max(differences)/Args.numBins:max(differences);
    [N] = histcie(differences,bins);
    bar(bins,N,1,Args.Color)
    xlim([0 max(differences)])
    xlabel('Difference in Firing Rate')
    ylabel('Number of Occurances')
end