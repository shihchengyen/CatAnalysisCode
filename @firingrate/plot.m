function [obj, varargout] = plot(obj,varargin)
%@firingrate/plot Plot function for the firingrate object.
%   OBJ = plot(OBJ) creates a firing rate plot of the neuronal
%   response.

Args = struct('LabelsOFF',0,'Color',[0 0 1],'SubPlot',0,...
    'numBins',100,'HistRate',0,'Normalize',0,'TossZeros',0,'loglog',0,...
    'semilogy',0,'Fit',0,'BootStrap',0,'BootDist',0,'Surrogates',0,...
    'PlotNumber',0,'PSTH',0,'HoldON',0,'Hetero',0,'Mean',0,'FitDist',0);

Args.flags = {'LabelsOFF','SubPlot','Normalize','TossZeros',...
    'loglog','semilogy','Fit','BootStrap','Surrogates',...
    'PSTH','BootDist','HoldON','Hetero','Mean','FitDist'};

Args = getOptArgs(varargin,Args);

if isempty(Args.NumericArguments)
    n = 1:obj.data.numSets;
else
    n = Args.NumericArguments{1};
end

%%%%% Plot the Data %%%%%%%%%%
numSets = obj.data.numSets;
%%%% Ploting Arguments %%%%
if Args.SubPlot
    for i = 1:numSets
        if Args.Fit
            y = obj.data.Fit(i).N;
            % x = obj.data.Fit(i).bins(2:end);
            x = obj.data.Fit(i).bins;
        else
            fr_values = obj.data.firingRate(:,i);
            if Args.TossZeros
                fr_values = fr_values(find(fr_values>0));
            end
            if Args.Normalize
                %fr_values = fr_values/nanmean(fr_values);
                fr_values = fr_values/nanmax(fr_values);
            end
            max_value = nanmax(fr_values); min_value = nanmin(fr_values);
            bins = (min_value:(max_value-min_value)/Args.numBins:max_value)';
            [y] = histcie(fr_values,bins,'DropLast'); x = bins(2:end);
        end
        nptSubplot(numSets,i);
        if(~isempty(x))
            if Args.Fit
                non_zeros = find(y>0);
                y = log10(y(non_zeros));
                x = x(non_zeros);
                Args.FitName1 = 'poly1'; % Exp Fit
                Args.FitName2 = fittype('a + b*log(x)'); % Power Fit
                nptSubplot(numSets,i);
                plot(x,y,'Color',Args.Color,'LineStyle','-','Marker','*')
                xlim([min(x) max(x)])
                xx = get(gca,'XTickLabel');
                xxsize = size(xx,1);
                clear xx1
                xx1{xxsize} = xx(xxsize,:);
                set(gca,'XTickLabel',xx1);
                a = axis; %% Get the axis
            else
                plot(x,y,'Color',Args.Color,'LineStyle','-','Marker','*')
                xlim([min(x) max(x)])
                if Args.loglog
                    set(gca,'XScale','log','YScale','log')
                elseif Args.semilogy
                    set(gca,'YScale','log')
                end
                a = axis; %% Get the axis
            end
            if Args.Fit
                nptSubplot(numSets,i);
                [fitobj,goodness,output] = fit(x,y,Args.FitName1);
                exp = (round((goodness.adjrsquare*100)))/100;
                hold on; plot(fitobj); legend off
                [fitobj,goodness,output] = fit(x,y,Args.FitName2);
                power = (round((goodness.adjrsquare*100)))/100;
                hold on; plot(fitobj); legend off
                %title(['Exp: ',num2str(exp) '  Power: ',num2str(power)])
                ylabel('')
                xlabel('')
            end
        end
        axis(a)
        hold off
    end
elseif(Args.PSTH)
    % mcounts = obj.data.meancounts * 1000/obj.data.binsize;
    h = stairs(obj.data.timebins,obj.data.meancounts);
    set(h,'Color',Args.Color);
elseif Args.BootStrap
    ec=0;
    pc=0;
    for ii = 1:length(n)
        hold on
        Mean_R_Values = obj.data.Fit(n(ii)).BootStrap_Mean_R_Values;
        Data_R_Values = obj.data.Fit(n(ii)).Firing_Rate_R_Values;
        H_Value = obj.data.Fit(n(ii)).BootStrap_H_Value;
        if ~isempty(Data_R_Values)
            if Args.BootDist
                R_Dist = obj.data.Fit(n(ii)).BootStrap_R_Values;
                for rr = 1:size(R_Dist,1)
                    plot(R_Dist(rr,2),R_Dist(rr,1),'*b')
                end
                plot(Data_R_Values(2),Data_R_Values(1),'or','LineWidth',2)
                line([0 1],[0 1],'Color','k')
            else
                plot(Mean_R_Values(2),Mean_R_Values(1),'ok','LineWidth',2)
                if H_Value == 1
                    if max(Mean_R_Values) > .75
                        plot(Mean_R_Values(2),Mean_R_Values(1),'or','LineWidth',2)
                    end
                end
            end
            if H_Value == 1
                if Mean_R_Values(1) > .75
                    if Mean_R_Values(1) > Mean_R_Values(2)
                        ec=ec+1;
                    end
                end
                if Mean_R_Values(2) > .75
                    if Mean_R_Values(2) > Mean_R_Values(1)
                        pc=pc+1;
                    end
                end
            end
        end
    end
    line([0 1],[0 1],'Color','k')
    title(['Power: ',num2str(pc), ' / Exp: ',num2str(ec)])
elseif(Args.Hetero)
    data = heterogeneity(obj,varargin{:});
    bar(data','grouped')
elseif Args.FitDist
    ec=0;
    pc=0;
    for ii = 1:length(n)
        hold on
        Data_R_Values = obj.data.Fit(n(ii)).Firing_Rate_R_Values;
        if ~isempty(Data_R_Values)
            plot(Data_R_Values(2),Data_R_Values(1),'ok','LineWidth',2)
            if max(Data_R_Values) > .75
                plot(Data_R_Values(2),Data_R_Values(1),'or','LineWidth',2)
            end
            if max(Data_R_Values) > .75
                if Data_R_Values(1) > Data_R_Values(2)
                    ec=ec+1;
                elseif Data_R_Values(2) > Data_R_Values(1)
                    pc=pc+1;
                end
            end
        end
    end
    line([0 1],[0 1],'Color','k')
    title(['Power: ',num2str(pc), ' / Exp: ',num2str(ec)])
else
    fr_values = obj.data.firingRate(:,n);
    if Args.Normalize
        fr_values = fr_values/nanmax(fr_values);
    end
    max_value = nanmax(fr_values);
    min_value = nanmin(fr_values);
    if Args.HistRate > 0
        bins = (min_value:Args.HistRate:max_value+Args.HistRate)';
    else
        bins = (min_value:(max_value-min_value)/Args.numBins:max_value)';
    end
    [y] = histcie(fr_values,bins,'DropLast'); x = bins(2:end);
    % Remove the zero values to create a continuous firingrate function
    non_zeros = find(y>0);
    y = y(non_zeros);
    x = x(non_zeros);
    if Args.semilogy
        y = log(y);
    elseif Args.loglog
        x = log(x);
    end
    if Args.HoldON
        hold on
    end
    plot(x,y,'Color',Args.Color,'LineStyle','-','Marker','*')
    xlim([min(x) max(x)])
    a = axis; %% Get the axis
    if Args.Fit
        hold on
        Args.FitName1 = 'poly1'; % Exp Fit
        Args.FitName2 = fittype('a + b*log(x)'); % Power Fit
        [fitobj,goodness,output] = fit(x,y,Args.FitName1);
        exp = (round((goodness.adjrsquare*100)))/100;
        plot(fitobj,'g');
        [fitobj,goodness,output] = fit(x,y,Args.FitName2);
        power = (round((goodness.adjrsquare*100)))/100;
        plot(fitobj);
        title(['E: ',num2str(exp) '  P: ',num2str(power)])
        legend('Data','Exponential','Power')
        ylabel('')
        xlabel('')
    end
    axis(a)
    hold off
end

if Args.SubPlot
    nptSubplot(numSets,'BottomLeft');
    if Args.Fit
        ylabel('Log (Number of Occurences)')
    else
        if Args.loglog
            ylabel('Log (Number of Occurences)')
        elseif Args.semilogy
            ylabel('Log (Number of Occurences)')
        else
            ylabel('Number of Occurences')
        end
    end
    if Args.loglog
        xlabel('Log (Spike Rate or Spike Count)')
    else
        xlabel('Spike Rate or Spike Count')
    end
elseif(Args.PSTH)
    xlabel('Time (ms)')
    ylabel('Rate (Hz)')
elseif Args.BootStrap
    xlabel('R^2 Power')
    ylabel('R^2 Exponential')
elseif(Args.Hetero)
    xlabel('Recording Sites')
    if(Args.Mean)
        ylabel('Mean Firing Rate (spikes/s)')
    else
        ylabel('Peak Firing Rate (spikes/s)')
    end
elseif Args.FitDist
    xlabel('R^2 Power')
    ylabel('R^2 Exponential')
else
    if ~Args.LabelsOFF
        if Args.Fit
            ylabel('Log (Number of Occurences)')
        else
            if Args.loglog
                ylabel('Log (Number of Occurences)')
            elseif Args.semilogy
                ylabel('Log (Number of Occurences)')
            else
                ylabel('Number of Occurences')
            end
        end
        if Args.loglog
            xlabel('Log (Spike Rate or Spike Count)')
        else
            xlabel('Spike Rate or Spike Count')
        end
    end
end

