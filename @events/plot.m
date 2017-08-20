function [obj, varargout] = plot(obj,varargin)
%@events/plot Plot function for events object.
%   OBJ = plot(OBJ) creates plot of the events analysis result
%
Args = struct('showTitle',1,'chunkSize',[],'showChunk',0,'GroupPlots',0,'GroupPlotIndex',0, ...
    'Color','b','xlabel',0,'stimInfoDir',['..' filesep '..' ],'linkedZoom',1,'showRaster',0, ...
    'showPSTH',0,'showEvents',0,'showThreshold',0,'showDurations',0,'showProbability',0,'numBins',25, ...
    'XMax',1000);

Args = getOptArgs(varargin,Args,'flags',{'showTitle','showChunk','xlabel','linkedZoom',...
    'showRaster','showPSTH','showEvents','showThreshold','showDurations','showProbability'});    

if isempty(Args.NumericArguments)
    n=1:get(obj,'Number');
else
    n=Args.NumericArguments{1};
end
if isempty(varargin)
    Args.showEvents=1;
end

if ~(Args.showProbability | Args.showDurations)
    eval(['cd ' obj.data.setNames{n}])
    if Args.showRaster
        as = adjspikes('auto');
        raster = as.data.raster';
    end
    if isempty(Args.chunkSize)
        sp=ispikes('auto');
        chunkSize = sp.data.chunkSize;
    else
        chunkSize = Args.chunkSize;
    end
    %%% Load the Stimulus Info %%%%
    pdir = pwd;
    cd(Args.stimInfoDir)
    stimInfo = stiminfo('auto');
    cd(pdir)

    time = [(n-1)*chunkSize ...
        n*chunkSize]*1000;
    rd = stimInfo.data.catInfo.repetition_duration;
    x=rem(time,rd);
    y=floor(time/rd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showPSTH
    hold on
    %%%% Plot the PSTH from the events analysis %%%%
    index = max(obj.data.events.PSTHBins(:,n));
    plot(obj.data.events.PSTHBins(1:index,n),obj.data.events.PSTH(1:index,n),'Color',Args.Color);
    if Args.xlabel
        xlabel('Time (msec)')
    end
    ylabel('Spikes / Second')
    xlim([0 rd])
    %ylim([0 max(obj.data.events.PSTH(:,n))]);
    if Args.showTitle
        title(obj.data.setNames{n})
    end
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showEvents
    hold on
    %index = max(obj.data.events.PSTH(:,n));
    index = 1;
    for r = 1:length(obj.data.events.StartEvent(:,n))
        line([obj.data.events.StartEvent(r,n) obj.data.events.EndEvent(r,n)],[index index],'Color',Args.Color,'LineWidth',3)
    end
    xlim([0 rd])
    %ylim([0 index+3])
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showThreshold
    hold on
    %% Upper Threshold %%%
    line([0 rd],[obj.data.events.Thresholds(n).meanupperthreshold obj.data.events.Thresholds(n).meanupperthreshold],'Color',Args.Color,'LineStyle','-','LineWidth',2)
    %% Lower Threshold %%%
    line([0 rd],[obj.data.events.Thresholds(n).meanlowerthreshold obj.data.events.Thresholds(n).meanlowerthreshold],'Color',Args.Color,'LineStyle','--','LineWidth',2)
    xlim([0 rd])
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showRaster
    NumReps = size(raster,1);
    r = transpose([1:size(raster(1:NumReps,:),1)])*ones(1,size(raster,2));
    hold on
    plot(raster,r+max(obj.data.events.PSTH(:,n))+3,'LineStyle','none','Marker','.','Color',Args.Color);
    ylim([0 max(obj.data.events.PSTH(:,n))+NumReps])
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showDurations
    %%%% Reshape the EventDurations into a single vector
    obj.data.events.EventDurations = reshape(obj.data.events.EventDurations,1,prod(size(obj.data.events.EventDurations)));
    binedges = 0:max(obj.data.events.EventDurations)/Args.numBins:max(obj.data.events.EventDurations);
    n = histcie(obj.data.events.EventDurations,binedges);
    bar(binedges,n,1,Args.Color);
    set(gca,'FontSize',16)
    ylabel('Number of Occurences')
    axis([0-(max(obj.data.events.EventDurations)/Args.numBins)/2 Args.XMax 0 max(n)])
    xlabel('Event Duration (msec)')
    string = ['EVENT DURATIONS -- ' num2str(get(obj,'Number')) ' Cells  Mean: ' num2str(nanmean(obj.data.events.EventDurations(:))) '  Median: ' num2str(nanmedian(obj.data.events.EventDurations(:)))];
    if Args.showTitle
        title(string)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Args.showProbability
    binedges = 0:1/Args.numBins:1;
    n = histcie(obj.data.events.EventProbability,binedges);
    bar(binedges,n,1,Args.Color);
    set(gca,'FontSize',16)
    ylabel('Number of Cells')
    axis([0-(1/Args.numBins)/2 1 0 max(n)])
    xlabel('Probability')
    string = ['EVENT PROBABILITIES -- ' num2str(get(obj,'Number')) ' Cells  Mean: ' num2str(nanmean(obj.data.events.EventProbability(:))) ' Median: ' num2str(nanmedian(obj.data.events.EventProbability(:)))]
    if Args.showTitle
        title(string)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






