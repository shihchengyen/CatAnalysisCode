function [obj, varargout] = plot(obj,varargin)
%@adjspikespth/plot Plot function for adjspikespth object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.
%

Args = struct('showTitle',1,'chunkSize',[],'showChunk',0, ...
    'Color','b','xlabel',0,'stimInfoDir',['..' filesep '..' ],...
    'linkedZoom',1,'plotType','line');

Args = getOptArgs(varargin,Args,'flags',{'showTitle','showChunk','xlabel','linkedZoom'});

if isempty(Args.NumericArguments)
    n=1;
else
    n=Args.NumericArguments{1};
end
cd(obj.data.setNames{n})
% cwd = pwd;
% eval(['cd ' Args.stimInfoDir]);
stimInfo = stiminfo('auto');
% eval(['cd ' pwd]);
if isempty(Args.chunkSize)
    sp=ispikes('auto');
    chunkSize = sp.data.chunkSize;
else
    chunkSize = Args.chunkSize;
end

time = [(n-1)*chunkSize ...
    n*chunkSize]*1000;
rd = stimInfo.data.catInfo.repetition_duration;
x=rem(time,rd);
y=floor(time/rd);

%Plot the PSTH
hold on
if strcmp(Args.plotType,'stairs')
    stairs(obj.data.plotVector(n,:),obj.data.psth(n,:),Args.Color);
else
    plot(obj.data.plotVector(n,:),obj.data.psth(n,:),'Color',Args.Color);
end
if Args.xlabel
    xlabel('Time (seconds)','FontSize',18)
end
ylabel('Spikes / Second','FontSize',18)
set(gca,'FontSize',16)
xlim([0 obj.data.plotVector(end)])
set(gca,'XTick',[0:2000:stimInfo.data.catInfo.repetition_duration+stimInfo.data.iniInfo.spontaneous_activity_duration_in_milliseconds+1000])
set(gca,'XTickLabel',[0:2:100])
xlim([-10 stimInfo.data.catInfo.repetition_duration+stimInfo.data.iniInfo.spontaneous_activity_duration_in_milliseconds+10])

if Args.showChunk
    a=axis;
    if x(2)>x(1)
        h=patch([x(1) x(1) x(2) x(2)],([y(1)+.5  y(1)+1.5  y(1)+1.5  y(1)+.5]*(abs(a(3)-a(4))/100)),'b', 'FaceAlpha', .3);
    else
        h=patch([x(1) x(1) a(2) a(2)],([y(1)+.5  y(1)+1.5  y(1)+1.5  y(1)+.5]*(abs(a(3)-a(4))/100)),'b', 'FaceAlpha', .3);
        h=[h patch([0 0 x(2) x(2)],([y(2)+.5  y(2)+1.5  y(2)+1.5  y(2)+.5]*(abs(a(3)-a(4))/100)),'b', 'FaceAlpha', .3)];
    end
end
hold off


