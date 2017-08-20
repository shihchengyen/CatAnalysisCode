function r = get(s,varargin)
%ISPIKES/GET Returns object properties
%   VALUE = GET(OBJ,PROP_NAME) returns an object 
%   property. PROP_NAME can be one of the following:
%      'SessionName' - name of session.
%      'Channel' - signal number inside streamer file.
%      'Duration'
%      'MinDuration'
%      'TotalClusters'
%      'MeansPresent'
%      'TrialMean','Trial',t
%      'TrialThreshold','Trial',t
%      'TrialClusterSpikeCount','Trial',t,'Cluster',c
%      'TrialClusterSpikeTime','Trial',t,'Cluster',c,'Spike',s
%
%   Dependencies: None.

Args = struct('SessionName',0,'GroupName',0,'GroupPlotProperties',0,'chunkSize',0, ...
    'showEvents',0,'showPSTH',0);
Args = getOptArgs(varargin,Args,'flags',{'SessionName','GroupName','showEvents','showPSTH'});

if(Args.SessionName)
    r = s.data.sessionname;
elseif(Args.GroupName)
    r = s.data.groupname;
elseif(Args.GroupPlotProperties>0)
    if(Args.GroupPlotProperties>1)
        if Args.showEvents
            r.separate = 'Vertical';
        elseif Args.showPSTH
            r.separate = 'No'
        else
            r.separate = 'No';
        end
    else
        r.separate = 'No';
    end
elseif(Args.chunkSize)
    r=getChunkInfo(s.data.sessionname,Args.chunkSize);
else
    r = get(s.nptdata,varargin{:});
end
