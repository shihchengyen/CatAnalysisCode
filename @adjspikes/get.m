function r = get(s,varargin)
%adjspikes/GET Returns object properties
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

Args = struct('SessionName',0,'GroupName',0,'CellName',0, ...
	'GroupPlotProperties',0,'chunkSize',0,'SeparateGroupPlot',0);
Args.flags = {'SessionName','GroupName','CellName','SeparateGroupPlot'};
Args = getOptArgs(varargin,Args);

if(Args.SessionName)
   r = s.data.sessionname;
elseif(Args.GroupName)
    r = s.data.groupname;
elseif(Args.CellName)
    r = s.data.cellname;
elseif(Args.GroupPlotProperties>0)
    if(Args.SeparateGroupPlot)
        r.separate = 'Vertical';
    else
        r.separate = 'No';
    end
elseif(Args.chunkSize)
    r=getChunkInfo(s.data.sessionname,Args.chunkSize);
else
    r = get(s.nptdata,varargin{1});
end
