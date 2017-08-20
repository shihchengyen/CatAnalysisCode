function [obj, varargout] = plot(obj,varargin)


Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Fano', 0);
Args.flags = {'LabelsOff','ArgsOnly','Fano'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly

    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
else
	% plot all data
	n = 1;
end
% cd /opt/data/cat;
% load ndobjs;
% cd (ndcc.SessionDirs{n});
cd(obj.nptdata.SessionDirs{n})
load ('adjspikes.mat');
xvals = as.data.adjFramePoints(1: length(obj.data.index(n):(obj.data.index(n+1)-1)));
if(xvals(1) ~= 0)
    xvals = xvals - xvals(1);
end
xvals = xvals + (xvals(2) - xvals(1))/2;
if(Args.Fano)
    fanopre = obj.data.fanopre(obj.data.index(n):(obj.data.index(n+1)-1));
    bar(xvals,fanopre,1);
    title('window precision with Fano factor');
    ylim([0 10]);
    xlabel('Time(ms)');
    ylabel('Precision Level');
else
    trepre = obj.data.trepre(obj.data.index(n):(obj.data.index(n+1)-1));
    bar(xvals,trepre,1);
    title('window precision with TRE score');
    ylim([0 10]);
    xlabel('Time(ms)');
    ylabel('Precision Level');
end

rvarl = length(Args.ReturnVars);
if(rvarl>0)
    % assign requested variables to varargout
    for rvi = 1:rvarl
     	 rvidx = rvi * 2;
         varargout{1}{rvidx-1} = Args.ReturnVars{rvi};
         varargout{1}{rvidx} = eval(Args.ReturnVars{rvi});
     end
end


zoom on
