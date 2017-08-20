function [obj, varargout] = plot(obj,varargin)
%REFRACTORY/plot - Plots data from REFRACTORY object
%   OBJ = PLOT(OBJ,N,VARARGIN) plots either the ISI histogram
%   or the recovery function for cell N.
%
%   The optional arguments are:
%      recovery - flag used to plot the recovery function instead
%                 of the ISI histogram.
%      qt - flag used to plot the free firing rate.
%      Wt - flag used to plot the spike probability function.
%      rt - flag used to plot the firing rate.
%      xMax - the data will be plotted from 0 ms to xMax ms
%             (default is 15 ms). This option is valid for both the
%             ISI histogram and the recovery function.
%      legFontSize - font size used for legend (default is 7 pts).
%                    This option is valid for both the ISI histogram 
%                    and the recovery function.
%      yScale - scale for y-axis [ {log} | linear ]. Valid only for
%               ISI histogram.
%      xScale - scale for x-axis [ log | {linear} ]. Valid only for
%               ISI histogram.
%
%   OBJ = PLOT(OBJ,N,VARARGIN)


Args = struct('xMax',15, ...
			  'yScale','log', ...
			  'xScale','linear', ...
			  'legFontSize',7, ...
			  'recovery',0, ...
			  'qt',0, ...
			  'Wt',0, ...
			  'rt',0, ...
              'ShowLegend',0, ...
              'Surrogates',0, ...
              'SurrogateSetsPerFile',100, ...
              'SurrogateBasename','framesg', ...
              'LabelsOff',0, ...
              'ReturnVars', {''}, ...
              'ArgsOnly',0);
Args.flags = {'recovery','qt','Wt','rt','ShowLegend','Surrogates','LabelsOff', ...
    'ArgsOnly'};
Args = getOptArgs(varargin,Args);

if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
    n = Args.NumericArguments{1};
else
    n = 1;
end

if (Args.recovery)
	% plot the fitted recovery function along with data points
	rfs = obj.data.wt{n};
	% use rtEdges instead of edges since function can be computed at 
	% binsize different from ISI binsize
	rfsx = obj.data.rtEdges{n}(1:length(rfs));
	plot(rfsx,rfs,'r-');
	hold on
	plot(obj.data.edges(1:obj.data.maxi(n)),obj.data.rfd{n},'d');
	hold off	
	xlim([0 Args.xMax])
    if(Args.ShowLegend)
    	l = legend(obj.data.cellname,0);
		set(l,'FontSize',Args.legFontSize);
    end
	title(getDataDirs('ShortName','DirString',obj.data.cellname));
    if(~Args.LabelsOff)
		xlabel('Time (ms)');
		ylabel('Probability');
	end
elseif(Args.qt | Args.Wt | Args.rt)
	if (Args.qt)
		% plot free firing rate, q(t), along with W(t) and r(t)
		% add a zero at the end of qt so it will match up with rtEdges and 
		% plot up to obj.data.duration rather than just 1 bin before
		stairs(obj.data.rtEdges{n},[obj.data.qt{n}; 0],'r');
		title(getDataDirs('ShortName','DirString',obj.data.cellname))
		if(~Args.LabelsOff)
			xlabel('Time (ms)');
			ylabel('Free Rate');
		end
		% turn hold on in case we want to plot multiple of these
		hold on
    end
	if (Args.Wt)
		stairs(obj.data.rtEdges{n},[obj.data.Wt{n}; 0],'g');
		title(getDataDirs('ShortName','DirString',obj.data.cellname))
		if(~Args.LabelsOff)
			xlabel('Time (ms)');
			ylabel('Probability');
		end
		% turn hold on in case we want to plot multiple of these
		hold on
    end
	if (Args.rt)
		stairs(obj.data.rtEdges{n},[obj.data.rt{n}; 0],'b');
		title(getDataDirs('ShortName','DirString',obj.data.cellname))
		if(~Args.LabelsOff)
			xlabel('Time (ms)');
			ylabel('Firing Rate');
		end
		% turn hold on in case we want to plot multiple of these
		hold on
    end
	% turn hold off now
	hold off
elseif(Args.Surrogates)
	% get number of plot
	n = Args.NumericArguments{1} - 1;
	% figure which file to load
	filen = floor(n/Args.SurrogateSetsPerFile) + 1;
	% figure which set in the file to plot
	setn = rem(n,Args.SurrogateSetsPerFile) + 1;
	% load surrogates
	sptrain = readSurrogateBin([Args.SurrogateBasename num2str(filen) '.bin']);
    cla
	plotRasters(sptrain{setn})
else
	% plot up to xMax since there are so many points otherwise
	% compute index corresponding to xMax
    xindices = 1:(Args.xMax/obj.data.isibinsize);
	plot(obj.data.edges(xindices),obj.data.hcounts{n}(xindices),'d')
	xlim([0 Args.xMax])
	set(gca,'XScale',Args.xScale,'YScale',Args.yScale);
    if(Args.ShowLegend)
    	l = legend(obj.data.cellname,0);
		set(l,'FontSize',Args.legFontSize);
    end
% 	hold on
% 	pfstart = obj.data.maxi(n);
% 	pfend = pfstart + obj.data.pfpts;
% 	xvals = obj.data.edges(pfstart:pfend);
% 	plot(xvals,exp(obj.data.pf(n,2)).*exp(obj.data.pf(n,1)*xvals),'r-')
% 	hold off
	title(getDataDirs('ShortName','DirString',obj.data.cellname));
	if(~Args.LabelsOff)
		xlabel('Interval (ms)');
		ylabel('Frequency');
	end
end

rvarl = length(Args.ReturnVars);
if(rvarl>0)
     % assign requested variables to varargout
     for rvi = 1:rvarl
         varargout{1}{rvi} = Args.ReturnVars{rvi};
         varargout{1}{rvi*2} = eval(Args.ReturnVars{rvi});
     end
end
