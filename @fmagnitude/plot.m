function obj = plot(obj,varargin)
%@fmagnitude/plot Plot function for fmagnitude objects.
%   OBJ = plot(OBJ,VARARGIN) plots the FFT amplitude.
%
%   These are the optional input arguments:
%      ShowDC - flag specifying that the DC amplitude should be
%               plotted as well.
%      NormRange - array specifying the min and max frequencies to be 
%                  used to normalize population plots of the mean. The 
%                  normalization values are computed for each column by 
%                  taking the mean of the first value smaller than the 
%                  min to the first value larger than the max. If 
%                  NormRange is empty (i.e. []), no normalization is 
%                  applied. The default is a NormRange of [5 6].
%
%   obj = plot(obj,'ShowDC','NormRange',[5 6]);

Args = struct('ShowDC',0,'XLim',[],'NormRange',[5 6],'TimeFrequency',0);
Args = getOptArgs(varargin,Args,'flags',{'ShowDC','TimeFrequency'});

% get length of frequencies
fl = length(obj.data.f);
if(Args.ShowDC)
	plotpts = 1:fl;
elseif(~isempty(Args.XLim))
    % make sure XLim has 2 entries
    if(length(Args.XLim)~=2)
        error('Error: XLim should have 2 numbers!');
    else
        % find first value in obj.data.f that exceeds XLim(1)
        fimin = find(obj.data.f>Args.XLim(1));
        % find first value in obj.data.f that exceed XLim(2)
        fimax = find(obj.data.f>Args.XLim(2));
        if(isempty(fimax))
            plotpts = fimin(1):fl;
        else
            plotpts = fimin(1):(fimax(1)-1);
        end
    end
else
	% plot without DC
	plotpts = 2:fl;
end
x = obj.data.f(plotpts);

if(isempty(Args.NumericArguments))
	% plot population data
	if(~isempty(Args.NormRange))
		% plot all means by normalizing by values in NormRange
		% find indexes in f that correspond to NormRange
		% find first index in f that is smaller than NormRange(1)
		fimin = find(obj.data.f>Args.NormRange(1));
		normfmin = fimin(1) - 1;
		% find first index in f that is larger than NormRange(2)
		fimax = find(obj.data.f>Args.NormRange(2));
		normfmax = fimax(1);
		% get the mean values in the range of normfmin to normfmax
		meanvalues = mean(obj.data.mean(normfmin:normfmax,:));
		% divide obj.data.mean by meanvalues
		y = obj.data.mean(plotpts,:) ./ repmat(meanvalues,size(plotpts,2),1);
	else
		y = obj.data.mean(plotpts,:);
	end
	plot(x,y)
    title(['FFT Magnitude ' num2str(obj.data.nCells) ' Cells']);
    xlabel('Frequency (Hz)')
    ylabel('FFT Magnitude')
	if(~isempty(Args.XLim))
		xlim(Args.XLim);
	end
else
	% plot data from individual cell
	n = Args.NumericArguments{1};
	if(Args.TimeFrequency)
		% get column indices for the fftMag matrix
		cstart = obj.data.colIndex(n) + 1;
		cend = obj.data.colIndex(n+1);
		imagesc(obj.data.timeVec,x,obj.data.fftMag(plotpts,cstart:cend));
		if(~isempty(Args.XLim))
			ylim([Args.XLim])
		end
        xlabel('Time (ms)')
        ylabel('Frequency (Hz)')
	else
		y = obj.data.mean(plotpts,n);
		y1 = obj.data.stdev(plotpts,n);
		plot(x,y)
		hold on
		plot(x,y+y1,'g')
		plot(x,y-y1,'g')
		hold off
		if(~isempty(Args.XLim))
			xlim(Args.XLim);
		end
        xlabel('Frequency (Hz)')
        ylabel('FFT Magnitude')
	end
	title(obj.data.cellname{n});
end

zoom xon
