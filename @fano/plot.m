function obj = plot(obj,varargin)
%@fano/plot Plot function for FANO object.
%   OBJ = plot(OBJ,VARARGIN) plots the variance versus the mean of the 
%   spike counts for all frames. The frames that have Fano Factors 
%   lower than 95% of the surrogates are plotted in red. These are the 
%   optional input arguments:
%      Numeric - plots the data for cell N.
%      'Indices' - plots the data for selected frames.
%
%   OBJ = plot(OBJ,'SurrogatePercent',VARARGIN) plots the percentage of
%   surrogates with higher Fano for each frame in cell N. The 
%
%   OBJ = plot(OBJ,'ZScores',VARARGIN) plots the z-scores of the Fano 
%   as compared to the surrogate data for the same frame.
%
%   OBJ = plot(OBJ) plots the variance versus the mean of the spike 
%   counts for each frame in all cells.
%
%   OBJ = plot(OBJ,'Indices',INDICES,VARARGIN) plots the variance 
%   versus the mean of the spike counts for the frames specified in 
%   INDICES. The title of the plot can be passed to the function if
%   it is prefixed by 'TitleString'.
%
%   OBJ = plot(OBJ,'SurrogatePercent',VARARGIN) plots the percentage of
%   surrogates with higher Fano for each frame in all cells.
%
%   OBJ = plot(OBJ,'SurrogatePercent','Indices',INDICES,VARARGIN) plots 
%   the percentage of surrogates with higher Fano for the frames
%   specified in INDICES. The title of the plot can be passed to the 
%   function if it is prefixed by 'TitleString'.
%
%

Args = struct('SurrogatePercent',0,'SurrPercentThreshold',0.95, ...
	'Indices',[],'TitleString','','HistCStep',0.05,'Hist',0,'ZScores',0, ...
	'SurrZScoreThreshold',3);
Args = getOptArgs(varargin,Args,'flags',{'SurrogatePercent','Hist','ZScores'});

if(~isempty(Args.Indices))
	indices = Args.Indices;	
	titlestr = Args.TitleString;
elseif(~isempty(Args.NumericArguments))
	n = Args.NumericArguments{1};
	% get the indices that go with cell n
	indices = find(obj.data.cellid==n);
	titlestr = obj.data.cellname{n};
else
	% default is to plot all data
	indices = 1:length(obj.data.scmean);
	titlestr = ['Data from all ' num2str(length(obj.data.cellname)) ' cells'];
end

if(Args.SurrogatePercent)
	surrPercent = obj.data.surrPercent(indices);
	% find values above threhold
	sigindices = find(surrPercent>Args.SurrPercentThreshold);
	if(Args.Hist)
		histcedges = 0:Args.HistCStep:1;
		n = histcie(surrPercent,histcedges);
		logbar(histcedges,n);
        if(~isempty(sigindices))
			hold on
			n2 = histcie(surrPercent(sigindices),histcedges);
			hb2 = logbar(histcedges,n2,'histc');
			set(hb2,'FaceColor','r')
            hold off
        end
		xlim([-0.01 1.01])
		xlabel('Percent of surrogates with higher Fano Factor')
		ylabel('Number of Frames')
		title(titlestr)
	else
		plot(surrPercent,'.')
		% get x values
		xvals = 1:length(surrPercent);
		hold on
		% plot in red
		plot(xvals(sigindices),surrPercent(sigindices),'r.')
		hold off
		ylabel('Percent of surrogates with higher Fano Factor')
		xlabel('Frame Number')
		title(titlestr)
	end
elseif(Args.ZScores)
	surrZScores = obj.data.surrZScores(indices);
	% find values above threhold
	sigindices = find(-surrZScores>Args.SurrZScoreThreshold);
	if(Args.Hist)
		histcedges = (floor(min(surrZScores)/Args.HistCStep)*Args.HistCStep):Args.HistCStep:(ceil(max(surrZScores)/Args.HistCStep)*Args.HistCStep);
		n = histcie(surrZScores,histcedges);
		logbar(histcedges,n,'histc');
        if(~isempty(sigindices))
			hold on
			n2 = histcie(surrZScores(sigindices),histcedges);
			hb2 = logbar(histcedges,n2,'histc');
			set(hb2,'FaceColor','r')
            hold off
        end
		xlabel('Fano ZScore')
		ylabel('Number of Frames')
		title(titlestr)
	else
		plot(surrZScores,'.')
		% get x values
		xvals = 1:length(surrZScores);
		hold on
		% plot in red
		plot(xvals(sigindices),surrZScores(sigindices),'r.')
		hold off
		ylabel('Fano ZScore')
		xlabel('Frame Number')
		title(titlestr)
	end
else
	% plot the variance versus the mean
	plot(obj.data.scmean(indices),obj.data.scstd(indices).^2, ...
		'Color',[0.5 0.5 0.5],'Marker','.','LineStyle','none');
	hold on
	% find points where the surrPercent is higher than the threshold
	sigindices = intersect(indices,find(obj.data.surrPercent>Args.SurrPercentThreshold));
	plot(obj.data.scmean(sigindices),obj.data.scstd(sigindices).^2,'r.')
	% get axis limits
	ax1 = axis;
	% find the smaller of xmax or ymax
	lmax = min([ax1(2) ax1(4)]);
	% draw y = x line
	line([0 lmax],[0 lmax],'Color','k');
	% draw min variance scallops
	x = 0:0.1:1;
	minvar = x .* (1 - x);
	for i = 0:(ax1(2)-1)
		plot(i+x,minvar,'k')
	end
	hold off
	xlabel('Mean spike count')
	ylabel('Variance')
	title(titlestr)
end
