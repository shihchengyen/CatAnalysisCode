function obj = plot(obj,varargin)
%@trentropy/plot Plot function for TRENTROPY object.
%   OBJ = plot(OBJ,N) plots the percent of surrogates that have higher
%   entropies than each frame of the data.
%
%   OBJ = plot(OBJ,N) plots the 

Args = struct('SurrPercentThreshold',0.95,'Indices',[],'TitleString','', ...
	'HistCStep',0.05,'Hist',0,'ZScores',0,'SurrZScoreThreshold',3);
Args = getOptArgs(varargin,Args,'flags',{'Hist','ZScores'});

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
	indices = 1:length(obj.data.entropy);
	titlestr = ['Data from all ' num2str(length(obj.data.cellname)) ' cells'];
end

if(Args.ZScores)
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
		xlabel('Entropy ZScore')
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
		ylabel('Entropy ZScore')
		xlabel('Frame Number')
		title(titlestr)
	end
else
	surrPercent = obj.data.surrPercent(indices);
	% find values above threhold
	sigindices = find(surrPercent>Args.SurrPercentThreshold);
	if(Args.Hist)
		histcedges = 0:Args.HistCStep:1;
		[n,x] = histcie(surrPercent,histcedges);
		logbar(histcedges,n,'histc');
        if(~isempty(sigindices))
			hold on
			n2 = histcie(surrPercent(sigindices),histcedges);
			hb2 = logbar(histcedges,n2,'histc');
			set(hb2,'FaceColor','r')
            hold off
        end
		xlim([-0.01 1.01])
		xlabel('Percent of surrogates with higher Entropies')
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
		ylabel('Percent of surrogates with higher Entropies')
		xlabel('Frame Number')
		title(titlestr)
	end
end