function r = getRefreshReversals(exp,nreversals,varargin)
%getRefreshReversals Returns number of refreshes between direction reversals
%   R = getRefreshReversals(EXP,NREVERSALS,VARARGIN) returns a vector
%   of length NREVERSALS containing number of refreshes between
%   direction reversals. The frequency, f, of the reversals are 
%   distributed as 1/(f^EXP).
%
%   The optional input arguments are:
%      'showplot' - plots data for visualization.
%      'refreshrate' - refresh rate of the monitor in Hz (default 
%                      is 150).
%      'filename' - name of file to which R is saved.
%
%   r = getRefreshReversals(exp,nreversals,'showplot','refreshrate',150,'filename','')

Args = struct('showplot',0, ...
			  'refreshrate',150, ...
			  'filename','');
			  
Args = getOptArgs(varargin,Args,'flags',{'showplot'});

% add 1 to exp since the equation for the beta distribution is x^(exp-1)
br = betarnd(exp+1,1,nreversals,1);

% use ceil to round to integer number of refreshes so minimum number of 
% refreshes will be 1
r = ceil(br*Args.refreshrate);

if(~isempty(Args.filename))
	fprintf(fopen(Args.filename,'wt'),'%d\n',r);
end

if(Args.showplot)
	[hbr,nbr] = hist(r);
	xvals = Args.refreshrate./nbr;
	h1 = plot(xvals,hbr,'o-');
	hold on
	h2 = plot(xvals,hbr(end)*1./(xvals.^exp),'r-');
	legend([h1 h2],{'data' ['1/f\^' num2str(exp)]});
	hold off
	t1 = ['Total duration: ' num2str(sum(r)/Args.refreshrate/60,'%.1f') ' mins'];
	t2 = sprintf(' Refresh min: %d max: %d',min(r),max(r));
	title([t1 t2]);
	xlabel('Reversal Frequency (Hz)')
	ylabel('Occurrences')
end
