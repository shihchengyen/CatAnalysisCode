function r = getDirectionReversals(a,b,p,duration,varargin)
%getDirectionReversals Get Buracas-like direction reversals
%   R = getDirectionReversals(A,B,P,DURATION,VARARGIN) returns a list of
%   number of refreshes between direction changes of a grating ala 
%   Buracas et al. Direction reversals from the anti-preferred to the
%   preferred direction are followed by A refreshes in the preferred
%   direction, which are followed by B refreshes in the anti-preferred
%   direction, after which the direction switches to the preferred
%   direction with probability P at each refresh. DURATION specifies 
%   the duration of the stimulus in seconds that the algorithm will
%   attempt to approximate.
%
%   These are the optional input arguments:
%      'FrameRate' - frame rate of the stimulus in Hz (default is 25).
%      'FileName' - name of file to use to output the sequence of 
%                   refreshes (default is '');
%      'ShowPlot' - plots the phase sequence and amplitude spectrum.
%      'ShowGrating' - animates a grating using the phase sequence.
%      'ShowAll' - equivalent to using both 'ShowPlot' and 'ShowGrating'.
%      'GratingSize' - defines the size of the grating in pixels
%                      (default is 100).
%      'GratingSF' - defines the number of pixels that make up one
%                    cycle of the grating (default is 10).
%      'GratingTF' - defines the desired temporal frequency of the
%                    grating in cycles/second (default is 1).
%      'FramePause' - defines how long to pause between animating 
%                     frames of the grating in seconds (default is 0.1).
%
%   r = getDirectionReversals(a,b,p,duration,'framerate',25,'filename,'',...
%           'showplot','showgrating','gratingsize',100,'gratingsf',10,...
%           'gratingtf',1,'framepause',0.1);

Args = struct('FrameRate',25,...
			  'FileName','',...
			  'ShowPlot',0, ...
			  'ShowGrating',0, ...
			  'GratingSize',100, ...
			  'GratingSF',10, ...
			  'GratingTF',1, ...
			  'FramePause',0.1);

Args = getOptArgs(varargin,Args,'flags',{'ShowPlot','ShowGrating'},...
					'aliases',{'ShowAll',{'ShowPlot','ShowGrating'}});

% get total number of frames in duration
nrefs = duration*Args.FrameRate;

% get random numbers
probs = rand(1,nrefs);

% find values smaller than p
pind = find(probs<p);

% find frames between pind
pframes = [pind(1) diff(pind)];
pfl = size(pframes,2);

% create vector of a frames
avec = repmat(a,1,pfl);

% list of refreshes can be returned by reshaping pframes and avec into
% a column vector since reshape grabs values columnwise
seq = reshape([avec; (pframes-1+b)],[],1);

% find cummulative sum to trim sequence to duration
cseq = cumsum(seq);
% find first index that exceeds specified duration
csi = find(cseq>nrefs);
% return vector before csi
rl = csi(1) - 1;
r = seq(1:rl);

if(~isempty(Args.FileName))
	fid = fopen(Args.FileName,'wt');
	fprintf(fid,'# generated %s by %s using:\n',date,mfilename);
	fprintf(fid,'# a: %d, b: %d, p: %f, duration: %f seconds, frame rate: %d Hz\n',...
					a,b,p,duration,Args.FrameRate);
	fprintf(fid,'# first number corresponds to number of refreshes where the grating is drifting in the preferred direction\n');
	fprintf(fid,'%d\n',r);
	fclose(fid);
end

if(Args.ShowPlot)
	% generate vectors of 1's and 0's that are half of rl and then
	% reshape to get a vector of alternating 1's and 0's
	rl2 = ceil(rl/2);
	vec10 = [ones(1,rl2); zeros(1,rl2)];
	vec = reshape(vec10,[],1);
	% if rl is even add 1 at the end
	if(rem(rl,2)==0)
		% add another 1 at the end of vec
		vec = [vec; 1];
	end
	% plot sequence
	stairs([1; cseq(1:rl)],vec);
	ylim([-0.05 1.05]);
	xlabel('Frame Number');
	set(gca,'YTick',0:1,'YTickLabel',['Anti-pref';'Preferred']);
	title(['a = ' num2str(a) ', b = ' num2str(b) ', p = ' num2str(p,'%.2f')])
end

if(Args.ShowGrating)
	% create grating
	figure
	a = 1:Args.GratingSize;
	a1 = repmat(a,Args.GratingSize,1);
	% wrap around spatial frequency, add 1 and shift since we want to
	% start from 0
	a2 = circshift(rem(a1,Args.GratingSF) + 1,[0 1]);
	
	% create colormap
	% theta corresponds to spatial phase
	twoPi = 2*pi;
	tstep = twoPi/Args.GratingSF;
	theta = (0:tstep:(twoPi-tstep))';
	% get number of phases
	nphases = round(Args.FrameRate/Args.GratingTF);
	% phi corresponds to temporal phase
	pstep = twoPi/nphases;
	phi = 0:pstep:(twoPi-pstep);
	% create matrices for matrix multiplication
	m1 = [theta ones(Args.GratingSF,1)];
	m2 = [ones(1,nphases); phi];
	m = m1 * m2;
	% compute grayscale in cycle and set range to [0 1]
	cmap = (cos(m)+1)*0.5;
	% repmat so we set rgb values to all be the same
	colormaps = repmat(cmap,[1 1 3]);
	% rearrange so that rgb values are more easily accessible
	cmaps = permute(colormaps,[1 3 2]);
	% set up palette sequence
	cseq = (1:nphases)';
	
	% show grating
	h = imshow(a2,cmaps(:,:,1),'notruesize');
	set(h,'EraseMode','xor');
	set(gcf,'DoubleBuffer','on');
	
	% let user hit return to start animation
	display('Hit return to start animation');
	pause
	
	% start animation
	palettedir = -1;
	for ci = 1:length(r)
		for ri = 1:r(ci)
			cseq = circshift(cseq,palettedir);
			colormap(cmaps(:,:,cseq(1)));
			pause(Args.FramePause);
		end
		palettedir = -palettedir;
	end
end
