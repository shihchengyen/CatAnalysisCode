function combineSurrSynchrony(repframes,surrSets)

surrdatafile = 'asxcsurr';

% needed to compile standalone executable
if(ischar(repframes))
	repframes = str2num(repframes);
end
if(ischar(surrSets))
	surrSets = str2num(surrSets);
end

% allocate memory
xc = zeros(surrSets,repframes);
for idx = 1:repframes
	l = load([surrdatafile num2str(idx,'%04d')]);
	xc(:,idx) = l.xc;
end
save(surrdatafile,'xc');
