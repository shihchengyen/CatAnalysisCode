function gdf = adjs2gdf

% add time shift to allow the unitary events analysis to work
% time shift is really 5000 ms but since we are using the raster
% field, we need to add the repetition time as well
tshift = 5000;
rmarker = 777;
fmarker = 888;

a = adjspikes('auto');
% get raster
raster = a.data.raster;
% get adjFramePoints
adjfp = a.data.adjFramePoints;
% arows will be the maximum number of spikes in a repetition
% acols will be the number of repetitions
[arows,acols] = size(raster);
% first point of adjFramePoints is not always 0 so we subtract by the first
% point to make sure everything is referenced to 0
lfpts = length(a.data.adjFramePoints) - 1;
% compute number of rows since acols is known
frows = lfpts / acols;
% reshape adjFramePoints into matrix so we can easily pull out the start of
% each repetition
fpts = reshape(a.data.adjFramePoints(1:lfpts),frows,acols) - adjfp(1);
% create time shift vector
tsvec = (0:(acols-1))*tshift;
% replicate it to crate a matrix of the same dimensions
tsmat = repmat(tsvec,frows,1);
% add two matrices together
sfpts = fpts + tsmat;
shiftedfpts = sfpts(:);
% grab start of each repetition with time shifted by tshift
rs = sfpts(1,:);
% replicate it and add to rasters
sr = raster + repmat(rs,arows,1);
repstart = rs';
% return only the non-nan entries in stimes
sr2 = sr(:);
sr2idx = isnan(sr2);
stimes = sr2(~sr2idx);
[g,c] = nptFileParts(pwd); 
gc = [g(end) c(end-2:end-1)]; 
gdf = [[repmat(str2num(gc),length(stimes),1) stimes]; [repmat(rmarker,acols,1) repstart]; [repmat(fmarker,lfpts,1) shiftedfpts]];
