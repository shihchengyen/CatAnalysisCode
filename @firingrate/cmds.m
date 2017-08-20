[ndcells,frcells] = ProcessDirs(ndcells,'Object','firingrate');
gpeak = heterogeneity(frcells,'NoSingles');
% sort the correlation coefficients within each site so the
% plot is more readable
% save the indices so we can use it later to mark the
% insignificant values
[gsdata,gsdi] = sort(gpeak,'descend');            
gdata = flipud(gsdata);
gsdi2 = flipud(gsdi);
[gdrows,gdcols] = size(gdata);
% add 2 rows of nan's so we can easily put a little bit of
% space between neighboring sites
gdata = [gdata; repmat(nan,2,gdcols)];
gdrows2 = gdrows + 2;
gdsize = [gdrows2,gdcols];
% find where the values are in gdata
gidx = ~isnan(gdata);
% count the number of pairs at each site
gdsum = sum(gidx);
% add 2 to account for the space between sites
gdnum2 = gdsum + 2;
gdc = 1:gdcols;
gdmax = nanmax(gdata(:));
fakeval = 2*gdmax;
% put a fake value in gdata so it will be extracted and thus
% create the spacing between sites
% use 1.1 since correlation coefficients cannot exceed 1
gdata(sub2ind(gdsize,gdsum+1,gdc)) = fakeval;
gdata(sub2ind(gdsize,gdnum2,gdc)) = fakeval;
% now find all non-nan values, which will include the 2 fake
% values added to each site
gidx2 = ~isnan(gdata);
% create a similar matrix to store the x values so we can plot
% neighboring sites with slightly different colors
gdx = gdata;
% put in the x values that will be used in the bar plot
gdx(gidx2) = 1:sum(gdnum2);
% plot the odd and even columns separately so we can use
% different colors
ocols = 1:2:gdcols;
ecols = 2:2:gdcols;
% get the odd columns first
% get the non-nan values
gidx2o = gidx2(:,ocols);
% get the data (e.g. correlation coefficients)
gdata1o = gdata(:,ocols);
% pull out the values into a column vector
gdata1 = gdata1o(gidx2o);
% convert the fake values back to nan
gdata1(gdata1==fakeval) = nan;
% get the x values for plotting
gdxo = gdx(:,ocols);
% pull out the x values into a column vector
gdx1 = gdxo(gidx2o);
% now create the bar plot
bar(gdx1,gdata1,1,'w');
hold on
% now do the same for even columns
gidx2e = gidx2(:,ecols);
gdata1e = gdata(:,ecols);
gdata1 = gdata1e(gidx2e);
gdata1(gdata1==fakeval) = nan;
gdxe = gdx(:,ecols);
gdx1 = gdxe(gidx2e);
bar(gdx1,gdata1,1,'FaceColor',repmat(1,1,3));
