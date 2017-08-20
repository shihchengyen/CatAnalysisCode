function fixadjxc

adjxc = adjspikesxcorr('auto'); 
if(~isfield(adjxc.data,'lags')) 
    fprintf('Adding lags field to saved adjspikesxcorr object...\n');
    ClusterDirs = getDataDirs('GetClusterDirs'); 
    cwd = pwd; 
    cd(ClusterDirs{1}); 
    adjs1 = adjspikes('auto');
    cd(cwd);
    reps = size(adjs1.data.raster,2); 
    framepoints = adjs1.data.adjFramePoints;
    frames = length(framepoints) - 1;
    repframes = frames/reps;
    fp2 = framepoints(1:(repframes+1)) - framepoints(1);
    [subbins,binSize,nsbins] = divideBins(fp2,'SubBinSize',1);
    l = load('asxcdata');
    adjxc.data.lags = subbins(l.startbins);
    save adjspikesxcorr adjxc
end
