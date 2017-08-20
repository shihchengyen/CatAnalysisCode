function job4

cdirs = getDataDirs('GetClusterDirs'); 
cwd = pwd; 
cd(cdirs{1}); 
adjs1 = adjspikes('auto'); 
reps = size(adjs1.data.raster,2); 
framepoints = adjs1.data.adjFramePoints; 
frames = length(framepoints) - 1; 
repframes = frames/reps; 
firstframe = repframes + 1; 
totalframes = repframes * 2; 
system(['asxcwrapper2 ' num2str(firstframe) ' ' num2str(totalframes) ' 100 asxcdata1.mat asxcsurrA &']);
