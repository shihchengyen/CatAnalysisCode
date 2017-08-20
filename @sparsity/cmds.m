% compute frame level sparseness 
[ndcells,sparseframe] = ProcessDirs(ndcells,'Object','sparsity','redo','Repetitions');
[ndresp,sparseframerate] = ProcessDirs(ndresp,'Object','sparsity','Repetitions','Rate');
% no different between using counts versus rate so we will just use spike
% counts
% compute sparseness using 3 different window sizes
[ndresp,sparse10] = ProcessDirs(ndresp,'Object','sparsity','redo','IgnoreFrames','WindowSize',10);
[ndresp,sparse50] = ProcessDirs(ndresp,'Object','sparsity','redo','IgnoreFrames','WindowSize',50);
[ndresp,sparse100] = ProcessDirs(ndresp,'Object','sparsity','redo','IgnoreFrames','WindowSize',100);
% combine data to do boxplot and compute significance
sparsedata = [sparse10.data.sparsity sparseframe.data.sparsity ...
        sparse50.data.sparsity sparse100.data.sparsity]; 
figure
boxplot(sparsedata,1)
% do lillietest to see if distributions are gaussian
[r,p] = lillietest(sparseframe.data.sparsity)
[r,p] = lillietest(sparse10.data.sparsity)
[r,p] = lillietest(sparse50.data.sparsity)
[r,p] = lillietest(sparse100.data.sparsity)
% all returned r = 1, p = NaN
% do non-parametric test
p = kruskalwallis(sparsedata)
% p = 0.11
% compute frame level sparseness 
[ndcells,sparseframe] = ProcessDirs(ndcells,'Object','sparsity','redo','Repetitions');
[ndcells,sparseframerate] = ProcessDirs(ndcells,'Object','sparsity','Repetitions','Rate');
% compute sparseness using 3 different window sizes
[ndcells,sparse10] = ProcessDirs(ndcells,'Object','sparsity','redo','IgnoreFrames','WindowSize',10);
[ndcells,sparse50] = ProcessDirs(ndcells,'Object','sparsity','redo','IgnoreFrames','WindowSize',50);
[ndcells,sparse100] = ProcessDirs(ndcells,'Object','sparsity','redo','IgnoreFrames','WindowSize',100);
% combine data to do boxplot and compute significance
sparsedata = [sparse10.data.sparsity sparseframe.data.sparsity ...
        sparse50.data.sparsity sparse100.data.sparsity]; 
figure
boxplot(sparsedata,1)
% do lillietest to see if distributions are gaussian
[r,p] = lillietest(sparseframe.data.sparsity)
[r,p] = lillietest(sparse10.data.sparsity)
[r,p] = lillietest(sparse50.data.sparsity)
[r,p] = lillietest(sparse100.data.sparsity)
% all returned r = 1, p = NaN
% do non-parametric test
p = kruskalwallis(sparsedata)
% p = 0.0356
% significant but only due to difference between 10 and 100 ms
p = kruskalwallis(sparsedata(:,1:3))
% p = 0.1204
p = kruskalwallis(sparsedata(:,2:4))
% p = 0.3446

% plot sparseness heterogeneity
plot(sparseframe,'Heterogeneity','IntraGroup','Bar');