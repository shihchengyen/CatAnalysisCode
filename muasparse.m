cwd=pwd;
ndcells = nptdata('SessionsFile','longmoviesessions.txt');
clist=get(ndcells,'SessionDirs')';
num_sessions=size(clist,1);
SS=[];
for a=1:num_sessions
    cd(clist{a});
    %WriteGroupMUAFiles;
    %StimInfo('auto','redo','save');
    S=ProcessSession(sparsity,'auto','redo','Rate','ISpikes');
    SS=[SS;S.data.sparsity];
    Values{a}=S.data.Values;
end
Sparse.Sparseness_Values=SS;
Sparse.Firing_Rates=Values;
Sparse.Sigma=3.5;
Sparse.SessionsList=clist;
cd(cwd);
save('Sparse3.5.mat','Sparse')