%Loads the waveforms files and creates the groups and multiunit clusters
%for the session without using mclust

% Session Directory
sdir=pwd;

% Look for the FD directory and move to the sort directory
d=nptDir('sort');
if isempty(d);
    mkdir('sort');
end
d=nptDir('FD');

if ~isempty(d);
    movefile('FD','sort');
else
    mkdir('FD')
end

% Groups of MUA
glist=nptDir('*.hdr');
wlist=nptDir('*waveforms.bin');
numGroups=size(glist,1);

% Write the ispikes files for the MUA
for ai = 1:numGroups
    
    sp=ispikes('Group',glist(ai).name(1:4),'auto','UseSort',0);    
    mkdir(['group' glist(ai).name(1:4) '\cluster01m']);
    filename = [sdir '\group' glist(ai).name(1:4) '\cluster01m\ispikes.mat'];
    save(filename,'sp');
    movefile(glist(ai).name,'sort');
    movefile(wlist(ai).name,'sort');
    
end

movefile('FD','sort')
