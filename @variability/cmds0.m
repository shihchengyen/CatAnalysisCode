[ndresp,sgisi] = ProcessDirs(ndresp,'nptDirCmd','a1 = readSurrogateBin(''framesg1.bin''); sgn = length(a1); repn = length(a1{1}); repsd = zeros(repn,1); for repi = 1:repn, repd = zeros(sgn,1); for sgi = 1:sgn, repd(sgi) = ~isempty(find(diff(a1{sgi}{repi})==0)); end; repsd(repi) = sum(repd); end; isireps = find(repsd); data = [data; [repmat(i,length(isireps),1) isireps repsd(isireps)]];');

% explanation of above command
% read surrogate file
a1 = readSurrogateBin('framesg1.bin'); 
% get number of surrogates
sgn = length(a1); 
% get number of repetitions
repn = length(a1{1}); 
% allocate memory
repsd = zeros(repn,1); 
% look over repetitions
for repi = 1:repn, 
    % allocate memory
    repd = zeros(sgn,1); 
    % loop over surrogates
    for sgi = 1:sgn, 
        % find inter-spike-intervals that are 0
        repd(sgi) = ~isempty(find(diff(a1{sgi}{repi})==0)); 
    end; % end loop over number of surrogates 
    % find number of intervals that are 0
    repsd(repi) = sum(repd); 
end; % end loop over number of repetitions
% find repetitions that have isi's that are 0
isireps = find(repsd); 
% return cell number, repetition number, and number of isi's that are 0
data = [data; [repmat(i,length(isireps),1) isireps repsd(isireps)]];
% find unique cell numbers
us = unique(sgisi(:,1));
% get session directories
s = ndresp.SessionDirs';
% find unique session directories
usd = s(us);
% remove path prefix
usd2 = strrep(usd,'/Users/syen/Documents/ShihCheng/Data/Neural/Cat/newcatdata/','');
