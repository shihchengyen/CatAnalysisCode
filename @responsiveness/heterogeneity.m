function H = heterogeneity(obj,varargin)
%responsiveness/heterogeneity calculates the heterogeneity between group values
%   G = groupDirs(OBJ) returns a matrix containing indices 
%   corresponding to different grouping of the session directories.
%   corresponding to each group in each column, padded with NaN's.

Args = struct('Pairs',0,'STD',0,'Diff',0,'Mean',0,'Bar',0);
Args = getOptArgs(varargin,Args,'flags',{'Pairs','STD','Diff','Mean','Bar'});

if(Args.Pairs)
	sgi = groupDirs(obj,'Pairs');
else
	sgi = groupDirs(obj);
end

if Args.STD
    % get entries that are NaN's
    sgnani = find(isnan(sgi));
    % replace nan with 1
    sgi(sgnani) = 1;
    % get valuess
    r = obj.data.responsiveness;
    % use index to get values
    avals = r(sgi);
    % put NaN's back
    avals(sgnani) = nan;
    % compute standard deviation
    H = nanstd(avals);
end

if Args.Diff
    % get valuess
    r = obj.data.responsiveness;
    % use index to get values
    avals = r(sgi);
    % compute the absolute differences 
    H = abs(diff(avals));
    if Args.Mean
        Mean_H=[];
        cellprefix = getDataDirs('CellPrefix');
        s = get(obj,'SessionDirs');
        % remove cluster01s etc so we get pathnames up to the group
        s1 = regexprep(s,[filesep cellprefix '.*'],'');
        % get the unique group directories and the corresponding indices
        [s2,s2i,s2j]  = unique(s1);
        s_index = sgi(1,1:end);
        for mi = 1:length(s2i)
            if mi == 1
                index = find(s_index <= s2i(mi));
            else
                index = find(s_index > s2i(mi-1) & s_index <= s2i(mi));
            end
            if ~isempty(index)
                Mean_H = [Mean_H mean(H(index))];
            end
        end      
        H = Mean_H;
    end    
end

if Args.Bar
    % get entries that are NaN's
    sgnani = find(isnan(sgi));
    % replace nan with 1
    sgi(sgnani) = 1;
    % get valuess
    r = obj.data.responsiveness;
    % use index to get values
    avals = r(sgi);
    % put NaN's back
    avals(sgnani) = nan;
    H = avals;
end
