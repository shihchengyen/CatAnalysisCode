function H = heterogeneity(obj,varargin)
%firingrate/heterogeneity calculates the heterogeneity between group values

Args = struct('Mean',0);
Args.flags = {'Mean'};
Args = getOptArgs(varargin,Args);

sgi = groupDirs(obj,varargin{:});

if(Args.Mean)
	% get mean for each column
	data = nanmean(obj.data.firingRate);
else
	% get max for each column
	data = max(obj.data.firingRate,[],1);
end

H = nanindex(data,sgi);
