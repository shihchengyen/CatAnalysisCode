function r = get(s,varargin)
%latency/GET Returns object properties
%   VALUE = GET(OBJ,PROP_NAME) returns an object 
%   property. PROP_NAME can be one of the following:

Args = struct('GroupPlotProperties',0);
Args.flags = {''};
Args = getOptArgs(varargin,Args);

if(Args.GroupPlotProperties>0)
	r.separate = 'Horizontal';
else
    r = get(s.nptdata,varargin{1});
end
