function p = get(obj,prop_name,varargin)
%catcell/get Returns appropriate object properties
%   VALUE = GET(OBJ,PROP_NAME,VARARGIN) returns the object property 
%   specified by PROP_NAME. PROP_NAME can be one of the following:
%      'Number'

Args = struct('Surrogates',0,'SurrogateSets',1000);

Args = getOptArgs(varargin,Args,'flags',{'Surrogates'});

if(Args.Surrogates & strcmp(lower(prop_name),'number'))
	p = Args.SurrogateSets;
else
	% call get function from parent class
	p = get(obj.nptdata,prop_name);
end
