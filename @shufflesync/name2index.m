function n = name2index(obj,name,varargin)
%shufflesync/name2index Method to find event number given input
%   N = name2index(OBJ,NAME) returns the index of the first entry
%   in OBJ.sessiondirs that matches NAME. If NAME is 'end', the last
%   index is returned.
%
%   n = name2index(obj,'name');

% if name is 'end', return the last index
if(strcmp(name,'end'))
	n = get(obj,'Number',varargin{:});
elseif(strncmpi('t:',name,2))
	% example of string will be t: 1655 which indicates we want to find the
	% window that contains the time 1655 ms
	% so first read the time
	wintime = sscanf(name,'t: %d');
	% now find the first window that has a start time larger than wintime
	wi = find(obj.data.WinStartTime>wintime);
	if(~isempty(wi))
		n = wi(1) - 1;
	else
		n = 0;
	end
else
    n = 0;
end
