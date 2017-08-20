function n = name2index(obj,name)
%trentropy/name2index Method to find session index given session name
%   N = name2index(OBJ,NAME) returns the index of the first entry
%   in OBJ.sessionname that matches NAME. If NAME is 'end', the last
%   index is returned.
%
%   n = name2index(obj,'name');

% if name is 'end', return the last index
if(strcmp(name,'end'))
	n = obj.data.cellid(end);
else
	% find entries in cell array that match name
	goodcells = regexp(obj.data.cellname,name);
	% find entries that had a match
	goodlist = find(~cellfun('isempty',goodcells));
	if(~isempty(goodlist))
		% if there are multiple matches, return the first one
		n = goodlist(1);
	else
		% if no matches, return 0
		n = 0;
	end
end
