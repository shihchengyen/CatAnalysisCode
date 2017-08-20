function robj = ProcessCell(obj,varargin)

redo = 0;

if(~checkMarkers(obj,redo,'cell'))
	robj = trentropy('auto',varargin{:});
	
	createProcessedMarker(obj,'cell');
else
	robj = trentropy;
end
