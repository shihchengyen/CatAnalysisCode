function robj = ProcessCell(obj,varargin)

redo = 0;

if(~checkMarkers(obj,redo,'cell'))
	robj = fano('auto',varargin{:});
	
	createProcessedMarker(obj,'cell');
else
	robj = fano;
end
