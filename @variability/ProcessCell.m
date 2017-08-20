function robj = ProcessCell(obj,varargin)

redo = 0;

if(~checkMarkers(obj,redo,'cell'))
	robj = variability('auto',varargin{:});
	
	createProcessedMarker(obj,'cell');
else
	robj = variability;
end
