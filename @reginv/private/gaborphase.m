function gv = gaborphase(numPhase,PhaseIndex,SpatialFreq,numPixels)
% input the number of phases for the gabor, the phase to start on, the
% Spatial frequency of the gabor and the number of pixels for the final gv
% square matrix.
phase = PhaseIndex/numPhase*pi*2;
numPix = numPixels/2;
index = 0;
for y = -numPix:numPix-1
	for x = -numPix:numPix-1
		wndw = exp(-0.5*(x*x+y*y)/SpatialFreq);
		grating = cos(2*pi*y/SpatialFreq - phase);
		gabor_values(index+1) = round((wndw*grating + 1)/2 * 255);
		index = index + 1;
	end
end
gv = reshape(gabor_values,numPixels,numPixels);