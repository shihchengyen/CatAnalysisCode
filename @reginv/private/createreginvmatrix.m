function [reginvmatrix] = createreginvmatrix(numPixels)
% Creates the Laplacian ID Matrix based on the numPixels^2.
%
% 2D Mexican Hat filter
% Ls = [0 -1 0
%      -1 4 -1
%       0 -1 0];

% Edges have issues to deal with
L_Edge = 2:1:numPixels-1;
R_Edge = (numPixels^2)-(numPixels-2):1:(numPixels^2)-1;
T_Edge = numPixels+1:numPixels:((numPixels^2)-(numPixels*2))+1;
B_Edge = numPixels*2:numPixels:numPixels*(numPixels-1);

% Load into Memory the Laplacian Square Matrix
reginvmatrix = zeros(numPixels^2);

% Loop to convolve the MexHat with each point on the Zero_Mat and create
% the reginvmatrix
for ii = 1:(numPixels^2)
    Z_Mat = zeros(numPixels);
    if ii == 1
        Z_Mat(ii)=4; Z_Mat(ii+numPixels)=-1; Z_Mat(ii+1)=-1;
    elseif ii == numPixels
        Z_Mat(ii)=4; Z_Mat(ii+numPixels)=-1; Z_Mat(ii-1)=-1;
    elseif ii == (numPixels^2)-(numPixels-1)
        Z_Mat(ii)=4; Z_Mat(ii+1)=-1; Z_Mat(ii-numPixels)=-1;
    elseif ii == (numPixels^2)
        Z_Mat(ii)=4; Z_Mat(ii-1)=-1; Z_Mat(ii-numPixels)=-1;
    elseif ismember(ii,L_Edge)
        Z_Mat(ii)=4; Z_Mat(ii-1)=-1; Z_Mat(ii+1)=-1; Z_Mat(ii+numPixels)=-1;
    elseif ismember(ii,R_Edge)
        Z_Mat(ii)=4; Z_Mat(ii-1)=-1; Z_Mat(ii+1)=-1; Z_Mat(ii-numPixels)=-1;
    elseif ismember(ii,T_Edge)
        Z_Mat(ii)=4; Z_Mat(ii-numPixels)=-1; Z_Mat(ii+1)=-1; Z_Mat(ii+numPixels)=-1;
    elseif ismember(ii,B_Edge)
        Z_Mat(ii)=4; Z_Mat(ii-numPixels)=-1; Z_Mat(ii-1)=-1; Z_Mat(ii+numPixels)=-1;
    else
        Z_Mat(ii)=4; Z_Mat(ii-numPixels)=-1; Z_Mat(ii+numPixels)=-1; Z_Mat(ii-1)=-1; Z_Mat(ii+1)=-1; 
    end  
    reginvmatrix(ii,:) = reshape(Z_Mat,1,numPixels^2);   
end