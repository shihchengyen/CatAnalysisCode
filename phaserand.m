% phaserand.m - randomizes phase spectrum of an image
%
% function imr = phaserand(im,fig)
%
% im: input image
% fig: flag for displaying result

function imr = phaserand(im,fig)

[n m]=size(im);

% windowing function
[x y]=meshgrid([-m/2:m/2-1]/m,[-n/2:n/2-1]/n);
sigma=0.25;
W=exp(-0.5*(x.^2+y.^2)/sigma^2);

% get the 2D fft of the image
% note: im has to be of type double in order for the operation .* to work
imF = fft2(W.*im);

% get amplitude and phase
a = abs(imF);
phi = angle(imF);

% reshapes phi into a row vector
nm = n * m;
phi1 = reshape(phi,nm,1);
% scrambles the phase using randperm function
rphi = phi1(randperm(nm));
% reshape rphi back into original dimensions
rand_phi = reshape(rphi,n,m);
rand_phi(1)=0; % ensures positive DC component

% compute phase randomized image
imr=real(ifft2(a.*exp(j*rand_phi)));
imr=imr-mean(imr(:))+mean(im(:));

if exist('fig')
   
    minp=min(im(:));
    maxp=max(im(:));

    figure(fig)
    colormap gray
    subplot(121)
    imagesc(im,[0 255]), axis image
    title('original')
    subplot(122)
    imagesc(imr,[0 255]), axis image
    title('phase randomized')
    drawnow
    
    figure(fig+1)
    subplot(121)
    hist(im(:));
    xlabel('pixel value')
    title('original histogram')
    subplot(122)
    hist(imr(:));
    xlabel('pixel value')
    title('phase randomized histogram')
    
end
