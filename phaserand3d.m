function imr = phaserand3d(im,fig)
% phaserand3d.m - randomizes phase spectrum of a movie
%
% function imr = phaserand3d(im,fig)
%
% im: input movie
% fig: flag for displaying result

[n m o]=size(im);

% windowing function - note that n is rows which is y and m is columns
% which is x
[x y z]=meshgrid([-m/2:m/2-1]/m,[-n/2:n/2-1]/n,[-o/2:o/2-1]/o);
sigma=0.25;
W=exp(-0.5*(x.^2+y.^2+z.^2)/sigma^2);

% get the 2D fft of the image
% note: im has to be of type double in order for the operation .* to work
imF = fftn(W.*im);

% get amplitude and phase
a = abs(imF);
phi = angle(imF);

% reshapes all elements of phi except the DC component into a row vector
nmo = n * m * o - 1;
phi1 = reshape(phi(2:end),nmo,1);
% scrambles the phase using randperm function
rphi = phi1(randperm(nmo));
% add DC component and reshape rphi back into original dimensions
rand_phi = reshape([phi(1,1,1); rphi],n,m,o);
% rand_phi(1,1,1)=0; % ensures positive DC component

% compute phase randomized image
imr=real(ifftn(a.*exp(j*rand_phi)));
% subtract overall mean from output and add overall mean of input
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
