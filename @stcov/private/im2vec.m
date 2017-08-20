function vecs=im2vec(im_frames)
% function vecs=im2vec(im_frames)
% im_frames: m x n x N array of N image frames, each m x n in size
% vecs: structure
% vecs.v: N x (m times n) matrix, where i-th row is a vecorized form of the
%         i-th image
% vecs.size: [m n N] of the original image frame
% 
% use vecs.v for computations like cov or gmm models;
% use vecs in vec2im - to transform results back to image frames

[m,n,N] = size(im_frames);
vecs.v=reshape(im_frames, m*n,N)';
vecs.size=[m n N];