function im_frames=vec2im(vecs)
% function im_frames=vec2im(vecs)
% vecs: structure
% vecs.v: N x (m times n) matrix, where i-th row is a vecorized form of the
%         i-th image
% vecs.size: [m n] of the original image frame
% im_frames: m x n x N array of N image frames, each m x n in size

m=vecs.size(1); n=vecs.size(2); N = vecs.size(3);
im_frames=reshape(vecs.v',m,n,N);