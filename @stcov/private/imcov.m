function [evals, eframes]=imcov(frames)
% function [evals, eframes]=imcov(frames)
% takes a sequence of response-conditioned frames and returns the
% eigenvalues and eigenvectors of the covariance matrix of those frames.
% The covariance is computed on the vectorized frames, where each frames
% is stretched into a column vector. The eigenvectors are transformed back
% to the frame format so they correspond spatially to the frames.
% 
% input: frames, m x n x N samples
% output:
%    evals - list of eigenvalues, sorted, highest first
%    eframes - list of corresponding eigenvectors, in the same format
% as frames: m x n x (mn), ready for frame-wise visualization.
%

vframes = im2vec(frames);

[evecs,evals] = eigs(cov(vframes.v));
[evals,edx] = sort(-diag(evals));
% sort in descending order. V.7 allows for an option, but don't know which
% previous versions support it.
N = min(size(frames,3),length(evals));
% max num of evecs to keep; if nm > N, the remaining ev's will be 0, so no
% need to send them.
evals = -evals(1:N);
vframes.v = evecs(:,edx(1:N))'; 
vframes.size(3)=N;
% replace the vectorized frames with the set of eigenvectors


eframes=vec2im(vframes);
