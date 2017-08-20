function newstim = framesexpand(stim2, R,type)
% function newstim = framesexpand(stim2, R, 'type')
% takes an image stack   stim2 and expands it R times by  "type"
% interpolation. R can be non-integer.
% 'type' can be 'nearest' (default), 'bilinear' or 'bicubic'
% interpolation. Each side of the output is R times larger than the
% original.

if nargin<3
	type='nearest';
end
for i1 = size(stim2,3):-1:1
	newstim(:,:,i1) = imresize(stim2(:,:,i1),R,type);
end
