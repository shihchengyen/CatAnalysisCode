load refractory
% break cells into groups for qsub processing
cstart = (@ - 1) * 2 + 1;
cend = @ * 2;
surrogates = generateSurrogates(rf,'cells',cstart:cend,'cstep',100);
save surrogates@ surrogates
exit
