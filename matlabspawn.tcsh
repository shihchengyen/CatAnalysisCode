# script to spawn multiple Matlab processes on bigdog
@ celln = 2
while (${celln} < 41)
	echo "unset DISPLAY; sed "s/@/${celln}/g" event.m | matlab -nojvm > mout${celln} 2>&1" | batch 
	@ celln++
end

# script to spawn multiple Matlab processes on ogier
@ celln = 2
while (${celln} < 41)
	echo "cd matlab/Cat; unset DISPLAY; sed "s/@/${celln}/g" event.m | matlab -nojvm > mout${celln} 2>&1" | qsub -V
	@ celln++
end

# script to spawn multiple Matlab processes on ogier using swarm
@ celln = 2
while (${celln} < 41)
	echo "hostname; unsetenv DISPLAY; sed 's/@/${celln}/g' event.m | matlab -nojvm >&! mout${celln}" > sm${celln}
	swarm -f sm${celln}
	@ celln++
end
