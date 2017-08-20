phaserand: phaserand.c
	cc -64 -O2 phaserand.c -o phaserand -I/ccb/home/syen/lib/fftw-3.0.1/api \
	-L/ccb/home/syen/lib/fftw-3.0.1/.libs -L/ccb/home/syen/lib/fftw-3.0.1/threads/.libs \
	-lfftw3_threads -lfftw3 -lpthread -lmalloc -lm

phaserandf: phaserandf.c
	cc -O phaserandf.c -o phaserandf -I/Users/syen/Documents/ShihCheng/lib/fftw-3.0.1-fma/api \
	-L/Users/syen/Documents/ShihCheng/lib/fftw-3.0.1-fma/.libs -lfftw3f -lm -lc
