#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <fftw3.h>

/* Constants */
#define NROWS 380
#define NCOLS 620
#define SIGMA 0.25
#define OUTPUTMEAN 128

struct SPseudoElement {
       long         randvalue;
       int           index;
};

/* Function prototypes */
int sort_SPseudoElement(const void *a, const void *b);

/* Main function */
main(int argc, char **argv)
{
	/* Local variables*/
	FILE *ifp;
	int startframe,endframe,nframes,frame,framesize,nelts,felts,relts,ncols2;
	int i,j,k,index,indim,indim2,n,indexi,indexj;
	float *in,*mag,*phase;
	float sigma2,iwin,jwin,kwin,r,im,m,ph;
	/* float max,min; */
	fftwf_plan p,p2;
	char filename[80];
	unsigned char *image;
	struct SPseudoElement* randelts;

	/* get arguments */
	startframe = atoi(argv[1]);
	endframe = atoi(argv[2]);
	
	/* initialize random seed */
	srandom((unsigned int) time(NULL));
	
	/* compute some variables */
	nframes = endframe - startframe + 1;
	framesize = NROWS * NCOLS;
	nelts = framesize * nframes;
	ncols2 = NCOLS/2 + 1;
	indim = 2 * ncols2;
	indim2 = NROWS * indim;
	sigma2 = pow(SIGMA,2);
	felts = nframes * NROWS * ncols2; /* # of non-redundant fft values */
	relts = felts - 1; /* we skip the DC when scrambling phases */
	
	/* allocate memory */
	in = fftwf_malloc(sizeof(float) * nframes * indim2);
	image = (unsigned char *) malloc(framesize);
	mag = (float *) malloc(sizeof(float) * felts);
	phase = (float *) malloc(sizeof(float) * felts);
	randelts = (struct SPseudoElement *) malloc(sizeof(struct SPseudoElement)*relts);
	
	/* read in images */
	for(frame=startframe,i=0;frame<=endframe;frame++,i++)
	{
		/* calculate i-th index */
		indexi = i * indim2;
		/* calculate i-th component of windowing function */
		iwin = pow(((float)i)/nframes - 0.5,2);
		sprintf(filename,"Cats%05d.pgm",frame);
		printf("Reading %s\n",filename);
		ifp = fopen(filename,"r");
		fscanf(ifp,"%*s %*d %*d %d",&n);
		if(n!=255)
			printf("Warning: Input file format not expected\n");
		n = fread(image,1,framesize,ifp);
		if(n!=framesize)
			printf("Warning: # of elts read not expected\n");
		fclose(ifp);
		for(j=0,index=0;j<NROWS;j++)
		{
			/* calculate j-th index */
			indexj = j * indim;
			/* calculate j-th component of windowing function */
			jwin = pow(((float)j)/NROWS - 0.5,2);
			for(k=0;k<NCOLS;k++,index++)
			{
				/* calculate k-th component of windowing function */
				kwin = pow(((float)k)/NCOLS - 0.5,2);
				/* modulate with windowing function */
				in[indexi+indexj+k] = (float) image[index] * exp(-0.5*(iwin+jwin+kwin)/sigma2);
			}
		}
	}
	
	/* execute fft
	p = fftwf_plan_dft_r2c_3d(int nx, int ny, int nz,
                               double *in, fftwf_complex *out,
                               unsigned flags); */
    printf("Executing fft...");
    p = fftwf_plan_dft_r2c_3d(nframes,NROWS,NCOLS,in,(fftwf_complex *)in,FFTW_ESTIMATE);
	fftwf_execute(p);
	printf("\n");
	
	/* compute magnitude and phase */
	for(i=0;i<felts;i++)
	{
		indexi = i * 2;
		r = in[indexi];
		im = in[indexi+1];
		mag[i] = sqrt(pow(r,2) + pow(im,2));
		phase[i] = atan2(im,r);
	}
                               
	/* initialize structure for random permutation */
	for(i=0;i<relts;i++) 
	{
		randelts[i].randvalue = random();     //get a random value
		randelts[i].index = i+1;
	}
	
	qsort((void *) randelts, relts, sizeof(randelts[0]), sort_SPseudoElement);
		
	/* scramble phases except phase of DC */
	for(i=1;i<felts;i++)
	{
		indexi = i * 2;
		m = mag[i];
		ph = phase[randelts[i-1].index];
		in[indexi] = m*cos(ph);
		in[indexi+1] = m*sin(ph);
	}
	
	/* execute ifft */
    printf("Executing ifft...");
    p2 = fftwf_plan_dft_c2r_3d(nframes,NROWS,NCOLS,(fftwf_complex *)in,in,FFTW_ESTIMATE);
	fftwf_execute(p2);
	printf("\n");

	/* write output images */
	for(frame=startframe,i=0;frame<=endframe;frame++,i++)
	{
		indexi = i * indim2;
		/* for(j=0,index=0,max=0,min=in[0];j<NROWS;j++) */
		for(j=0,index=0;j<NROWS;j++)
		{
			indexj = j * indim;
			for(k=0;k<NCOLS;k++,index++)
			{
				/* divide by nelts to normalize DFT */
				/* if(in[indexi+indexj+k]>max)
					max = in[indexi+indexj+k];
				if(in[indexi+indexj+k]<min)
					min = in[indexi+indexj+k]; */
				image[index] = (unsigned char) (in[indexi+indexj+k]/nelts + OUTPUTMEAN);
			}
		}
		sprintf(filename,"PinkCats%05d.pgm",frame);
		printf("Writing %s\n",filename);
		ifp = fopen(filename,"w");
		fprintf(ifp,"P5\n%d %d 255\n",NCOLS,NROWS);
		n = fwrite(image,1,framesize,ifp);
		fclose(ifp);
		if(n!=framesize)
			printf("Warning: # of elts written not expected\n");
		/* printf("Min: %f Max: %f\n",(float) min/nelts,(float) max/nelts); */
	}

	/* release memory */
	fftwf_destroy_plan(p2);
	fftwf_destroy_plan(p);
	free(randelts);
	free(phase);
	free(mag);
	free(image);
	fftwf_free(in);
}

int sort_SPseudoElement(const void *a, const void *b)
{
    if( (((struct SPseudoElement*) a)->randvalue) > (((struct SPseudoElement*) b)->randvalue) )
        return 1;
    else if( (((struct SPseudoElement*) a)->randvalue) < (((struct SPseudoElement*) b)->randvalue) )
        return -1;
    else
        return 0;
}
