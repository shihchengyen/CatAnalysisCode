#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <malloc.h>
#include <fftw3.h>

/* Constants */
#define NROWS 380
#define NCOLS 620
#define SIGMA 0.25

struct SPseudoElement {
       long         randvalue;
       unsigned int index;
};

/* Function prototypes */
int sort_SPseudoElement(const void *a, const void *b);

/* Main function */
main(int argc, char **argv)
{
	/* Local variables*/
	FILE *ifp;
	unsigned int startframe,endframe,nframes,frame,framesize,nelts,felts,relts,ncols2;
	unsigned int i,j,k,index,indim,indim2,n,indexi,indexj;
	/* long is 64 bits on SGI 64-bit systems so it should be enough to *
	 * store the running sum of all the gray levels in order to compute *
	 * both the input and output mean gray level. */
	unsigned long inmean,outmean;
	double *in,*mag,*phase,*dptr;
	double sigma2,iwin,jwin,kwin,r,im,m,ph;
	double max,min;
	fftw_plan p,p2;
	char filename[80];
	unsigned char *image,ival;
	struct SPseudoElement* randelts;
	double scale;
	
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
	printf("Size of in %ld\n",sizeof(double)*nframes*indim2);
	in = fftw_malloc(sizeof(double) * nframes * indim2);
	if(in==NULL)
	{
		printf("Error allocating memory for in\n");
		return(NULL);
	}
	printf("Size of image %u\n",framesize);
	image = (unsigned char *) malloc(framesize);
	if(image==NULL)
	{
		printf("Error allocating memory for image\n");
		fftw_free(in);
		return(NULL);
	}		
	printf("Size of mag %ld\n",sizeof(double)*felts);
	mag = (double *) malloc(sizeof(double) * felts);
	if(mag==NULL)
	{
		printf("Error allocating memory for mag\n");
		free(image);
		fftw_free(in);
		return(NULL);
	}
	printf("Size of phase %ld\n",sizeof(double)*felts);
	phase = (double *) malloc(sizeof(double) * felts);
	if(phase==NULL)
	{
		printf("Error allocating memory for phase\n");
		free(mag);
		free(image);
		fftw_free(in);
		return(NULL);
	}
	printf("Size of randelts %ld\n",sizeof(struct SPseudoElement)*relts);
	randelts = (struct SPseudoElement *) malloc(sizeof(struct SPseudoElement)*relts);
	if(randelts==NULL)
	{
		printf("Error allocating memory for randelts\n");
		free(phase);
		free(mag);
		free(image);
		fftw_free(in);
		return(NULL);
	}
	
	/* loop to catch instances when min is some large negative number */
	/* initialize scale to 0 so we will enter the loop at least once */
	/* scale gets initialized to 1 later when we actually look for the */
	/* proper scale */
	for(scale=0;scale<0.1;)
	{
		
		/* read in images */
		for(frame=startframe,i=0,inmean=0;frame<=endframe;frame++,i++)
		{
			/* calculate i-th index */
			indexi = i * indim2;
			/* calculate i-th component of windowing function */
			iwin = pow(((double)i)/nframes - 0.5,2);
			sprintf(filename,"Cats%05d.pgm",frame);
			printf("Reading %s\n",filename);
			ifp = fopen(filename,"r");
			fscanf(ifp,"%*s %*d %*d %d%*c",&n);
			if(n!=255)
				printf("Warning: Input file format not expected\n");
			n = fread(image,1,framesize,ifp);
			if(n!=framesize)
				printf("Warning: # of elts read not expected\n");
			/* try reading one more byte to make sure we are at the end of the file */
			if(fread(&ival,1,1,ifp)!=0)
				printf("Extraneous bytes in input file\n");
			fclose(ifp);
			for(j=0,index=0;j<NROWS;j++)
			{
				/* calculate j-th index */
				indexj = indexi + (j * indim);
				/* calculate j-th component of windowing function */
				jwin = iwin + pow(((double)j)/NROWS - 0.5,2);
				for(k=0;k<NCOLS;k++,index++)
				{
					/* calculate k-th component of windowing function */
					kwin = pow(((double)k)/NCOLS - 0.5,2);
					/* modulate with windowing function */
					ival = image[index];
					inmean += ival;
					in[indexj+k] = (double) ival * exp(-0.5*(jwin+kwin)/sigma2);
				}
			}
		}
		
		/* compute input mean */
		inmean = rint((double) inmean/nelts);
		
		/* initialize threads */
		if(fftw_init_threads()==0)
		{
			printf("Error initializing threads\n");
			free(randelts);
			free(phase);
			free(mag);
			free(image);
			fftw_free(in);
			return(NULL);
		}
				
		/* set number of threads to use */
		printf("Setting up threads\n");
		fftw_plan_with_nthreads(8);
		
		/* execute fft
		p = fftw_plan_dft_r2c_3d(int nx, int ny, int nz,
								   double *in, fftw_complex *out,
								   unsigned flags); */
		printf("Planning fft...\n");
		p = fftw_plan_dft_r2c_3d(nframes,NROWS,NCOLS,in,(fftw_complex *)in,FFTW_ESTIMATE);
		if(p==NULL)
		{
			printf("Error creating FFT plan\n");
			/* release memory */
			fftw_destroy_plan(p);
			fftw_cleanup_threads();
			free(randelts);
			free(phase);
			free(mag);
			free(image);
			fftw_free(in);
			return(NULL);
		}
		else
			printf("Plan created, executing fft...");
		fftw_execute(p);
		printf("Done!\n");
		
		/* compute magnitude and phase */
		for(i=0;i<felts;i++)
		{
			indexi = i * 2;
			r = in[indexi];
			im = in[indexi+1];
			mag[i] = sqrt(pow(r,2) + pow(im,2));
			phase[i] = atan2(im,r);
			if(i<relts)
			{
				/* initialize structure for random permutation */
				randelts[i].randvalue = random();     //get a random value
				randelts[i].index = i+1;
			}
		}
								   
		qsort((void *) randelts, relts, sizeof(randelts[0]), sort_SPseudoElement);
		
		/* scramble phases except phase of DC */
		for(i=1;i<felts;i++)
		{
			indexi = i * 2;
			m = mag[i];
			ph = phase[randelts[i-1].index];
			/* ph = phase[i]; */
			in[indexi] = m*cos(ph);
			in[indexi+1] = m*sin(ph);
		}
		
		/* execute ifft */
		printf("Planning ifft...\n");
		p2 = fftw_plan_dft_c2r_3d(nframes,NROWS,NCOLS,(fftw_complex *)in,in,FFTW_ESTIMATE);
		if(p2==NULL)
		{
			printf("Error creating FFT plan\n");
			/* release memory */
			fftw_destroy_plan(p2);
			fftw_destroy_plan(p);
			fftw_cleanup_threads();
			free(randelts);
			free(phase);
			free(mag);
			free(image);
			fftw_free(in);
			return(NULL);
		}
		else
			printf("Plan created, executing ifft...");
		fftw_execute(p2);
		printf("Done!\n");
		
		/* compute output mean */
		for(frame=startframe,i=0,outmean=0,max=0,min=in[0],index=0;frame<=endframe;frame++,i++)
		/* for(frame=startframe,i=0,outmean=0;frame<=endframe;frame++,i++) */
		{
			indexi = i * indim2;
			for(j=0;j<NROWS;j++)
			{
				indexj = indexi + (j * indim);
				for(k=0;k<NCOLS;k++,index++)
				{
					/* get pointer so we don't have to keep doing pointer arithmetric */
					dptr = &in[indexj+k];
					/* divide by nelts to normalize DFT */
					*dptr /= nelts;
					outmean += *dptr;
					if(*dptr>max)
						max = *dptr;
					if(*dptr<min)
						min = *dptr;
				}
			}
		}
		printf("Computed output mean\n");
		
		/* compute output mean */
		outmean = rint((double) outmean/nelts);
	
		/* compute how much to change output values */
		ival = inmean - outmean;
		
		/* initialize scale values to 1 since if they get computed they will be *
		 * smaller than 1 */
		scale = 1; 
		 
		/* if after shifting the mean, min is still negative, we have to *
		 * scale the values to make sure it is not smaller than 0 */
		if(min+ival<0)
		{
			/* min should be negative so ival will be positive */
			ival = -min;
			/* take absolute value to make sure scale does not come out negative */
			scale = fabs(inmean/(outmean-min));
		}
		/* if after shifting the mean, max is larger than 255, we have to *
		 * scale the values to make sure it is not larger than 255
		if(max+ival>255)
		{
			scale2 = (255-ival)/max;
		} */
		/* take the smaller of the two scales */
		/* if(scale1<scale2)
			scale = scale1;
		else
			scale = scale2; */

		printf("Mean (In: %lu, Out: %lu), Max: %lf, Min: %lf, Scale: %lf\n",\
			inmean,outmean,max,min,scale);
	}

	/* write output images */
	for(frame=startframe,i=0;frame<=endframe;frame++,i++)
	{
		indexi = i * indim2;
		for(j=0,index=0;j<NROWS;j++)
		{
			indexj = indexi + (j * indim);
			for(k=0;k<NCOLS;k++,index++)
			{
				image[index] = (unsigned char) ((in[indexj+k] + ival) * scale);
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
	}

	/* release memory */
	fftw_destroy_plan(p2);
	fftw_destroy_plan(p);
	fftw_cleanup_threads();
	free(randelts);
	free(phase);
	free(mag);
	free(image);
	fftw_free(in);
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
