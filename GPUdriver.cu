#include <stdio.h>
#include <stdlib.h>
#include <vector_types.h>
#include <math.h>
#include <cutil.h>
#include "precision.h"

// includes, kernels
#include <HCq_kernel.cu>

extern "C" {
	
	int gethostname(char *name, size_t len);
	void bzero(void *s, size_t n);
	char *strcpy(char *dest, const char *src);
	
#define MAXDEVICES 4
#define MAXNODES 16
	int cstGPU_init (float* obs_h, int mObs)
	{
		int myproc=0;
		int mydevice=0;
		
		int i, dc;
		cudaError_t cerr ;
		char hostname[64] ;
		struct cudaDeviceProp dp ;
		//  manage devices if multiheaded
		cudaGetDeviceCount( &dc ) ;
		if ( dc > MAXDEVICES ) 
		{ fprintf(stderr, "warning: more than %d devices on node (%d)\n", MAXDEVICES, dc ) ; dc = MAXDEVICES ; }
		fprintf(stderr,"Number of devices on this node: %d\n", dc) ;
		
		// i = *myproc % dc ;
		
		i = mydevice ;
		if ( dc > 0 ) 
		{
			cerr = cudaSetDevice( i );
			if ( cerr ) {
				fprintf(stderr,"    non-zero cerr %d\n",cerr) ;
			}
		}
		gethostname( hostname, 64 ) ;
		fprintf(stderr,"Setting device %02d for task %03d on host %s\n",i,myproc,hostname) ;
		
		cerr = cudaGetDeviceProperties( &dp, i );
		if ( cerr ) {
			fprintf(stderr,"Device %02d: cerr = %d\n",i,cerr) ;
		} else {
			fprintf(stderr,"Device %02d: name %s\n",i,dp.name) ;
			fprintf(stderr,"Device %02d: mem       %d\n",i,(int)dp.totalGlobalMem) ;
			fprintf(stderr,"Device %02d: smem      %d\n",i,(int)dp.sharedMemPerBlock) ;
			fprintf(stderr,"Device %02d: nreg      %d\n",i,dp.regsPerBlock) ;
			fprintf(stderr,"Device %02d: warp      %d\n",i,dp.warpSize) ;
			fprintf(stderr,"Device %02d: pitch     %d\n",i,(int)dp.memPitch) ;
			fprintf(stderr,"Device %02d: maxthrds  %d\n",i,dp.maxThreadsPerBlock) ;
			fprintf(stderr,"Device %02d: maxtdim   %d %d %d\n",i,dp.maxThreadsDim[0]
					,dp.maxThreadsDim[1]
					,dp.maxThreadsDim[2]) ;
			fprintf(stderr,"Device %02d: maxgdim   %d %d %d\n",i,dp.maxGridSize[0]
					,dp.maxGridSize[1]
					,dp.maxGridSize[2]) ;
			fprintf(stderr,"Device %02d: clock     %d\n",i,dp.clockRate) ;
			fprintf(stderr,"Device %02d: talign    %d\n",i,(int)dp.textureAlignment) ;
		}
				
		float BoundaryConditions[9][4] = 
		//	0		1		M-1		M
		{{	-4,		-1,		-1,		-4 },
		{	0,		1,		1,		0 },
		{	2,		-1,		-1,		2 },
		{   -4,     -1,     1,      0 },
		{   -4,     -1,     -1,     2 },
		{   0,      1,      -1,     -4 },
		{   0,      1,      -1,     2 },
		{   2,      -1,     -1,     -4 },
		{   2,      -1,     1,      0 }};

		CUDA_SAFE_CALL(cudaMemcpyToSymbol(BC, BoundaryConditions, sizeof(BoundaryConditions)));

		CUDA_SAFE_CALL(cudaMalloc((void **)&obs_d,mObs*9*sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpy(obs_d,obs_h,mObs*9*sizeof(float),cudaMemcpyHostToDevice)) ;

		CUDA_SAFE_CALL(cudaMalloc((void **)&HCq_d,mObs*sizeof(float)));
				
		return(0) ;
	}
	
	// Preload the coefficients 
	void loadSplineCoeffs_GPU(float* coeffHost, int numCoeffs)
	{
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(coeffDevice,coeffHost,numCoeffs*sizeof(float))) ;
	}
	
	void cstGPU_finalize()
	{
		// clean up memory
		CUDA_SAFE_CALL(cudaFree(obs_d));
		CUDA_SAFE_CALL(cudaFree(HCq_d));
	}
	
	// Try and evaluate a spline on the GPU
	void HCq_GPU(int mObs, float rmax, float rmin, float zmax, float zmin, float* HCq_h, int pState, int zState)
	{
		
		// Constants and boundary conditions
		int R = pState-1;
		float dr = (rmax - rmin) / R;
		float drrecip = 1./dr;
		int Z = zState-1;
		float dz = (zmax - zmin) / Z;
		float dzrecip = 1./dz;
		float onesixth = 1./6.;
		
		/* create and start timer
		unsigned int timer = 0;
		CUT_SAFE_CALL(cutCreateTimer(&timer));
		CUT_SAFE_CALL(cutStartTimer(timer)); */
		
		// setup execution parameters
		int rem = mObs%BLOCKSIZE != 0 ? 1 : 0;
		dim3 threads(BLOCKSIZE);
		dim3 grid(mObs / BLOCKSIZE + rem);
		
		// execute the kernel
		HCq_kernel<<< grid, threads >>>(obs_d, HCq_d, R, Z, rmin, dr, drrecip, zmin, dz, dzrecip, onesixth);
		
		// check if kernel execution generated and error
		CUT_CHECK_ERROR("Kernel execution failed");
		
		// copy result from device to host
		CUDA_SAFE_CALL(cudaMemcpy(HCq_h, HCq_d, mObs*sizeof(float), cudaMemcpyDeviceToHost));
		
		/* stop and destroy timer
		CUT_SAFE_CALL(cutStopTimer(timer));
		printf("Processing time: %f (ms) \n", cutGetTimerValue(timer));
		CUT_SAFE_CALL(cutDeleteTimer(timer)); */
		
	}
	
}
