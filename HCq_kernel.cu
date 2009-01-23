#include <stdio.h>
#ifndef _CSTGPU_KERNEL_H
#define _CSTGPU_KERNEL_H
#define BLOCKSIZE 6

// Declare some variables
__constant__ float BC[9][4];
__constant__ float coeffDevice[8500];
float* obs_d;
float* HCq_d;

// Basis Functions
__device__ float Basis(int m, float x, int M, float xmin, 
						float DX, float DXrecip, float ONESIXTH, int C)
{
	
	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
	
	if (z < 2.0)
	{
		z = 2 - z;
		b = (z*z*z) * ONESIXTH;
		z -= 1.0;
		if (z > 0)
			b -= (z*z*z) * 4 * ONESIXTH;
	}
	
	// Boundary conditions, if any, are an additional addend.
	if (m == 0 || m == 1) {
		float l = 0;
		xm = xmin + (-1 * DX);
		delta = (x - xm) * DXrecip;
		z = fabsf(delta);
		
		if (z < 2.0)
		{
			z = 2 - z;
			l = (z*z*z) * ONESIXTH;
			z -= 1.0;
			if (z > 0)
				l -= (z*z*z) * 4 * ONESIXTH;
		}
		b += BC[C][m] * l;
	} else if (m == M-1 || m == M) {
		float r = 0;
		xm = xmin + ((M+1) * DX);
		delta = (x - xm) * DXrecip;
		z = fabsf(delta);
		
		if (z < 2.0)
		{
			z = 2 - z;
			r = (z*z*z) * ONESIXTH;
			z -= 1.0;
			if (z > 0)
				r -= (z*z*z) * 4 * ONESIXTH;
		}
		b += BC[C][m+3-M] * r;
	}
	return b;
}

__device__ float DBasis(int m, float x, int M, float xmin, 
						float DX, float DXrecip, float ONESIXTH, int C)
{
	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
	
	if (z < 2.0)
	{
		z = 2.0 - z;
		b = (z*z) * ONESIXTH;
		z -= 1.0;
		if (z > 0)
			b -= (z*z) * 4 * ONESIXTH;
		b *= ((delta > 0) ? -1.0 : 1.0) * 3.0 / DX;
	}
	
	// Boundary conditions, if any, are an additional addend.
	if (m == 0 || m == 1) {
		float l = 0;
		xm = xmin + (-1 * DX);
		delta = (x - xm) * DXrecip;
		z = fabsf(delta);
		
		if (z < 2.0)
		{
			z = 2 - z;
			l = (z*z) * ONESIXTH;
			z -= 1.0;
			if (z > 0)
				l -= (z*z) * 4 * ONESIXTH;
			l *= ((delta > 0) ? -1.0 : 1.0) * 3.0 / DX;
		}
		
		b += BC[C][m] * l;
	} else if (m == M-1 || m == M) {
		float r = 0;
		xm = xmin + ((M+1) * DX);
		delta = (x - xm) * DXrecip;
		z = fabsf(delta);
		
		if (z < 2.0)
		{
			z = 2 - z;
			r = (z*z) * ONESIXTH;
			z -= 1.0;
			if (z > 0)
				r -= (z*z) * 4 * ONESIXTH;
			r *= ((delta > 0) ? -1.0 : 1.0) * 3.0 / DX;	
		}
		b += BC[C][m+3-M] * r;
	}
	return b;
}

__global__ void HCq_kernel(float* obs_d,float* HCq_d, int R, int Z, float rmin, float DR, float DRrecip, 
						   float zmin, float DZ, float DZrecip, float ONESIXTH)
{

	// Block and thread indices
	//int bx = blockIdx.x;
	//int tx = threadIdx.x;
	int xi = blockIdx.x*BLOCKSIZE + threadIdx.x;
	int mi = xi*9;
	int R1 = R+1;
	float w1 = obs_d[mi];
	float w2 = obs_d[mi+1];
	float w3 = obs_d[mi+2];
	float w4 = obs_d[mi+3];
	float w5 = obs_d[mi+4];
	float w6 = obs_d[mi+5];
	float radius = obs_d[mi+6];
	float height = obs_d[mi+7];
	float invRadius = 1./radius;
	float HCq = 0;
	int m = (int)((radius - rmin)*DRrecip);
	int n = (int)((height - zmin)*DZrecip);
	float bz = 0;
	float br = 0;
	float bzp = 0;
	float brp = 0;
	int bc = 1;
	// rhoV = BC_LZERO_RSECOND, r & BC_ZERO_SECOND, z
	for (int r = m-1; r <= m+2; ++r) {
		for (int z = n-1; z <= n+2; ++z) {				
			if ((r < 0) or (r > R) or (z < 0) or (z > Z)) continue;
			if ((r > 1) and (r < R-1) and (z > 1) and (z < Z-1)) {
				// No BCs to worry about, calculate the basis once
				br = Basis(r, radius, R, rmin, DR, DRrecip, ONESIXTH, 2);
				bz = Basis(z, height, Z, zmin, DZ, DZrecip, ONESIXTH, 2);
				brp = DBasis(r, radius, R, rmin, DR, DRrecip, ONESIXTH, 4);
				bzp = DBasis(z, height, Z, zmin, DZ, DZrecip, ONESIXTH, 4);
				//printf("%d,  %d, %f, %f\n", r, z, br, bz);
				bc = 0;
				HCq += coeffDevice[z*5*R1 + r*5] * br * bz * w1;
				float coeff = coeffDevice[z*5*R1 + r*5 +1];
				HCq += coeff * br * (-bzp) * w2 * 1e3 * invRadius;
				HCq += coeff * brp * bz * w3 * invRadius;
				HCq += coeffDevice[z*5*R1 + r*5 +2] * br * bz * w4;
				HCq += coeffDevice[z*5*R1 + r*5 +3] * br * bz * w5;
				HCq += coeffDevice[z*5*R1 + r*5 +4] * br * bz * w6;
			} else {
				if (w1) {
					if (bc) { 
						br = Basis(r, radius, R, rmin, DR, DRrecip, ONESIXTH, 4);
						bz = Basis(z, height, Z, zmin, DZ, DZrecip, ONESIXTH, 2);
					}
					HCq += coeffDevice[z*5*R1 + r*5] * br * bz * w1;
				}
				float coeff;
				if (w2 or w3) coeff = coeffDevice[z*5*R1 + r*5 +1];
				if (w2) {
					if (bc) {
						br = Basis(r, radius, R, rmin, DR, DRrecip, ONESIXTH, 4);
						bzp = DBasis(z, height, Z, zmin, DZ, DZrecip, ONESIXTH, 4);
					}
					HCq += coeff * br * (-bzp) * w2 * 1e3 * invRadius;
				}
				if (w3) {
					if (bc) {
						brp = DBasis(r, radius, R, rmin, DR, DRrecip, ONESIXTH, 4);
						bz = Basis(z, height, Z, zmin, DZ, DZrecip, ONESIXTH, 4);
					}
					HCq += coeff * brp * bz * w3 * invRadius;
				}
				if (w4 or w5 or w6) {
					if (bc) {
						br = Basis(r, radius, R, rmin, DR, DRrecip, ONESIXTH, 2);
						bz = Basis(z, height, Z, zmin, DZ, DZrecip, ONESIXTH, 2);
					}
				}
				if (w4)
					HCq += coeffDevice[z*5*R1 + r*5 +2] * br * bz * w4;
				if (w5)
					HCq += coeffDevice[z*5*R1 + r*5 +3] * br * bz * w5;
				if (w6)
					HCq += coeffDevice[z*5*R1 + r*5 +4] * br * bz * w6;
			}
		}
	}
	HCq_d[xi] = HCq;
	//printf("%d : %f\n", xi,  HCq);
}

#endif
