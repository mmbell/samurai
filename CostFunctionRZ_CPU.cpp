/*
 *  CostFunctionRZ_CPU.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionRZ_CPU.h"
#include <cmath>

/* First, compute HCq and Hxb via recursive filter on q, followed by loop over M observations
 yielding 2 y-estimate M x 1 vectors. Solve for innovation vector d by subtracting from y.
 For function, subtract from HCq and take inner product (with R-1), add to i.p. of q vector.
 For gradient, loop over M again, this time summing product of operator*basis*oberror-1*y-hat
 and operator*basis*oberror1*d for each coefficient N, yielding an N x 1 vector. Filter this
 vector and the innovation vector. Add 1+vector1+vector2 to get the gradient. Treat the penalty
 constraints as observations. */

CostFunctionRZ_CPU::CostFunctionRZ_CPU(const int& numObs, const int& stateSize)
	: CostFunction(numObs, stateSize)
{
}

CostFunctionRZ_CPU::~CostFunctionRZ_CPU()
{
}

void CostFunctionRZ_CPU::finalize()
{	

	delete iFilter;
	delete jFilter;
	delete[] currState;
	delete[] currGradient;
	delete[] innovation;
	delete[] tempState;
	delete[] tempGradient;
	delete[] HCq;
	delete[] xt;
	delete[] df;
	delete[] stateA;
	delete[] stateB;
	delete[] CTHTd;
	delete[] iTemp;
	delete[] jTemp;
	
}

void CostFunctionRZ_CPU::initialize(const real& imin, const real& imax, const int& idim,
									const real& jmin, const real& jmax, const int& jdim,
									const real* ia, const real* ja, real* bgB, real* obs, SplineD* zs,
									const unsigned int* IXdim, const vector<real>* IX,
									const unsigned int* JXdim, const vector<real>* JX)
{

	// Initialize number of control variables -- one less than physical variables
	varDim = 5;
	
	// Initialize background errors and filter scales
	bgError[0] = 10.;
	bgError[1] = 500.;
	bgError[2] = 300.;
	bgError[3] = 20.;
	bgError[4] = 0.01;	
	int iFilterScale = 3;
	int jFilterScale = 2;

	// Assign local object pointers
	bgState = bgB;
	obsVector = obs;
	iMin = imin;
	iMax = imax;
	iDim = idim;
	jMin = jmin;
	jMax = jmax;
	jDim = jdim;
	iA = ia;
	jA = ja;
	IXDim = IXdim;
	IXform = IX;
	if (JXdim) JXDim = JXdim;
	if (JX) JXform = JX;
	DI = (iMax - iMin) / (iDim - 1);
	DIrecip = 1./DI;
	DJ = (jMax - jMin) / (jDim - 1);
    DJrecip = 1./DJ;
	
	// Set up the recursive filter
	iFilter = new RecursiveFilter(4,iFilterScale);
	jFilter = new RecursiveFilter(4,jFilterScale);
	//iFilter = new RecursiveFilter(4,4);
	//jFilter = new RecursiveFilter(4,2);

	
	/*iFilterCoeff = iFilter->getFilterCoefficients();
	iFilterBC = iFilter->getSn();
	jFilterCoeff = jFilter->getFilterCoefficients();
	jFilterBC = jFilter->getSn(); */
	
	// Allocate memory for the needed arrays
	// These are common to all CostFunctions
	currState = new real[nState];
	currGradient = new real[nState];
	tempState = new real[nState];
	tempGradient = new real[nState];
	xt = new real[nState];
	df = new real[nState];
	innovation = new real[mObs];	

	// These are local to this one
	HCq = new real[mObs];
	stateB = new real[nState];
	stateA = new real[nState];
	CTHTd = new real[nState];
	iTemp = new real[iDim];
	jTemp = new real[jDim];
	zSpline = zs;
	initState();
}	

void CostFunctionRZ_CPU::initState()
{
	
	// Compute and display the variable RMS and BG errors for reference
	for (int var = 0; var < varDim; var++) {
		double varScale = 0;
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				varScale += bgState[varDim*iDim*jIndex +varDim*iIndex + var] * bgState[varDim*iDim*jIndex +varDim*iIndex	+ var];
			}
		}
		varScale = sqrt(varScale/(iDim*jDim));
		if (varScale) {
			double errPct = 100*bgError[var]/varScale;
			cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var] 
			<< " ( " << errPct << " %)" << endl;
		} else {
			cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var] 
			<< " ( Infinite! %)" << endl;
		}
	}
	
	// Clear the state vector
	cout << "Initializing state vector..." << endl;
	for (int n = 0; n < nState; n++) {
		currState[n] = 0.0;
		currGradient[n] = 0.0;
		tempState[n] = 0.0;
		tempGradient[n] = 0.0;
		xt[n] = 0.0;
		df[n] = 0.0;
		stateA[n] = 0.0;
		stateB[n] = 0.0;
	}
		
	// SA transform = bg B's -> bg A's
	SAtransform(bgState, stateA);
	
	// d = y - HXb
	calcInnovation();
	
	// HTd
	calcHTranspose(innovation, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtransform(stateA, stateB);

	// D^T
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] *= bgError[var];
			}
		}
	}
	
	for (int var = 0; var < varDim; var++) {
		//FI
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iTemp[iIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			iFilter->filterArray(iTemp, iDim);
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex]; 
			}
		}
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				CTHTd[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
	}
	
}	

double CostFunctionRZ_CPU::funcValue(double* state)
{
	// Update the Y hat vector
	updateHCq(state);

	double qIP, obIP;
	qIP = 0;
	obIP = 0;
	// Compute inner product of state vector
	for (int n = 0; n < nState; n++) {
		qIP += state[n]*state[n];
	}
		
	// Subtract d from HCq to yield mObs length vector and compute inner product
	#pragma omp parallel for reduction(+:obIP)
	for (int m = 0; m < mObs; m++) {
		obIP += (HCq[m]-innovation[m])*(obsVector[m*10+1])*(HCq[m]-innovation[m]);
	}
	double J = 0.5*qIP + 0.5*obIP;
	return J;
	
}

void CostFunctionRZ_CPU::funcGradient(double* state, double* gradient)
{
	
	updateHCq(state);
	
	// HTHCq
	calcHTranspose(HCq, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtransform(stateA, stateB);
	
	// D^T
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] *= bgError[var];
			}
		}
	}
	
	for (int var = 0; var < varDim; var++) {
		//FI
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iTemp[iIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			iFilter->filterArray(iTemp, iDim);
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex]; 
			}
		}
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
	}
	
	for (int n = 0; n < nState; n++) {
		gradient[n] = state[n] + stateB[n] - CTHTd[n];
	}
	
	
}

void CostFunctionRZ_CPU::updateHCq(double* state)
{
		
	for (int var = 0; var < varDim; var++) {
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = state[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
		//FI
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iTemp[iIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			iFilter->filterArray(iTemp, iDim);
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex]; 
			}
		}
	}
	
	// D
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] *= bgError[var];
			}
		}
	}

	// S (SA transform) yield A's, put it in the tempState
	SAtransform(stateB, stateA);

	// H
	#pragma omp parallel for
	for (int m = 0; m < mObs; m++) {
		int mi = m*10;
		real w1 = obsVector[mi+2];
		real w2 = obsVector[mi+3];
		real w3 = obsVector[mi+4];
		real w4 = obsVector[mi+5];
		real w5 = obsVector[mi+6];
		real w6 = obsVector[mi+7];		
		real i = obsVector[mi+8];
		real j = obsVector[mi+9];
		real invI = 1./i;
		real tempsum = 0;
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		real ibasis = 0;
		real jbasis = 0;
		real idbasis = 0;
		real jdbasis = 0;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {				
				if ((iNode < 0) or (iNode >= iDim) or (jNode < 0) or (jNode >= jDim)) continue;
				// If conditions may cause parallel warps to diverge -- maybe actually faster to just calculate them all but need to test
				if(w1) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * w1;
				}
				if(w2 or w3) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					idbasis = DBasis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jdbasis = DBasis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					float coeff = stateA[varDim*iDim*jNode +varDim*iNode + 1];
					tempsum += coeff * ibasis * (-jdbasis) * w2 * 1e3 * invI;
					tempsum += coeff * idbasis * jbasis * w3 * invI;
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis * w4;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * w5;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * w6;
				}
			}
		}
		HCq[m] = tempsum;
	}	

	
}

// Basis Functions
float CostFunctionRZ_CPU::Basis(int m, float x, int M, float xmin, 
					   float DX, float DXrecip, int C)
{

	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
	float ONESIXTH = 0.16666666666666666666666666667;
	float BC[9][4] = 
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
	
	if (z < 2.0)
	{
		z = 2 - z;
		b = (z*z*z) * ONESIXTH;
		z -= 1.0;
		if (z > 0)
			b -= (z*z*z) * 4 * ONESIXTH;
	}
	
	// Boundary conditions, if any, are an additional addend.
	if (C >= 0) {
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
	}
	return b;
}

float CostFunctionRZ_CPU::DBasis(int m, float x, int M, float xmin, 
						float DX, float DXrecip, int C)
{
	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
	float ONESIXTH = 0.16666666666666666666666666667;
	float BC[9][4] = 
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
	if (C >= 0) {
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
	}
	return b;
}


void CostFunctionRZ_CPU::updateBG()
{

	for (int var = 0; var < varDim; var++) {
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = currState[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
		//FI
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iTemp[iIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			iFilter->filterArray(iTemp, iDim);
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex]; 
			}
		}
	}
	
	// D
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] *= bgError[var];
			}
		}
	}
	
	// No SA transform on BG update since we are directly summing B's
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				bgState[varDim*iDim*jIndex +varDim*iIndex + var] += stateB[varDim*iDim*jIndex +varDim*iIndex + var];	
			}
		}
	}

}

void CostFunctionRZ_CPU::calcInnovation()
{
	// Initialize and fill the innovation vector
	cout << "Initializing innovation vector..." << endl;
	for (int m = 0; m < mObs; m++) {
		HCq[m] = 0.0;
		innovation[m] = obsVector[m*10];
	}
	
	#pragma omp parallel for
	for (int m = 0; m < mObs; m++) {
		int mi = m*10;
		real w1 = obsVector[mi+2];
		real w2 = obsVector[mi+3];
		real w3 = obsVector[mi+4];
		real w4 = obsVector[mi+5];
		real w5 = obsVector[mi+6];
		real w6 = obsVector[mi+7];		
		real i = obsVector[mi+8];
		real j = obsVector[mi+9];
		real invI = 1./i;
		real tempsum = 0;
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		real ibasis = 0;
		real jbasis = 0;
		real idbasis = 0;
		real jdbasis = 0;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {				
				if ((iNode < 0) or (iNode >= iDim) or (jNode < 0) or (jNode >= jDim)) continue;
				// If conditions may cause parallel warps to diverge -- maybe actually faster to just calculate them all but need to test
				if(w1) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * w1;
				}
				if(w2 or w3) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					idbasis = DBasis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jdbasis = DBasis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					float coeff = stateA[varDim*iDim*jNode +varDim*iNode + 1];
					tempsum += coeff * ibasis * (-jdbasis) * w2 * 1e3 * invI;
					tempsum += coeff * idbasis * jbasis * w3 * invI;
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis * w4;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * w5;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * w6;
				}
			}
		}
		innovation[m] -= tempsum;
	}
}	

void CostFunctionRZ_CPU::calcHTranspose(const real* yhat, real* Astate) 
{
	
	// Clear the Astate
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = 0.;
			}
		}
	}
	
	// Calculate H Transpose	
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			#pragma omp parallel for
			for (int m = 0; m < mObs; m++) {
				// Sum over obs this time
				// Multiply state by H weights
				int mi = m*10;
				real invError = obsVector[mi+1];
				real w1 = obsVector[mi+2];
				real w2 = obsVector[mi+3];
				real w3 = obsVector[mi+4];
				real w4 = obsVector[mi+5];
				real w5 = obsVector[mi+6];
				real w6 = obsVector[mi+7];		
				real i = obsVector[mi+8];
				real j = obsVector[mi+9];
				real invI = 1./i;
				int ii = (int)((i - iMin)*DIrecip);
				if ((iIndex < ii-1) or (iIndex > ii+2)) continue;
				int jj = (int)((j - jMin)*DJrecip);
				if ((jIndex < jj-1) or (jIndex > jj+2)) continue;
				real ibasis = 0;
				real jbasis = 0;
				real idbasis = 0;
				real jdbasis = 0;
				
				if(w1) {
					ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 2);
					Astate[varDim*iDim*jIndex +varDim*iIndex] 
						+= yhat[m] * ibasis * jbasis * w1 * invError;
				}
				if(w2 or w3) {
					ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 4);
					idbasis = DBasis(iIndex, i, iDim-1, iMin, DI, DIrecip, 4);
					jdbasis = DBasis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 4);
					Astate[varDim*iDim*jIndex +varDim*iIndex + 1] 
						+= yhat[m] * ibasis * (-jdbasis) * w2 * 1e3 * invI * invError;
					Astate[varDim*iDim*jIndex +varDim*iIndex + 1]
						+= yhat[m] * idbasis * jbasis * w3 * invI * invError;
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 2);
					jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 2);
					Astate[varDim*iDim*jIndex +varDim*iIndex + 2] 
						+= yhat[m] * ibasis * jbasis * w4 * invError;
					Astate[varDim*iDim*jIndex +varDim*iIndex + 3] 
						+= yhat[m] * ibasis * jbasis * w5 * invError;
					Astate[varDim*iDim*jIndex +varDim*iIndex + 4]
						+= yhat[m] * ibasis * jbasis * w6 * invError;
				}
				
			}
		}
	}
	
}

bool CostFunctionRZ_CPU::SAtransform(real* Bstate, real* Astate)
{
	
	unsigned int bands = 3;
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			real* jB = new real[jDim];
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jB[jIndex] = Bstate[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			
			// Forward substitution to find y.  The diagonals of the lower
			// triangular matrix are taken to be 1.
			unsigned int i, j;
			unsigned int M = jDim;
			real sum;
			for (i = 2; i <= M; ++i)
			{
				sum = jB[i-1];
				for (j = (i > bands) ? i-bands : 1; j < i; ++j)
				{
					sum -= jA[(i-1)*jDim+j-1]*jB[j-1];
				}
				jB[i-1] = sum;
			}
			
			// Now for the backward substitution
			jB[M-1] /= jA[(M-1)*jDim+M-1];
			for (i = M-1; i >= 1; --i)
			{
				if (jA[(i-1)*jDim+i-1] == 0)	// oops!
					return false;
				sum = jB[i-1];
				for (j = i+1; (j <= M) && (j <= i+bands); ++j)
				{
					sum -= jA[(i-1)*jDim+j-1]*jB[j-1];
				}
				jB[i-1] = sum / jA[(i-1)*jDim+i-1];
			}
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = jB[jIndex]; 
			}
			delete[] jB;
		}
		
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			real* iB = new real[iDim];
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iB[iIndex] = Astate[varDim*iDim*jIndex +varDim*iIndex + var];
			}

			// Forward substitution to find y.  The diagonals of the lower
			// triangular matrix are taken to be 1.
			unsigned int i, j;
			unsigned int M = iDim;
			real sum;
			for (i = 2; i <= M; ++i)
			{
				sum = iB[i-1];
				for (j = (i > bands) ? i-bands : 1; j < i; ++j)
				{
					sum -= iA[(i-1)*iDim+j-1]*iB[j-1];
				}
				iB[i-1] = sum;
			}
			
			// Now for the backward substitution
			iB[M-1] /= iA[(M-1)*iDim+M-1];
			for (i = M-1; i >= 1; --i)
			{
				if (iA[(i-1)*iDim+i-1] == 0)	// oops!
					return false;
				sum = iB[i-1];
				for (j = i+1; (j <= M) && (j <= i+bands); ++j)
				{
					sum -= iA[(i-1)*iDim+j-1]*iB[j-1];
				}
				iB[i-1] = sum / iA[(i-1)*iDim+i-1];
			}
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = iB[iIndex]; 
			}
			delete[] iB;
		}		
	}

	return true;
}		



