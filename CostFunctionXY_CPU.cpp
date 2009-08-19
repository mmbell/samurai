/*
 *  CostFunctionXY_CPU.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionXY_CPU.h"
#include "MetObs.h"
#include <cmath>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QFile>
#include <QDir>

/* First, compute HCq and Hxb via recursive filter on q, followed by loop over M observations
 yielding 2 y-estimate M x 1 vectors. Solve for innovation vector d by subtracting from y.
 For function, subtract from HCq and take inner product (with R-1), add to i.p. of q vector.
 For gradient, loop over M again, this time summing product of operator*basis*oberror-1*y-hat
 and operator*basis*oberror1*d for each coefficient N, yielding an N x 1 vector. Filter this
 vector and the innovation vector. Add 1+vector1+vector2 to get the gradient. Treat the penalty
 constraints as observations. */

CostFunctionXY_CPU::CostFunctionXY_CPU(const int& numObs, const int& stateSize)
	: CostFunction(numObs, stateSize)
{
}

CostFunctionXY_CPU::~CostFunctionXY_CPU()
{
}

void CostFunctionXY_CPU::finalize()
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
	delete[] stateC;
	delete[] CTHTd;
	delete[] iTemp;
	delete[] jTemp;
	delete[] iL;
	delete[] jL;
	delete[] stateU;
	delete[] Uprime;
}

void CostFunctionXY_CPU::initialize(const real& imin, const real& imax, const int& idim,
									const real& jmin, const real& jmax, const int& jdim,
									const real* ia, const real* ja, real* bgU, real* obs,
									const unsigned int* IXdim, const vector<real>* IX,
									const unsigned int* JXdim, const vector<real>* JX)
{

	// Initialize number of control variables -- one less than physical variables
	varDim = 5;
	
	// Initialize background errors and filter scales
	bgError[0] = 0.02;
	bgError[1] = 0.1;
	bgError[2] = 5.0;
	bgError[3] = 3.0;
	bgError[4] = 3.0;	

	LI = 2.5;
	LJ = 250;

	LF[0] = 0.1;
	LF[1] = 0.5;
	LF[2] = 0.5;
	LF[3] = 0.5;
	LF[4] = 0.5;
	
	/* bcLeft[0] = R1T0; bcRight[0] = R1T2; bcTop[0] = R1T2; bcBottom[0] = R1T2;
	bcLeft[1] = R1T0; bcRight[1] = R1T2; bcTop[1] = R1T0; bcBottom[1] = R1T2;
	bcLeft[2] = R1T2; bcRight[2] = R1T2; bcTop[2] = R1T2; bcBottom[2] = R1T2;
	bcLeft[3] = R1T2; bcRight[3] = R1T2; bcTop[3] = R1T2; bcBottom[3] = R1T2;
	bcLeft[4] = R1T2; bcRight[4] = R1T2; bcTop[4] = R1T2; bcBottom[4] = R1T2; */
		
	rhoBase = 1.156;
	rhoInvScaleHeight = 9.9504e-5;

	// Assign local object pointers
	bgFields = bgU;
	rawObs = obs;
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
	real iFilterScale = 4.;
	real jFilterScale = 4.;
	//bgErrorScale = 2 * 3.141592653589793 *iFilterScale*jFilterScale;
	bgErrorScale = 1.;
	iFilter = new RecursiveFilter(4,iFilterScale);
	jFilter = new RecursiveFilter(4,jFilterScale);

	/* Test the filter
	real* p = new real[iDim];
	for (int i = 0; i < iDim; i++) {
		p[i] = 0;
	}
	
	p[iDim/2] = 1.;
	real Pi = 3.141592653589793;
	iFilter->filterArray(p, iDim);
	for (int i = 0; i < iDim; i++) {
		real ibasis =  Basis(i, iDim/2, iDim-1, iMin, DI, DIrecip, 2);
		cout << i << "\t" << p[i] << "\t" << ibasis << endl;
	} */
	
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
	bState = iDim*jDim*varDim;
	CTHTd = new real[nState];
	stateU = new real[nState];
	Uprime = new real[nState];
	HCq = new real[mObs];
	obsVector = new real[mObs*11];
	
	bgState = new real[bState];
	stateB = new real[bState];
	stateA = new real[bState];
	stateC = new real[bState];
	iTemp = new real[iDim];
	jTemp = new real[jDim];
	
	// Set up the spline matrices
	setupSplines();

}	

void CostFunctionXY_CPU::initState()
{

	// Clear the state vector
	cout << "Initializing state vector..." << endl;
	for (int n = 0; n < nState; n++) {
		currState[n] = 0.0;
		currGradient[n] = 0.0;
		tempState[n] = 0.0;
		tempGradient[n] = 0.0;
		xt[n] = 0.0;
		df[n] = 0.0;
	}
	for (int b = 0; b < bState; b++) {
		bgState[b] = 0.0;
		stateA[b] = 0.0;
		stateB[b] = 0.0;
		stateC[b] = 0.0;
	}
	
	
	// SB Transform on the original bg fields
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
			for (int imu = -1; imu <= 1; imu += 2) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
				int ii = (int)((i - iMin)*DIrecip);
				for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
					for (int jmu = -1; jmu <= 1; jmu += 2) {
						real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5));
						int jj = (int)((j - jMin)*DJrecip);
						for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
							if ((iNode < 0) or (iNode >= iDim)) continue;
							for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
								if ((jNode < 0) or (jNode >= jDim)) continue;
								int iBCL, jBCL;
								if (var == 0) {
									iBCL = R2T20;
									jBCL = R1T2;
								} else if (var == 1) {									
									iBCL = R2T20;
									jBCL = R2T20;
								} else {
									iBCL = R1T2;
									jBCL = R1T2;
								}
								if (((var == 0) and iNode) or ((var == 1) and iNode and jNode) or (var > 1)) {
									real im = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL, R1T2);
									real jm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL, R1T2);
									int bgJ = jIndex*2 + (jmu+1)/2;
									int bgI = iIndex*2 + (imu+1)/2;
									
									// Reset Psi to zero for multiple iterations!!
									if (var == 1) {
										bgFields[varDim*(iDim-1)*2*bgJ +varDim*bgI + var] = 0;
									}
									
									stateB[varDim*iDim*jNode +varDim*iNode + var] += 
										0.25 * bgFields[varDim*(iDim-1)*2*bgJ +varDim*bgI + var] * im * jm;
								}	
							}
						}	
					}
				}
			}
		}	
	}
	
	// Compute and display the variable RMS and BG errors for reference
	for (int var = 0; var < varDim; var++) {
		double varScale = 0;
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				varScale += stateB[varDim*iDim*jIndex +varDim*iIndex + var] * stateB[varDim*iDim*jIndex +varDim*iIndex + var];
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
	
	// Estimate the mixing lengths for the variables here?
	
	// SA transform = bg B's -> bg A's
	SAtransform(stateB, bgState);
	
	// Load the obs locally and weight the nonlinear observation operators by interpolated bg fields
	for (int m = 0; m < mObs; m++) {
		int mi = m*11;
		for (int ob = 0; ob < 11; ob++) {
			obsVector[mi+ob] = rawObs[mi+ob];
		}
		real type = obsVector[mi+10];
		if (type <= 1) continue; 

		real i = obsVector[mi+8];
		real invI;
		if (i != 0) {
			invI = 1./i;
		} else {
			invI = 0.;
		}
		real j = obsVector[mi+9];
		real vBG = 0.;
		real uBG = 0.;
		real rhoprime = 0.;
		real qvprime = 0.;
		
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		real ibasis = 0.;
		real jbasis = 0.;
		real idbasis = 0.;
		real jdbasis = 0.;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {				
				if ((iNode < 0) or (iNode >= iDim) or (jNode < 0) or (jNode >= jDim)) continue;
				if (iNode) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					vBG += bgState[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * invI * 1.e3;
				}
				if (iNode and jNode) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R2T20, R1T2);
					idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, R2T20, R1T2);
					jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, R2T20, R1T2);
					float coeff = bgState[varDim*iDim*jNode +varDim*iNode + 1];
					uBG += coeff * ibasis * (-jdbasis) * invI * 1.e5;
				}
				ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
				jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
				qvprime += bgState[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis;
				rhoprime += bgState[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis;
			}
		}
		
		real rhoBar = rhoBase*exp(-rhoInvScaleHeight*j);
		real qBar = 19.562 - 0.004066*j + 7.8168e-7*j*j;
		real rhoa = rhoBar + rhoprime / 100;
		real rhoq = (qBar + qvprime) * rhoa / 1000.;
		real rhoBG = rhoa + rhoq;
		
		// Adjust the relevant fields
		if (type == MetObs::sfmr) {
			//real wsBG = vBG*vBG; // + uBG*uBG;
			obsVector[mi] *= rhoBG;
			//obsVector[mi+2] = 2.*vBG/wsBG;
			//obsVector[mi+3] = 0.; //2.*uBG/wsBG;
		}
		if (type == MetObs::radar) {
			obsVector[mi] *= rhoBG;
		}
	}
	
	// d = y - HXb
	calcInnovation();

	// Output the original background field
	outputAnalysis("background", bgState , false);

	// HTd
	calcHTranspose(innovation, stateC);
	
	SCtranspose(stateC, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtransform(stateA, stateB);
	
	SBtranspose(stateB, CTHTd);
			
}	

double CostFunctionXY_CPU::funcValue(double* state)
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
	//#pragma omp parallel for reduction(+:obIP)
	for (int m = 0; m < mObs; m++) {
		obIP += (HCq[m]-innovation[m])*(obsVector[m*11+1])*(HCq[m]-innovation[m]);
	}
	double J = 0.5*qIP + 0.5*obIP;
	return J;
	
}

void CostFunctionXY_CPU::funcGradient(double* state, double* gradient)
{
	
	updateHCq(state);
	
	// HTHCq
	calcHTranspose(HCq, stateC);
	
	SCtranspose(stateC, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtransform(stateA, stateB);
	
	SBtranspose(stateB, stateU);
	
	for (int n = 0; n < nState; n++) {
		gradient[n] = state[n] + stateU[n] - CTHTd[n];
	}
	
	
}

void CostFunctionXY_CPU::updateHCq(double* state)
{

	// SB transform from the q's
	SBtransform(state, stateB);
	
	// S (SA transform) yield A's, put it in the tempState
	SAtransform(stateB, stateA);
	
	SCtransform(stateA, stateC);
	
	// H
	#pragma omp parallel for
	for (int m = 0; m < mObs; m++) {
		int mi = m*11;
		real w1 = obsVector[mi+2];
		real w2 = obsVector[mi+3];
		real w3 = obsVector[mi+4];
		real w4 = obsVector[mi+5];
		real w5 = obsVector[mi+6];
		real w6 = obsVector[mi+7];		
		real i = obsVector[mi+8];
		real invI;
		if (i != 0) {
			invI = 1./i;
		} else {
			invI = 0.;
		}
		real j = obsVector[mi+9];
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
					if (iNode) {
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
						tempsum += stateC[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * w1 * invI * 1.e3;
					}
				}
				if(w2 or w3) {
					if (iNode and jNode) {
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R2T20, R1T2);
						idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, R2T20, R1T2);
						jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, R2T20, R1T2);
						float coeff = stateC[varDim*iDim*jNode +varDim*iNode + 1];
						tempsum += coeff * ibasis * (-jdbasis) * w2 * invI * 1.e5;
						tempsum += coeff * idbasis * jbasis * w3 * invI * 1.e2;
					}
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					tempsum += stateC[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis * w4;
					tempsum += stateC[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * w5;
					tempsum += stateC[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * w6;
				}
			}
		}
		HCq[m] = tempsum;
		//cout << m << "\t" << HCq[m] << endl;
	}	
	//cout << endl;
	
}

void CostFunctionXY_CPU::updateBG()
{

	// SB transform from the q's
	SBtransform(currState, stateB);
	
	// S (SA transform) yield A's
	SAtransform(stateB, stateA);
	
	SCtransform(stateA, stateC);
	
	outputAnalysis("increment", stateC, false);	
	
	// In BG update we are directly summing C + A
	ofstream cstream("CoeffAnalysis.out");
	cstream << "Variable\tI\tJ\tBackground\tAnalysis\tIncrement\n";
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				cstream << var << "\t" << iIndex << "\t" << jIndex << "\t"; 
				cstream << bgState[varDim*iDim*jIndex +varDim*iIndex + var] << "\t";
				bgState[varDim*iDim*jIndex +varDim*iIndex + var] += stateC[varDim*iDim*jIndex +varDim*iIndex + var];
				cstream << bgState[varDim*iDim*jIndex +varDim*iIndex + var] << "\t"; 
				cstream << stateC[varDim*iDim*jIndex +varDim*iIndex + var] << endl;
			}
		}
	}
	
	// S (SA transform) yield A's
	//SAtransform(bgState, stateA);
	
	outputAnalysis("analysis", bgState, true);
	
}

void CostFunctionXY_CPU::calcInnovation()
{
	// Initialize and fill the innovation vector
	cout << "Initializing innovation vector..." << endl;
	for (int m = 0; m < mObs; m++) {
		HCq[m] = 0.0;
		innovation[m] = obsVector[m*11];
	}
	
	real innovationRMS = 0;
	#pragma omp parallel for reduction(+:innovationRMS)
	for (int m = 0; m < mObs; m++) {
		int mi = m*11;
		real w1 = obsVector[mi+2];
		real w2 = obsVector[mi+3];
		real w3 = obsVector[mi+4];
		real w4 = obsVector[mi+5];
		real w5 = obsVector[mi+6];
		real w6 = obsVector[mi+7];		
		real i = obsVector[mi+8];
		real invI;
		if (i != 0) {
			invI = 1./i;
		} else {
			invI = 0.;
		}
		real j = obsVector[mi+9];
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
					if (iNode) {
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
						tempsum += bgState[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * w1 * invI * 1.e3;
					}
				}
				if(w2 or w3) {
					if (iNode and jNode) {
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R2T20, R1T2);
						idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, R2T20, R1T2);
						jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, R2T20, R1T2);
						float coeff = bgState[varDim*iDim*jNode +varDim*iNode + 1];
						tempsum += coeff * ibasis * (-jdbasis) * w2 * invI * 1.e5;
						tempsum += coeff * idbasis * jbasis * w3 * invI * 1.e2;
					}
				}
				if (w4 or w5 or w6) {
					//ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 2);
					//jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					tempsum += bgState[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis * w4;
					tempsum += bgState[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * w5;
					tempsum += bgState[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * w6;
				}
			}
		}
		innovation[m] -= tempsum;
		innovationRMS += (innovation[m]*innovation[m]);
	}
	if (mObs) innovationRMS /= mObs;
	innovationRMS = sqrt(innovationRMS);
	cout << "Innovation RMS : " << innovationRMS << endl;
	
}	

void CostFunctionXY_CPU::calcHTranspose(const real* yhat, real* Astate) 
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
				int mi = m*11;
				real qhat = yhat[m];
				real invError = obsVector[mi+1];
				real w1 = obsVector[mi+2];
				real w2 = obsVector[mi+3];
				real w3 = obsVector[mi+4];
				real w4 = obsVector[mi+5];
				real w5 = obsVector[mi+6];
				real w6 = obsVector[mi+7];		
				real i = obsVector[mi+8];
				real invI;
				if (i != 0) {
					invI = 1./i;
				} else {
					invI = 0.;
				}
				real j = obsVector[mi+9];
				int ii = (int)((i - iMin)*DIrecip);
				if ((iIndex < ii-1) or (iIndex > ii+2)) continue;
				int jj = (int)((j - jMin)*DJrecip);
				if ((jIndex < jj-1) or (jIndex > jj+2)) continue;
				real ibasis = 0;
				real jbasis = 0;
				real idbasis = 0;
				real jdbasis = 0;
				
				if(w1) {
					if (iIndex) {
						ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
						jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
						#pragma omp atomic
						Astate[varDim*iDim*jIndex +varDim*iIndex] 
							+= qhat * ibasis * jbasis * w1 * invError * invI * 1.e3;
						//cout << m << "\t" << Astate[varDim*iDim*jIndex +varDim*iIndex] << endl;
					}
				}
				if(w2 or w3) {
					if (iIndex and jIndex) {
						ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
						jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 0, R2T20, R1T2);
						idbasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 1, R2T20, R1T2);
						jdbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 1, R2T20, R1T2);
						#pragma omp atomic
						Astate[varDim*iDim*jIndex +varDim*iIndex + 1] 
							+= qhat * ibasis * (-jdbasis) * w2 * invError * invI* 1.e5
							 + qhat * idbasis * jbasis * w3 * invError * invI * 1.e2;
					}
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					#pragma omp atomic
					Astate[varDim*iDim*jIndex +varDim*iIndex + 2] 
						+= qhat * ibasis * jbasis * w4 * invError;
					#pragma omp atomic
					Astate[varDim*iDim*jIndex +varDim*iIndex + 3] 
						+= qhat * ibasis * jbasis * w5 * invError;
					#pragma omp atomic
					Astate[varDim*iDim*jIndex +varDim*iIndex + 4]
						+= qhat * ibasis * jbasis * w6 * invError;
				}
			}
		}
	}
	
}

bool CostFunctionXY_CPU::SAtransform(real* Bstate, real* Astate)
{
	int k;
	for (int var = 0; var < varDim; var++) {
		real* jB = new real[jDim];
		real* x = new real[jDim];
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jB[jIndex] = Bstate[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			// Solve for A's using compact storage
			real sum = 0;
			for (int j = 0; j < jDim; j++) {
				for (sum=jB[j], k=-1;k>=-3;k--) {
					if ((j+k >= 0) and ((j*4-k) >= 0))
						sum -= jL[jDim*4*var + j*4-k]*x[j+k];
				}
				x[j] = sum/jL[jDim*4*var + j*4];
			}	
			for (int j=jDim-1;j>=0;j--) {
				for (sum=x[j], k=1;k<=3;k++) {
					if ((j+k < jDim) and (((j+k)*4+k) < jDim*4))
						sum -= jL[jDim*4*var + (j+k)*4+k]*x[j+k];
				}
				x[j] = sum/jL[jDim*4*var + j*4];
			}
			
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = x[jIndex]; 
			}
		}
		delete[] jB;
		delete[] x;
		
		real* iB = new real[iDim];
		x = new real[iDim];
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iB[iIndex] = Astate[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			// Solve for A's using compact storage
			real sum = 0;
			for (int i = 0; i < iDim; i++) {
				for (sum=iB[i], k=-1;k>=-3;k--) {
					if ((i+k >= 0) and ((i*4-k) >= 0))
						sum -= iL[iDim*4*var + i*4-k]*x[i+k];
				}
				x[i] = sum/iL[iDim*4*var + i*4];
			}	
			for (int i=iDim-1;i>=0;i--) {
				for (sum=x[i], k=1;k<=3;k++) {
					if ((i+k < iDim) and (((i+k)*4+k) < iDim*4))
						sum -= iL[iDim*4*var + (i+k)*4+k]*x[i+k];
				}
				x[i] = sum/iL[iDim*4*var + i*4];
			}
			
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = x[iIndex]; 
			}
			
		}
		delete[] iB;
		delete[] x;
	}

	return true;
}		

void CostFunctionXY_CPU::SBtransform(real* Ustate, real* Bstate)
{
	// Clear the Bstate
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				Bstate[varDim*iDim*jIndex +varDim*iIndex + var] = 0.;
			}
		}
	}

	/* Clear the Uprime
	for (int n = 0; n < nState; n++) {
		Uprime[n] = 0.;
	}
	int iSpread = int(LI*2/DI);
	int jSpread = int(LJ*2/DJ);

	#pragma omp parallel for
	for (int var = 0; var < varDim; var++) {
		real errorscale = 1;
		for (int iIndex = 0; iIndex < (iDim-1)*2; iIndex++) {
			real i1 = iMin + DI * (iIndex/2 + (0.5*sqrt(1./3.) * pow(-1.,(iIndex % 2)+1.) + 0.5));
			
			for (int jIndex = 0; jIndex < (jDim-1)*2; jIndex++) {
				real j1 = jMin + DJ * (jIndex/2 + (0.5*sqrt(1./3.) * pow(-1.,(jIndex % 2)+1.) + 0.5));
				
				for (int iPrime = iIndex-iSpread; iPrime <= iIndex+iSpread; iPrime++) {
					if ((iPrime < 0) or (iPrime >= (iDim-1)*2)) continue;

					real i2 = iMin + DI * (iPrime/2 + (0.5*sqrt(1./3.) * pow(-1.,(iPrime % 2)+1.) + 0.5));
					for (int jPrime = jIndex-jSpread; jPrime <= jIndex+jSpread; jPrime++) {
						if ((jPrime < 0) or (jPrime >= (jDim-1)*2)) continue;
						real j2 = jMin + DJ * (jPrime/2 + (0.5*sqrt(1./3.) * pow(-1.,(jPrime % 2)+1.) + 0.5));

						real f1 = bgFields[varDim*(iDim-1)*2*jIndex +varDim*iIndex + var];
						real f2 = bgFields[varDim*(iDim-1)*2*jPrime +varDim*iPrime + var];
						if (var <= 1) {
							// Scale the BG error by radius
							errorscale = i1;
						} else {
							errorscale = 1.;
						}
						real correlation = exp( -0.5*((i1-i2)*(i1-i2)/(LI*LI) 
													  + (j1-j2)*(j1-j2)/(LJ*LJ) 
													  + (f1-f2)*(f1-f2)/(LF[var]*LF[var])) );
						Uprime[varDim*(iDim-1)*2*jIndex +varDim*iIndex + var] += 
						Ustate[varDim*(iDim-1)*2*jPrime +varDim*iPrime + var] * correlation * bgError[var] * errorscale;
					}
				}
			}
		}
	} */
					
	//#pragma omp parallel for
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
			for (int imu = -1; imu <= 1; imu += 2) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
				int ii = (int)((i - iMin)*DIrecip);
				for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
					for (int jmu = -1; jmu <= 1; jmu += 2) {
						real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5));
						int jj = (int)((j - jMin)*DJrecip);
						for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
							if ((iNode < 0) or (iNode >= iDim)) continue;
							for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
								if ((jNode < 0) or (jNode >= jDim)) continue;
								int iBCL, jBCL;
								if (var == 0) {
									iBCL = R2T20;
									jBCL = R1T2;
								} else if (var == 1) {									
									iBCL = R2T20;
									jBCL = R2T20;
								} else {
									iBCL = R1T2;
									jBCL = R1T2;
								}
								if (((var == 0) and iNode) or ((var == 1) and iNode and jNode) or (var > 1)) {
									real im = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL, R1T2);
									real jm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL, R1T2);
									int bgJ = jIndex*2 + (jmu+1)/2;
									int bgI = iIndex*2 + (imu+1)/2;
									Bstate[varDim*iDim*jNode +varDim*iNode + var] += 
										//0.25 * Uprime[varDim*(iDim-1)*2*bgJ +varDim*bgI + var] * im * jm; 
										0.25 * Ustate[varDim*(iDim-1)*2*bgJ +varDim*bgI + var] * im * jm; 
								}
							}
						}	
					}
				}
			}
		}	
	}
		
}


void CostFunctionXY_CPU::SBtranspose(real* Bstate, real* Ustate)
{
		
	// Clear the Ustate
	for (int n = 0; n < nState; n++) {
		Ustate[n] = 0;
		Uprime[n] = 0;
	}
		
	//#pragma omp parallel for
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
			for (int imu = -1; imu <= 1; imu += 2) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
				int ii = (int)((i - iMin)*DIrecip);
				for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
					for (int jmu = -1; jmu <= 1; jmu += 2) {
						real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5));
						int jj = (int)((j - jMin)*DJrecip);
						for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
							if ((iNode < 0) or (iNode >= iDim)) continue;
							for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
								if ((jNode < 0) or (jNode >= jDim)) continue;
								int iBCL, jBCL;
								if (var == 0) {
									iBCL = R2T20;
									jBCL = R1T2;
								} else if (var == 1) {									
									iBCL = R2T20;
									jBCL = R2T20;
								} else {
									iBCL = R1T2;
									jBCL = R1T2;
								}
								if (((var == 0) and iNode) or ((var == 1) and iNode and jNode) or (var > 1)) {
									real im = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL, R1T2);
									real jm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL, R1T2);
									int bgJ = jIndex*2 + (jmu+1)/2;
									int bgI = iIndex*2 + (imu+1)/2;
									Ustate[varDim*(iDim-1)*2*bgJ +varDim*bgI + var] += 
									//Uprime[varDim*(iDim-1)*2*bgJ +varDim*bgI + var] += 
										0.25 * Bstate[varDim*iDim*jNode +varDim*iNode + var] * im * jm;
								}
							}
						}	
					}
				}
			}
		}	
	}
	
	/*
	int iSpread = int(LI*2/DI);
	int jSpread = int(LJ*2/DJ);
	#pragma omp parallel for
	for (int var = 0; var < varDim; var++) {
		real errorscale = 1.;
		for (int iIndex = 0; iIndex < (iDim-1)*2; iIndex++) {
			real i1 = iMin + DI * (iIndex/2 + (0.5*sqrt(1./3.) * pow(-1.,(iIndex % 2)+1.) + 0.5));
			
			for (int jIndex = 0; jIndex < (jDim-1)*2; jIndex++) {
				real j1 = jMin + DJ * (jIndex/2 + (0.5*sqrt(1./3.) * pow(-1.,(jIndex % 2)+1.) + 0.5));

				for (int iPrime = iIndex-iSpread; iPrime <= iIndex+iSpread; iPrime++) {
					if ((iPrime < 0) or (iPrime >= (iDim-1)*2)) continue;
					real i2 = iMin + DI * (iPrime/2 + (0.5*sqrt(1./3.) * pow(-1.,(iPrime % 2)+1.) + 0.5));

					for (int jPrime = jIndex-jSpread; jPrime <= jIndex+jSpread; jPrime++) {
						if ((jPrime < 0) or (jPrime >= (jDim-1)*2)) continue;
						real j2 = jMin + DJ * (jPrime/2 + (0.5*sqrt(1./3.) * pow(-1.,(jPrime % 2)+1.) + 0.5));						

						real f1 = bgFields[varDim*(iDim-1)*2*jIndex +varDim*iIndex + var];
						real f2 = bgFields[varDim*(iDim-1)*2*jPrime +varDim*iPrime + var];
						if (var <= 1) {
							// Scale the BG error by radius
							errorscale = i1;
						} else {
							errorscale = 1.;
						}
						real correlation = exp( -0.5*((i1-i2)*(i1-i2)/(LI*LI) 
													  + (j1-j2)*(j1-j2)/(LJ*LJ) 
													  + (f1-f2)*(f1-f2)/(LF[var]*LF[var])) );
						Ustate[varDim*(iDim-1)*2*jIndex +varDim*iIndex + var] += 
						Uprime[varDim*(iDim-1)*2*jPrime +varDim*iPrime + var] * correlation * bgError[var] * errorscale;
					}
				}
			}
		}
	} */
	
}

void CostFunctionXY_CPU::SCtransform(real* Astate, real* Cstate)
{
	// Isoptropic Recursive filter for speed, no anisotropic "triad" working yet 
	for (int var = 0; var < varDim; var++) {
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = Astate[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
		//FI
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iTemp[iIndex] = Astate[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			iFilter->filterArray(iTemp, iDim);
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				// D
				real errorscale = 1.;
				if (var <= 1) {
					// Scale the BG error by radius
					errorscale = iIndex * DI * bgErrorScale;
				} else {
					errorscale = bgErrorScale;
				}
				//D
				Cstate[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex] * bgError[var] * errorscale; 
			}
		}
	}
}

void CostFunctionXY_CPU::SCtranspose(real* Cstate, real* Astate)
{
	
	// Isoptropic Recursive filter for speed, no anisotropic "triad" working yet 
	for (int var = 0; var < varDim; var++) {
		//FI & D
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				real errorscale = 1.;
				if (var <= 1) {
					// Scale the BG error by radius
					errorscale = iIndex * DI * bgErrorScale;
				} else {
					errorscale = bgErrorScale;
				}								
				iTemp[iIndex] = Cstate[varDim*iDim*jIndex +varDim*iIndex + var] * bgError[var] * errorscale;
			}
			iFilter->filterArray(iTemp, iDim);
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				Cstate[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex]; 
			}
		}
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = Cstate[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
	}
}


bool CostFunctionXY_CPU::setupSplines()
{
	
	// Do the spline via a Cholesky decomposition
	// and manipulate the DC Filter
	real Pi = 3.141592653589793;
	real** P = new real*[iDim];
	real* p = new real[iDim];
	iL = new real[varDim*iDim*4];
	real* iBL = new real[varDim*iDim*4];
	for (int i = 0; i < iDim; i++) {
		P[i] = new real[iDim];
		p[i] = 0.;
	}
		
	for (int i = 0; i < varDim*iDim*4; i++) {
		iL[i] = 0;
		iBL[i] = 0;
	}
	
	real cutoff_wl = 4;
	real eq = pow( (cutoff_wl/(2*Pi)) , 6);
	for (int var = 0; var < varDim; var++) {
		int iBCL;
		if (var <= 1) {
			iBCL = R2T20;
		} else {
			iBCL = R1T2;
		}
				
		for (int i = 0; i < iDim; i++) {
			for (int j = 0; j < iDim; j++) {
				P[i][j] = 0;
			}
		}
				
		for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
			for (int imu = -1; imu <= 1; imu += 2) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
				int ii = (int)((i - iMin)*DIrecip);
				for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
					if ((iNode < 0) or (iNode >= iDim)) continue;
					real pm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL, R1T2);
					real qm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 3, iBCL, R1T2);
					real pn, qn;
					P[iNode][iNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((iNode+1) < iDim) {
						pn = Basis(iNode+1, i, iDim-1, iMin, DI, DIrecip, 0, iBCL, R1T2);
						qn = Basis(iNode+1, i, iDim-1, iMin, DI, DIrecip, 3, iBCL, R1T2);
						P[iNode][iNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+1][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((iNode+2) < iDim) {
						pn = Basis(iNode+2, i, iDim-1, iMin, DI, DIrecip, 0, iBCL, R1T2);
						qn = Basis(iNode+2, i, iDim-1, iMin, DI, DIrecip, 3, iBCL, R1T2);
						P[iNode][iNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+2][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((iNode+3) < iDim) {
						pn = Basis(iNode+3, i, iDim-1, iMin, DI, DIrecip, 0, iBCL, R1T2);
						qn = Basis(iNode+3, i, iDim-1, iMin, DI, DIrecip, 3, iBCL, R1T2);
						P[iNode][iNode+3] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+3][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
				}
			}
		}
		
		/* for (int i = 0; i < iDim; i++) {
			for (int j = 0; j < iDim; j++) {
				cout << P[i][j] << " ";
			} cout << endl;
		} */	
				
		// Cholesky decomp of P+Q
		for (int i=0;i<iDim;i++) {
			for (int j=i;j<iDim;j++) {
				double sum=P[i][j];
				for (int k=i-1;k>=0;k--) {
					sum -= P[i][k]*P[j][k];
				}
				if (i == j) {
					if (sum <= 0.0) { 
						std::cout << "cholesky failed at i,j sum\n";
						break;
					} else {
						p[i] = sqrt(sum);
					}
				} else {
					P[j][i]=sum/p[i];
					if (p[i] == 0.) { 
						std::cout << "Problem! " << i << "\t" << j << "\n";
					}
				}
			}
		}
		
		for (int i = 0; i < iDim; i++) {
			iL[iDim*4*var + i*4] = p[i];
			//cout << iL[iDim*4*var + i*4] << " ";
			for (int n=1;n<4;n++) {
				if ((i-n) >= 0) {
					iL[iDim*4*var + i*4+n] = P[i][i-n];
				}
				//cout << iL[iDim*4*var + i*4 + n] << " ";
			} //cout << endl;
		} //cout << endl;
	}
	
	// Clear the temporary memory and reallocate for jDim
	for (int i = 0; i < iDim; i++) {
		delete[] P[i];
	}
	delete[] P;
	delete[] p;
	
	P = new real*[jDim];
	p = new real[jDim];
	jL = new real[varDim*jDim*4];
	for (int j = 0; j < jDim; j++) {
		P[j] = new real[jDim];
		p[j] = 0.;
	}
	
	for (int j = 0; j < varDim*jDim*4; j++) {
		jL[j] = 0;
	}
	
	cutoff_wl = 4;
	eq = pow( (cutoff_wl/(2*Pi)) , 6);
	for (int var = 0; var < varDim; var++) {
		int jBCL;
		if (var == 1) {
			jBCL = R2T20;
		} else {
			jBCL = R1T2;
		}

		for (int i = 0; i < jDim; i++) {
			for (int j = 0; j < jDim; j++) {
				P[i][j] = 0;
			}
		}		
		for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
			for (int jmu = -1; jmu <= 1; jmu += 2) {
				real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5));
				int jj = (int)((j - jMin)*DJrecip);
				for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
					if ((jNode < 0) or (jNode >= jDim)) continue;
					real pm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL, R1T2);
					real qm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 3, jBCL, R1T2);
					real pn, qn;
					P[jNode][jNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((jNode+1) < jDim) {
						pn = Basis(jNode+1, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL, R1T2);
						qn = Basis(jNode+1, j, jDim-1, jMin, DJ, DJrecip, 3, jBCL, R1T2);
						P[jNode][jNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+1][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((jNode+2) < jDim) {
						pn = Basis(jNode+2, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL, R1T2);
						qn = Basis(jNode+2, j, jDim-1, jMin, DJ, DJrecip, 3, jBCL, R1T2);
						P[jNode][jNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+2][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));						
					}
					if ((jNode+3) < jDim) {
						pn = Basis(jNode+3, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL, R1T2);
						qn = Basis(jNode+3, j, jDim-1, jMin, DJ, DJrecip, 3, jBCL, R1T2);
						P[jNode][jNode+3] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+3][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));						
					}
				}
			}
		}
		
		// Cholesky decomp
		for (int i=0;i<jDim;i++) {
			for (int j=i;j<jDim;j++) {
				double sum=P[i][j];
				for (int k=i-1;k>=0;k--) {
					sum -= P[i][k]*P[j][k];
				}
				if (i == j) {
					if (sum <= 0.0) { 
						std::cout << "cholesky failed at i,j sum\n";
						break;
					} else {
						p[i] = sqrt(sum);
					}
				} else {
					P[j][i]=sum/p[i];
					if (p[i] == 0.) { 
						std::cout << "Problem! " << i << "\t" << j << "\n";
					}
				}
			}
		}
		
		for (int j = 0; j < jDim; j++) {
			jL[jDim*4*var + j*4] = p[j];
			//cout << jL[jDim*4*var + j*4] << " ";
			for (int n=1;n<4;n++) {
				if ((j-n) >= 0) {
					jL[jDim*4*var + j*4+n] = P[j][j-n];
				}
				//cout << jL[jDim*4*var + j*4 + n] << " ";
			} //cout << endl;
		} //cout << endl;
	}
	
	// Clear the temporary memory
	for (int j = 0; j < jDim; j++) {
		delete[] P[j];
	}
	delete[] P;
	delete[] p;
	
	/*p = new real[iDim];
	//Analytic Test set
	for (int i = 0; i < iDim; i++) {
		p[i] = 0.;
	}
	int var = 3;
	for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
		for (int imu = -1; imu <= 1; imu += 2) {
			real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
			int ii = (int)((i - iMin)*DIrecip);
			for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
				if ((iNode < 0) or (iNode >= iDim)) continue;
				real pm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
				p[iNode] += 0.5*pm*exp(-(i-(iDim-5))*(i-(iDim-5))/25);
			}
		}
	}
	real* x = new real[iDim];
	real* y = new real[iDim];
	// Solve for A's
	real sum = 0;
	int k; */
	
	/* for (int i = 0; i < iDim; i++) {
		for (sum=p[i], k=i-1;k>=0;k--)
			sum -= P[i][k]*x[k];
		x[i] = sum/iL[iDim*4*var + i*4];
	}	
	for (int i=iDim-1;i>=0;i--) {
		for (sum=x[i], k=i+1;k<iDim;k++)
			sum -= P[k][i]*x[k];
		x[i] = sum/iL[iDim*4*var + i*4];
	} */
	 
	/*Solve for A's using compact storage
	sum = 0;
	for (int i = 0; i < iDim; i++) {
		for (sum=p[i], k=-1;k>=-3;k--) {
			if ((i+k >= 0) and ((i*4-k) >= 0))
				sum -= iL[iDim*4*var + i*4-k]*y[i+k];
		}
		y[i] = sum/iL[iDim*4*var + i*4];
	}	
	for (int i=iDim-1;i>=0;i--) {
		for (sum=y[i], k=1;k<=3;k++) {
			if ((i+k < iDim) and (((i+k)*4+k) < iDim*4))
				sum -= iL[iDim*4*var + (i+k)*4+k]*y[i+k];
		}
		y[i] = sum/iL[iDim*4*var + i*4];
	}
	//for (int i = 0; i < iDim; i++) {
	//	cout << i << "\t" << x[i] << "\t" << y[i] << endl;
	//} cout << endl;
	
	real si;
	for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
		for (real imu = 0; imu <= 1; imu+=2) {
//real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
			real i = iMin + DI * (imu + iIndex);
			int ii = (int)((i - iMin)*DIrecip);
			si = 0;
			for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
				if ((iNode < 0) or (iNode >= iDim)) continue;
				real pm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
				if (iNode >= 0)
					si += y[iNode]*pm; // + sin(2*i*2*Pi/iDim) + sin(4*i*2*Pi/iDim) + sin(8*i*2*Pi/iDim) + sin(16*i*2*Pi/iDim);
			}
			//real analytic = sin(i*2*Pi/iDim) + sin(2*i*2*Pi/iDim) + sin(4*i*2*Pi/iDim) + sin(8*i*2*Pi/iDim) + sin(16*i*2*Pi/iDim); 
			real analytic = exp(-(i-(iDim-5))*(i-(iDim-5))/25);
			cout << i << "\t" << p[iIndex] << "\t" << y[iIndex] << "\t" << si << "\t" << analytic << "\t" << si - analytic << endl;
		}
	} */
	
	return true;
	
}

bool CostFunctionXY_CPU::outputAnalysis(const QString& suffix, real* Astate, bool updateMish)
{
	
	real* fieldNodes = new real[12*iDim*jDim];
	
	// H --> to Mish for output
	QString fluxout = "Flux_" + suffix + ".out";
	ofstream fluxstream(fluxout.toAscii().data());
	fluxstream << "Radius\tHeight\trhoM\trhoE\tu\tv\tw\tPsi\tqv\trho\tT\tP\n";
	fluxstream.precision(10);
	real CoriolisF = 6e-5;
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		for (int ihalf = 0; ihalf <=1; ihalf++) {
			for (int imu = -ihalf; imu <= ihalf; imu++) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5*ihalf));
				real invI;
				if (i != 0) {
					invI = 1./i;
				} else {
					invI = 0.;
				}
				if (i > (iDim-1)*DI) continue;
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					for (int jhalf =0; jhalf <=1; jhalf++) {
						for (int jmu = -jhalf; jmu <= jhalf; jmu++) {
							real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5*jhalf));
							real rhoBar = rhoBase*exp(-rhoInvScaleHeight*j);
							real qBar = 19.562 - 0.004066*j + 7.8168e-7*j*j;
							real hBar = 3.5e5;
							if (j > (jDim-1)*DJ) continue;	
							int ii = (int)((i - iMin)*DIrecip);
							int jj = (int)((j - jMin)*DJrecip);
							real ibasis = 0.;
							real jbasis = 0.;
							real idbasis = 0.;
							real jdbasis = 0.;
							real rhov = 0.;
							real rhou = 0.;
							real rhow = 0.;
							real hprime = 0.;
							real qvprime = 0.;
							real rhoprime = 0.;
							real psi = 0.;
							
							for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
								for (int jNode = jj-1; jNode <= jj+2; ++jNode) {				
									if ((iNode < 0) or (iNode >= iDim) or (jNode < 0) or (jNode >= jDim)) continue;
									if (iNode) {
										ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
										jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
										rhov += Astate[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * invI * 1.e3;
									}
									
									if (iNode and jNode) {
										ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
										jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R2T20, R1T2);
										idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, R2T20, R1T2);
										jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, R2T20, R1T2);
										float coeff = Astate[varDim*iDim*jNode +varDim*iNode + 1];
										rhou += coeff * ibasis * (-jdbasis) * invI * 1.e5;
										rhow += coeff * idbasis * jbasis * invI * 1.e2;
										psi += coeff * ibasis * jbasis;
									}
									
									ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
									jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
									hprime += Astate[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis;
									qvprime += Astate[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis;
									rhoprime += Astate[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis;
								}
							}
							
							// Get it to relevant variables for the flux calculation
							real rhoa = rhoBar + rhoprime / 100;
							real rhoq = (qBar + qvprime) * rhoa / 1000.;
							real rho = rhoa + rhoq;
							real v = rhov / rho;
							real u = rhou / rho;
							real w = rhow / rho;
							real KE = 0.5*rho*(v*v + u*u + w*w);
							real rhoM = (i*rhov*1000 + 0.5*rho*CoriolisF*i*i*1e6);
							real rhoE = rho*(hBar + hprime*1.e3) + KE;
							real T = ((hBar + hprime*1.e3) - 2.501e3*(qBar + qvprime) - 9.81*j)/1005.7;
							real press = T*rhoa*287./100.;
							//real e = rhoE/rho;
							
							// Output it
							fluxstream << scientific << i << "\t" << j << "\t" << rhoM << "\t" << rhoE
							<< "\t" << u << "\t" << v << "\t" << w << "\t" << psi 
							<< "\t" << qvprime << "\t" << rhoprime << "\t" << T << "\t" << press << "\n";
							
							if (updateMish and imu and jmu and ihalf and jhalf) {
								int bgJ = jIndex*2 + (jmu+1)/2;
								int bgI = iIndex*2 + (imu+1)/2;
								bgFields[varDim*(iDim-1)*2*bgJ +varDim*bgI + 0] = rhov*i*1.e-3;
								bgFields[varDim*(iDim-1)*2*bgJ +varDim*bgI + 1] = psi;
								bgFields[varDim*(iDim-1)*2*bgJ +varDim*bgI + 2] = hprime;
								bgFields[varDim*(iDim-1)*2*bgJ +varDim*bgI + 3] = qvprime;
								bgFields[varDim*(iDim-1)*2*bgJ +varDim*bgI + 4] = rhoprime;
							}
								
							if (!ihalf and !jhalf){
								// On the node
								fieldNodes[12*iDim*jIndex + 12*iIndex + 0] = rhov*i*1.e-3;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 1] = hprime;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 2] = u;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 3] = v;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 4] = w;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 5] = psi;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 6] = qvprime;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 7] = rhoprime;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 8] = T;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 9] = press;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 10] = rhoM;
								fieldNodes[12*iDim*jIndex + 12*iIndex + 11] = rhoE;
							}
							
						}
					}
				}
			}
		}
	}
	
	
	// H -> to QC
	// Write the Obs to a summary text file
	QString qcout = "QC_" + suffix + ".out";
	ofstream qcstream(qcout.toAscii().data());
	ifstream obstream("./Observations.out");
	ostream_iterator<string> os(qcstream, "\t ");
	*os++ = "Observation";
	*os++ = "Inverse Error";
	*os++ = "Weight 1";
	*os++ = "Weight 2";
	*os++ = "Weight 3";
	*os++ = "Weight 4";
	*os++ = "Weight 5";
	*os++ = "Weight 6";
	*os++ = "r";
	*os++ = "z";	
	*os++ = "Type";
	*os++ = "Analysis";
	*os++ = "Background";
	qcstream << endl;
	double temp;
	
	ostream_iterator<double> od(qcstream, "\t ");
	for (int m = 0; m < mObs; m++) {
		int mi = m*11;
		real w1 = obsVector[mi+2];
		real w2 = obsVector[mi+3];
		real w3 = obsVector[mi+4];
		real w4 = obsVector[mi+5];
		real w5 = obsVector[mi+6];
		real w6 = obsVector[mi+7];		
		real i = obsVector[mi+8];
		real invI;
		if (i != 0) {
			invI = 1./i;
		} else {
			invI = 0.;
		}
		real j = obsVector[mi+9];
		real type = obsVector[mi+10];
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
					if (iNode) {
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
						tempsum += Astate[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * w1 * invI * 1.e3;
					}
				}
				if(w2 or w3) {
					if (iNode and jNode) {
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R2T20, R1T2);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R2T20, R1T2);
						idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, R2T20, R1T2);
						jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, R2T20, R1T2);
						float coeff = Astate[varDim*iDim*jNode +varDim*iNode + 1];
						tempsum += coeff * ibasis * (-jdbasis) * w2 * invI * 1.e5;
						tempsum += coeff * idbasis * jbasis * w3 * invI  * 1.e2;
					}
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					tempsum += Astate[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis * w4;
					tempsum += Astate[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * w5;
					tempsum += Astate[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * w6;
				}
			}
		}
		
		for (int t=0; t<11; t++) {
			//obstream >> temp;
			*od++ = obsVector[mi+t];
		}
		*od++ = tempsum;
		*od++ = obsVector[mi]-innovation[m];
		qcstream << endl;
		
	}
	QString fileName = "samurai_XY_" + suffix;
	QString outFileName;
	if(QDir::isAbsolutePath(fileName)) {
		outFileName = fileName;
	}
	else {
		outFileName = QDir::current().filePath(fileName);
	}
	
	// Write out the CAPPI to an asi file
	// Initialize header
	int id[511];
	for (int n = 1; n <= 510; n++) {
		id[n]=-999;
	}
	
	// Calculate headers
	QStringList fieldNames;
	fieldNames << "RV" << "HP" << "U" << "V" << "W" << "SF" << "QV" << "RO" << "T" << "P" << "RM" << "RE";
	id[175] = fieldNames.size();
    for(int n = 0; n < id[175]; n++) {
		QString name_1 = fieldNames.at(n).left(1);
        QString name_2 = fieldNames.at(n).mid(1,1);
		int int_1 = *name_1.toAscii().data();
		int int_2 = *name_2.toAscii().data();
		id[176 + (5 * n)] = (int_1 * 256) + int_2;
		id[177 + (5 * n)] = 8224;
		id[178 + (5 * n)] = 8224;
		id[179 + (5 * n)] = 8224;
		id[180 + (5 * n)] = 1;
	}
	
	// Cartesian file
	id[16] = 17217;
	id[17] = 21076;
	
	/* Lat and Lon
	 id[33] = (int)latReference;
	 id[34] = (int)((latReference - (float)id[33]) * 60.);
	 id[35] = (int)((((latReference - (float)id[33]) * 60.) - (float)id[34]) * 60.) * 100;
	 if (lonReference < 0) {
	 lonReference += 360.;
	 }
	 id[36] = (int)lonReference;
	 id[37] = (int)((lonReference - (float)id[36]) * 60.);
	 id[38] = (int)((((lonReference - (float)id[36]) * 60.) - (float)id[37]) * 60.) * 100; */
	
	id[33] = 0;
	id[34] = 0;
	id[35] = 0;
	id[36] = 0;
	id[37] = 0;
	id[38] = 0;
	id[40] = 90;
	
	// Scale factors
	id[68] = 100;
	id[69] = 64;
	
	// X Header
	id[160] = (int)iMin;
	id[161] = (int)iMax*100;
	id[162] = (int)iDim;
	id[163] = (int)(DI * 1000);
	id[164] = 1;
	
	// Y Header
	id[165] = (int)(jMin);
	id[166] = (int)(jMax/10);
	id[167] = (int)jDim;
	id[168] = (int)DJ;
	id[169] = 2;
	
	// Z Header
	id[170] = 1000;
	id[171] = 1000;
	id[172] = 1;
	id[173] = 1000;
	id[174] = 3;
	
	// Number of radars
	id[303] = 1;
	
	// Index of center
	id[309] = (int)(1);
	id[310] = (int)(1);
	id[311] = 0;
	
	// Write ascii file for grid2ps
	//Message::toScreen("Trying to write cappi to "+outFileName);
	outFileName += ".asi";
	QFile asiFile(outFileName);
	if(!asiFile.open(QIODevice::WriteOnly)) {
		cout << "Can't open CAPPI file for writing" << endl;
		return false;
	}
	
	QTextStream out(&asiFile);
	
	// Write header
    int line = 0;
	for (int n = 1; n <= 510; n++) {
		line++;
		out << qSetFieldWidth(8) << id[n];
		if (line == 10) {
			out << endl;
            line = 0;
		}
	}
		
	// Write data
	for(int k = 0; k < 1; k++) {
		out << reset << "level" << qSetFieldWidth(2) << k+1 << endl;
		for(int j = 0; j < jDim; j++) {
			out << reset << "azimuth" << qSetFieldWidth(3) << j+1 << endl;
			for(int n = 0; n < fieldNames.size(); n++) {
				out << reset << left << fieldNames.at(n) << endl;
				int line = 0;
				for (int i = 0; i < iDim;  i++){
					out << reset << qSetRealNumberPrecision(3) << scientific << qSetFieldWidth(10) << fieldNodes[fieldNames.size()*iDim*j + fieldNames.size()*i + n];
					line++;
					if (line == 8) {
						out << endl;
						line = 0;
					}
				}
				if (line != 0) {
					out << endl;
				}
			}
		}
	}
	
	delete[] fieldNodes;
	
	return true;
	
}     

// Basis Functions
float CostFunctionXY_CPU::BasisOri(int m, float x, int M, float xmin, 
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

float CostFunctionXY_CPU::DBasisOri(int m, float x, int M, float xmin, 
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
		b *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
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
				l *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
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
				r *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;	
			}
			b += BC[C][m+3-M] * r;
		}
	}
	return b;
}


float CostFunctionXY_CPU::DDBasisOri(int m, float x, int M, float xmin, 
								  float DX, float DXrecip, int C)
{
	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
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
		b = z;
		z -= 1.0;
		if (z > 0)
			b -= z * 4;
		b *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
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
				z = 2.0 - z;
				l = z;
				z -= 1.0;
				if (z > 0)
					l -= z * 4;
				l *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			}
			b += BC[C][m] * l;
		} else if (m == M-1 || m == M) {
			float r = 0;
			xm = xmin + ((M+1) * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			if (z < 2.0)
			{
				z = 2.0 - z;
				r = z;
				z -= 1.0;
				if (z > 0)
					r -= z * 4;
				r *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			}
			b += BC[C][m+3-M] * r;
		}
	}
	return b;
}

float CostFunctionXY_CPU::DDDBasisOri(int m, float x, int M, float xmin, 
								   float DX, float DXrecip, int C)
{
	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
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
	
	if ((z > 1.0) and (z < 2.0)) {
		b = 1;
	} else if (z < 1.0) {
		b = -3.;
	}
	b *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
	
	// Boundary conditions, if any, are an additional addend.
	if (C >= 0) {
		if (m == 0 || m == 1) {
			float l = 0;
			xm = xmin + (-1 * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			if ((z > 1.0) and (z < 2.0)) {
				l = 1;
			} else if (z < 1.0) {
				l = -3.;
			}
			l *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			
			b += BC[C][m] * l;
		} else if (m == M-1 || m == M) {
			float r = 0;
			xm = xmin + ((M+1) * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			if ((z > 1.0) and (z < 2.0)) {
				r = 1;
			} else if (z < 1.0) {
				r = -3.;
			}
			r *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			b += BC[C][m+3-M] * r;
		}
	}
	return b;
}


real CostFunctionXY_CPU::Basis(int m, real x, int M, real xmin, 
								real DX, real DXrecip, int derivative, 
								int BL, int BR, real lambda)
{
	
	real b = 0;
	real xm = xmin + (m * DX);
	real delta = (x - xm) * DXrecip;
	real z = fabsf(delta);
	real ONESIXTH = 0.16666666666666666666666666667;
	
	switch (derivative) {
		case 0:
			if (z < 2.0)
			{
				z = 2 - z;
				b = (z*z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					b -= (z*z*z) * 4 * ONESIXTH;
			}			
			break;
		case 1:
			if (z < 2.0)
			{
				z = 2.0 - z;
				b = (z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					b -= (z*z) * 4 * ONESIXTH;
				b *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
			}			
			break;
		case 2:
			if (z < 2.0)
			{
				z = 2.0 - z;
				b = z;
				z -= 1.0;
				if (z > 0)
					b -= z * 4;
				b *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			}
			break;
		case 3:
			if ((z > 1.0) and (z < 2.0)) {
				b = 1;
			} else if (z < 1.0) {
				b = -3.;
			}
			b *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			break;
	}
	
	// Add on the boundary conditions
	real bmod = 0;
	int node = -2;
	real coeffmod = 0.;
	if (m == 0) {
		// Left BC
		switch (BL) {
			case 0:
				node = -1;
				coeffmod = -4.;
				break;
			case 1:
				node = -1;
				coeffmod = 0.;
				break;
			case 2:
				node = -1;
				coeffmod = 2.;
				break;
			case 3:
				node = -1;
				coeffmod = -4./(3.*lambda + 1.);
				break;
			case 4:
				// There is no contribution from this node 
				return b;
			case 5:
				// There is no contribution from this node 
				return b;
			case 6:
				// There is no contribution from this node 
				return b;				
		}
	} else if (m == 1) {
		// Left BC
		switch (BL) {
			case 0:
				node = -1;
				coeffmod = -1.;
				break;
			case 1:
				node = -1;
				coeffmod = 1.;
				break;
			case 2:
				node = -1;
				coeffmod = -1.;
				break;
			case 3:
				node = -1;
				coeffmod = (3.*lambda - 1.)/(3.*lambda + 1.);
				break;
			case 4:
				node = -1;
				coeffmod = 1.;
				break;
			case 5:
				node = -1;
				coeffmod = -1.;
				break;				
			case 6:
				// There is no contribution from this node 
				return b;				
		}
				
	} else if (m == M) {
		// Right BC
		switch (BR) {
			case 0:
				node = M+1;
				coeffmod = -4.;
				break;
			case 1:
				node = M+1;
				coeffmod = 0.;
				break;
			case 2:
				node = M+1;
				coeffmod = 2.;
				break;
			case 3:
				node = M+1;
				coeffmod = -4./(3.*lambda + 1.);
				break;
			case 4:
				// There is no contribution from this node 
				return 0.;
			case 5:
				// There is no contribution from this node 
				return 0.;
			case 6:
				// There is no contribution from this node 
				return 0.;				
		}
	} else if (m == M-1) {
		// Left BC
		switch (BL) {
			case 0:
				node = M+1;
				coeffmod = -1.;
				break;
			case 1:
				node = M+1;
				coeffmod = 1.;
				break;
			case 2:
				node = M+1;
				coeffmod = -1.;
				break;
			case 3:
				node = M+1;
				coeffmod = (3.*lambda - 1.)/(3.*lambda + 1.);
				break;
			case 4:
				node = M+1;
				coeffmod = 1.;
				break;
			case 5:
				node = M+1;
				coeffmod = -1.;
				break;				
			case 6:
				// There is no contribution from this node 
				return 0.;				
		}
	}
	
	xm = xmin + (node * DX);
	delta = (x - xm) * DXrecip;
	z = fabsf(delta);
	switch (derivative) {
		case 0:
			if (z < 2.0)
			{
				z = 2 - z;
				bmod = (z*z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					bmod -= (z*z*z) * 4 * ONESIXTH;
			}			
			break;
		case 1:
			if (z < 2.0)
			{
				z = 2.0 - z;
				bmod = (z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					bmod -= (z*z) * 4 * ONESIXTH;
				bmod *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
			}			
			break;
		case 2:
			if (z < 2.0)
			{
				z = 2.0 - z;
				bmod = z;
				z -= 1.0;
				if (z > 0)
					bmod -= z * 4;
				bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			}
			break;
		case 3:
			if ((z > 1.0) and (z < 2.0)) {
				bmod = 1;
			} else if (z < 1.0) {
				bmod = -3.;
			}
			bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			break;
	}
	
	b += coeffmod * bmod;
	
	// R2 needs one more addition
	if ((m == 1) and (BL == 4)) {
		node = 0;
		coeffmod = -0.5;
	} else if ((m == M-1) and (BR == 4)) {
		node = M;
		coeffmod = -0.5;
	} else {
		return b;
	}
	
	xm = xmin + (node * DX);
	delta = (x - xm) * DXrecip;
	z = fabsf(delta);
	switch (derivative) {
		case 0:
			if (z < 2.0)
			{
				z = 2 - z;
				bmod = (z*z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					bmod -= (z*z*z) * 4 * ONESIXTH;
			}			
			break;
		case 1:
			if (z < 2.0)
			{
				z = 2.0 - z;
				bmod = (z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					bmod -= (z*z) * 4 * ONESIXTH;
				bmod *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
			}			
			break;
		case 2:
			if (z < 2.0)
			{
				z = 2.0 - z;
				bmod = z;
				z -= 1.0;
				if (z > 0)
					bmod -= z * 4;
				bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			}
			break;
		case 3:
			if ((z > 1.0) and (z < 2.0)) {
				bmod = 1;
			} else if (z < 1.0) {
				bmod = -3.;
			}
			bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			break;
	}
	
	b += coeffmod * bmod;
	
	return b;
	
}
