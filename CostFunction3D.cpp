/*
 *  CostFunction3D.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunction3D.h"
#include "MetObs.h"
#include <cmath>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QFile>
#include <QDir>

CostFunction3D::CostFunction3D(const int& numObs, const int& stateSize)
	: CostFunction(numObs, stateSize)
{
}

CostFunction3D::~CostFunction3D()
{
}

void CostFunction3D::finalize()
{	

	delete iFilter;
	delete jFilter;
	delete kFilter;
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

void CostFunction3D::initialize(const real& imin, const real& imax, const int& idim,
									const real& jmin, const real& jmax, const int& jdim,
									const real& kmin, const real& kmax, const int& kdim,
									real* bgU, real* obs)
{

	// Initialize number of control variables -- one less than physical variables
	varDim = 6;
	
	// Set the constant height for the XY case
	constHeight = 0.;
	
	// Initialize background errors and filter scales
	bgError[0] = 200.;
	bgError[1] = 200.;
	bgError[2] = 10.0;
	bgError[3] = 5.0;
	bgError[4] = 3.0;	
	bgError[5] = 3.0;
	
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
	kMin = kmin;
	kMax = kmax;
	kDim = kdim;
	DI = (iMax - iMin) / (iDim - 1);
	DIrecip = 1./DI;
	DJ = (jMax - jMin) / (jDim - 1);
    DJrecip = 1./DJ;
	DK = (kMax - kMin) / (kDim - 1);
    DKrecip = 1./DK;
	
	// Set up the recursive filter
	real iFilterScale = 8.;
	real jFilterScale = 8.;
	real kFilterScale = 4.;	
	bgErrorScale = 1.;
	iFilter = new RecursiveFilter(4,iFilterScale);
	jFilter = new RecursiveFilter(4,jFilterScale);
	kFilter = new RecursiveFilter(4,kFilterScale);
	
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
	bState = iDim*jDim*kDim*varDim;
	CTHTd = new real[nState];
	stateU = new real[nState];
	Uprime = new real[nState];
	HCq = new real[mObs];
	obsVector = new real[mObs*12];
	
	bgState = new real[bState];
	stateB = new real[bState];
	stateA = new real[bState];
	stateC = new real[bState];
	iTemp = new real[iDim];
	jTemp = new real[jDim];
	kTemp = new real[kDim];
	
	// Set up the spline matrices
	setupSplines();

}	

void CostFunction3D::initState()
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
	SBtransform(bgFields, stateB);
		
	// SA transform = bg B's -> bg A's
	SAtransform(stateB, bgState);
	
	// Compute and display the variable RMS and BG errors for reference
	for (int var = 0; var < varDim; var++) {
		double varScale = 0;
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					int bIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var;
					varScale += bgState[bIndex] * bgState[bIndex];
				}
			}
		}
		varScale = sqrt(varScale/(iDim*jDim*kDim));
		if (varScale) {
			double errPct = 100*bgError[var]/varScale;
			cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var] 
			<< " ( " << errPct << " %)" << endl;
		} else {
			cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var] 
			<< " ( Infinite! %)" << endl;
		}
	}	
	
	// Load the obs locally and weight the nonlinear observation operators by interpolated bg fields
	obAdjustments();
	
	// d = y - HXb
	calcInnovation();

	// Output the original background field
	outputAnalysis("background", bgState , true);

	// HTd
	calcHTranspose(innovation, stateC);
	
	SCtranspose(stateC, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtransform(stateA, stateB);
	
	SBtranspose(stateB, CTHTd);
			
}	

double CostFunction3D::funcValue(double* state)
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
		obIP += (HCq[m]-innovation[m])*(obsVector[m*12+1])*(HCq[m]-innovation[m]);
	}
	double J = 0.5*qIP + 0.5*obIP;
	return J;
	
}

void CostFunction3D::funcGradient(double* state, double* gradient)
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

void CostFunction3D::updateHCq(double* state)
{

	// SB transform from the q's
	SBtransform(state, stateB);
	
	// S (SA transform) yield A's, put it in the tempState
	SAtransform(stateB, stateA);
	
	SCtransform(stateA, stateC);
	
	// H
	#pragma omp parallel for
	for (int m = 0; m < mObs; m++) {
		int mi = m*12;
		real i = obsVector[mi+2];
		real j = obsVector[mi+3];
		real k = obsVector[mi+4];
		real tempsum = 0;
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		int kk = (int)((k - kMin)*DKrecip);
		real ibasis = 0;
		real jbasis = 0;
		real kbasis = 0;
		real idbasis = 0;
		real jdbasis = 0;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
					if ((iNode < 0) or (iNode >= iDim) or 
						(jNode < 0) or (jNode >= jDim) or
						(kNode < 0) or (kNode >= kDim)) continue;
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
					int cIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
					real basis3x = ibasis * jbasis * kbasis;
					for (int var = 0; var < varDim; var++) {
						tempsum += stateC[cIndex + var] *  basis3x * obsVector[mi+6 + var];
					}
				}
			}
		}
		HCq[m] = tempsum;
	}	
	
}

void CostFunction3D::updateBG()
{

	// SB transform from the q's
	SBtransform(currState, stateB);
	
	// S (SA transform) yield A's
	SAtransform(stateB, stateA);
	
	SCtransform(stateA, stateC);
	
	outputAnalysis("increment", stateC, false);	
	
	// In BG update we are directly summing C + A
	ofstream cstream("CoeffAnalysis.out");
	cstream << "Variable\tI\tJ\tK\tBackground\tAnalysis\tIncrement\n";
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					cstream << var << "\t" << iIndex << "\t" << jIndex << "\t" << kIndex << "\t";
					int bgIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var;
					cstream << bgState[bgIndex] << "\t";
					bgState[bgIndex] += stateC[bgIndex];
					cstream << bgState[bgIndex] << "\t"; 
					cstream << stateC[bgIndex] << endl;
				}
			}
		}
	}
		
	outputAnalysis("analysis", bgState, true);
	
}

void CostFunction3D::calcInnovation()
{
	// Initialize and fill the innovation vector
	cout << "Initializing innovation vector..." << endl;
	for (int m = 0; m < mObs; m++) {
		HCq[m] = 0.0;
		innovation[m] = obsVector[m*12];
	}
	
	real innovationRMS = 0;
	#pragma omp parallel for reduction(+:innovationRMS)
	for (int m = 0; m < mObs; m++) {
		int mi = m*12;
		real i = obsVector[mi+2];
		real j = obsVector[mi+3];
		real k = obsVector[mi+4];
		real tempsum = 0;
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		int kk = (int)((k - kMin)*DKrecip);
		real ibasis = 0;
		real jbasis = 0;
		real kbasis = 0;
		real idbasis = 0;
		real jdbasis = 0;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
					if ((iNode < 0) or (iNode >= iDim) or 
						(jNode < 0) or (jNode >= jDim) or
						(kNode < 0) or (kNode >= kDim)) continue;
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
					int stateIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
					real basis3x = ibasis * jbasis * kbasis;
					for (int var = 0; var < varDim; var++) {
						tempsum += bgState[stateIndex + var] * basis3x * obsVector[mi+6+var];
					}
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

void CostFunction3D::calcHTranspose(const real* yhat, real* Astate) 
{
	
	// Clear the Astate
	for (int b = 0; b < bState; b++) {
		Astate[b] = 0.;
	}
	
	// Calculate H Transpose	
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int kIndex = 0; kIndex < kDim; kIndex++) {
				#pragma omp parallel for
				for (int m = 0; m < mObs; m++) {
					// Sum over obs this time
					// Multiply state by H weights
					int mi = m*12;
					real qhat = yhat[m];
					real invError = obsVector[mi+1];
					real i = obsVector[mi+2];
					real j = obsVector[mi+3];
					real k = obsVector[mi+4];
					int ii = (int)((i - iMin)*DIrecip);
					if ((iIndex < ii-1) or (iIndex > ii+2)) continue;
					int jj = (int)((j - jMin)*DJrecip);
					if ((jIndex < jj-1) or (jIndex > jj+2)) continue;
					int kk = (int)((k - kMin)*DKrecip);
					if ((kIndex < kk-1) or (kIndex > kk+2)) continue;
					real ibasis = 0;
					real jbasis = 0;
					real kbasis = 0;
					ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					kbasis = Basis(kIndex, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
					int stateIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex;
					real qbasise = qhat* ibasis * jbasis * kbasis *invError;
					for (int var = 0; var < varDim; var++) {
						#pragma omp atomic
						Astate[stateIndex + var] += qbasise * obsVector[mi+6+var];
					}
				}
			}
		}
	}
}

bool CostFunction3D::SAtransform(real* Bstate, real* Astate)
{

	#pragma omp parallel for
	for (int var = 0; var < varDim; var++) {
		int l;
		real* kB = new real[kDim];
		real* x = new real[kDim];
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					kB[kIndex] = Bstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var];
				}
				// Solve for A's using compact storage
				real sum = 0;
				for (int k = 0; k < kDim; k++) {
					for (sum=kB[k], l=-1;l>=-3;l--) {
						if ((k+l >= 0) and ((k*4-l) >= 0))
							sum -= kL[kDim*4*var + k*4-l]*x[k+l];
					}
					x[k] = sum/kL[kDim*4*var + k*4];
				}	
				for (int k=kDim-1;k>=0;k--) {
					for (sum=x[k], l=1;l<=3;l++) {
						if ((k+l < kDim) and (((k+l)*4+l) < kDim*4))
							sum -= kL[kDim*4*var + (k+l)*4+l]*x[k+l];
					}
					x[k] = sum/kL[kDim*4*var + k*4];
				}
				
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = x[kIndex]; 
				}
			}
		}
		delete[] kB;
		delete[] x;
		
		real* jB = new real[jDim];
		x = new real[jDim];
		for (int kIndex = 0; kIndex < kDim; kIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					jB[jIndex] = Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var];
				}
				// Solve for A's using compact storage
				real sum = 0;
				for (int j = 0; j < jDim; j++) {
					for (sum=jB[j], l=-1;l>=-3;l--) {
						if ((j+l >= 0) and ((j*4-l) >= 0))
							sum -= jL[jDim*4*var + j*4-l]*x[j+l];
					}
					x[j] = sum/jL[jDim*4*var + j*4];
				}	
				for (int j=jDim-1;j>=0;j--) {
					for (sum=x[j], l=1;l<=3;l++) {
						if ((j+l < jDim) and (((j+l)*4+l) < jDim*4))
							sum -= jL[jDim*4*var + (j+l)*4+l]*x[j+l];
					}
					x[j] = sum/jL[jDim*4*var + j*4];
				}
				
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = x[jIndex]; 
				}
			}
		}
		delete[] jB;
		delete[] x;
		
		real* iB = new real[iDim];
		x = new real[iDim];
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int kIndex = 0; kIndex < kDim; kIndex++) {
				for (int iIndex = 0; iIndex < iDim; iIndex++) {
					iB[iIndex] = Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var];
				}
				// Solve for A's using compact storage
				real sum = 0;
				for (int i = 0; i < iDim; i++) {
					for (sum=iB[i], l=-1;l>=-3;l--) {
						if ((i+l >= 0) and ((i*4-l) >= 0))
							sum -= iL[iDim*4*var + i*4-l]*x[i+l];
					}
					x[i] = sum/iL[iDim*4*var + i*4];
				}	
				for (int i=iDim-1;i>=0;i--) {
					for (sum=x[i], l=1;l<=3;l++) {
						if ((i+l < iDim) and (((i+l)*4+l) < iDim*4))
							sum -= iL[iDim*4*var + (i+l)*4+l]*x[i+l];
					}
					x[i] = sum/iL[iDim*4*var + i*4];
				}
				
				for (int iIndex = 0; iIndex < iDim; iIndex++) {
					Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = x[iIndex]; 
				}
				
			}
		}
		delete[] iB;
		delete[] x;
		
	}

	return true;
}		

void CostFunction3D::SBtransform(real* Ustate, real* Bstate)
{
	// Clear the Bstate
	for (int b = 0; b < bState; b++) {
		Bstate[b] = 0.;
	}
					
	#pragma omp parallel for
	for (int var = 0; var < varDim; var++) {
		
		for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
			for (int imu = -1; imu <= 1; imu += 2) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
				int ii = (int)((i - iMin)*DIrecip);
				for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
					if ((iNode < 0) or (iNode >= iDim)) continue;
					
					for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
						for (int jmu = -1; jmu <= 1; jmu += 2) {
							real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5));
							int jj = (int)((j - jMin)*DJrecip);
							for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
								if ((jNode < 0) or (jNode >= jDim)) continue;
								
								for (int kIndex = 0; kIndex < (kDim-1); kIndex++) {
									for (int kmu = -1; kmu <= 1; kmu += 2) {
										real k = kMin + DK * (kIndex + (0.5*sqrt(1./3.) * kmu + 0.5));
										int kk = (int)((k - kMin)*DKrecip);
										for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
											if ((kNode < 0) or (kNode >= kDim)) continue;
											
											real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
											real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
											real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
											int bgI = iIndex*2 + (imu+1)/2;
											int bgJ = jIndex*2 + (jmu+1)/2;
											int bgK = kIndex*2 + (kmu+1)/2;
											#pragma omp atomic
											Bstate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + var] += 
												0.125 * Ustate[varDim*(iDim-1)*2*(jDim-1)*2*bgK +varDim*(iDim-1)*2*bgJ +varDim*bgI + var] 
													* ibasis * jbasis * kbasis; 
										}
									}	
								}
							}
						}
					}	
				}
			}
		}
	}
}


void CostFunction3D::SBtranspose(real* Bstate, real* Ustate)
{
		
	// Clear the Ustate
	for (int b = 0; b < bState; b++) {
		Ustate[b] = 0;
	}
	
	#pragma omp parallel for
	for (int var = 0; var < varDim; var++) {
		
		for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
			for (int imu = -1; imu <= 1; imu += 2) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
				int ii = (int)((i - iMin)*DIrecip);
				for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
					if ((iNode < 0) or (iNode >= iDim)) continue;
					
					for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
						for (int jmu = -1; jmu <= 1; jmu += 2) {
							real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5));
							int jj = (int)((j - jMin)*DJrecip);
							for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
								if ((jNode < 0) or (jNode >= jDim)) continue;
								
								for (int kIndex = 0; kIndex < (kDim-1); kIndex++) {
									for (int kmu = -1; kmu <= 1; kmu += 2) {
										real k = kMin + DK * (kIndex + (0.5*sqrt(1./3.) * kmu + 0.5));
										int kk = (int)((k - kMin)*DKrecip);
										for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
											if ((kNode < 0) or (kNode >= kDim)) continue;
																					
											real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
											real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
											real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
											int bgI = iIndex*2 + (imu+1)/2;
											int bgJ = jIndex*2 + (jmu+1)/2;
											int bgK = kIndex*2 + (kmu+1)/2;
											#pragma omp atomic										
											Ustate[varDim*(iDim-1)*2*(jDim-1)*2*bgK +varDim*(iDim-1)*2*bgJ +varDim*bgI + var] += 
											0.125 * Bstate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + var] 
											* ibasis * jbasis * kbasis; 
										}
									}	
								}
							}
						}
					}	
				}
			}
		}
	}
}

void CostFunction3D::SCtransform(real* Astate, real* Cstate)
{
	// Isoptropic Recursive filter for speed, no anisotropic "triad" working yet 
	for (int var = 0; var < varDim; var++) {
		//FK
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					kTemp[kIndex] = Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
				}
				kFilter->filterArray(kTemp, kDim);
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = kTemp[kIndex];
				}
			}
		}
		//FJ
		for (int kIndex = 0; kIndex < kDim; kIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					jTemp[jIndex] = Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
				}
				jFilter->filterArray(jTemp, jDim);
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
				}
			}
		}
		//FI
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int kIndex = 0; kIndex < kDim; kIndex++) {
				for (int iIndex = 0; iIndex < iDim; iIndex++) {
					iTemp[iIndex] = Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
				}
				iFilter->filterArray(iTemp, iDim);
				for (int iIndex = 0; iIndex < iDim; iIndex++) {
					// D
					Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex] * bgError[var]; 
				}
			}
		}		
	}
}

void CostFunction3D::SCtranspose(real* Cstate, real* Astate)
{
	
	// Isoptropic Recursive filter for speed, no anisotropic "triad" working yet 
	for (int var = 0; var < varDim; var++) {
		//FI & D
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int kIndex = 0; kIndex < kDim; kIndex++) {
				for (int iIndex = 0; iIndex < iDim; iIndex++) {
					iTemp[iIndex] = Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] * bgError[var];
				}
				iFilter->filterArray(iTemp, iDim);
				for (int iIndex = 0; iIndex < iDim; iIndex++) {
					Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex]; 
				}
			}
		}
		//FJ
		for (int kIndex = 0; kIndex < kDim; kIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					jTemp[jIndex] = Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
				}
				jFilter->filterArray(jTemp, jDim);
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
				}
			}
		}
		//FK
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					kTemp[kIndex] = Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
				}
				kFilter->filterArray(kTemp, kDim);
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = kTemp[kIndex];
				}
			}
		}
	}		
}


bool CostFunction3D::setupSplines()
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
					real pm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					real qm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 3, R1T2, R1T2);
					real pn, qn;
					P[iNode][iNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((iNode+1) < iDim) {
						pn = Basis(iNode+1, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
						qn = Basis(iNode+1, i, iDim-1, iMin, DI, DIrecip, 3, R1T2, R1T2);
						P[iNode][iNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+1][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((iNode+2) < iDim) {
						pn = Basis(iNode+2, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
						qn = Basis(iNode+2, i, iDim-1, iMin, DI, DIrecip, 3, R1T2, R1T2);
						P[iNode][iNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+2][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((iNode+3) < iDim) {
						pn = Basis(iNode+3, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
						qn = Basis(iNode+3, i, iDim-1, iMin, DI, DIrecip, 3, R1T2, R1T2);
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
					real pm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					real qm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 3, R1T2, R1T2);
					real pn, qn;
					P[jNode][jNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((jNode+1) < jDim) {
						pn = Basis(jNode+1, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
						qn = Basis(jNode+1, j, jDim-1, jMin, DJ, DJrecip, 3, R1T2, R1T2);
						P[jNode][jNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+1][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((jNode+2) < jDim) {
						pn = Basis(jNode+2, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
						qn = Basis(jNode+2, j, jDim-1, jMin, DJ, DJrecip, 3, R1T2, R1T2);
						P[jNode][jNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+2][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));						
					}
					if ((jNode+3) < jDim) {
						pn = Basis(jNode+3, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
						qn = Basis(jNode+3, j, jDim-1, jMin, DJ, DJrecip, 3, R1T2, R1T2);
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
	
	// Clear the temporary memory and reallocate for kDim
	for (int j = 0; j < jDim; j++) {
		delete[] P[j];
	}
	delete[] P;
	delete[] p;
	
	P = new real*[kDim];
	p = new real[kDim];
	kL = new real[varDim*kDim*4];
	for (int k = 0; k < kDim; k++) {
		P[k] = new real[kDim];
		p[k] = 0.;
	}
	
	for (int k = 0; k < varDim*kDim*4; k++) {
		kL[k] = 0;
	}
	
	cutoff_wl = 4;
	eq = pow( (cutoff_wl/(2*Pi)) , 6);
	for (int var = 0; var < varDim; var++) {
		
		for (int i = 0; i < kDim; i++) {
			for (int j = 0; j < kDim; j++) {
				P[i][j] = 0;
			}
		}		
		for (int kIndex = 0; kIndex < (kDim-1); kIndex++) {
			for (int kmu = -1; kmu <= 1; kmu += 2) {
				real k = kMin + DK * (kIndex + (0.5*sqrt(1./3.) * kmu + 0.5));
				int kk = (int)((k - kMin)*DKrecip);
				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
					if ((kNode < 0) or (kNode >= kDim)) continue;
					real pm = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
					real qm = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 3, R1T2, R1T2);
					real pn, qn;
					P[kNode][kNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((kNode+1) < kDim) {
						pn = Basis(kNode+1, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
						qn = Basis(kNode+1, k, kDim-1, kMin, DK, DKrecip, 3, R1T2, R1T2);
						P[kNode][kNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[kNode+1][kNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((kNode+2) < kDim) {
						pn = Basis(kNode+2, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
						qn = Basis(kNode+2, k, kDim-1, kMin, DK, DKrecip, 3, R1T2, R1T2);
						P[kNode][kNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[kNode+2][kNode] += 0.5 * ((pm * pn) + eq * (qm * qn));						
					}
					if ((kNode+3) < kDim) {
						pn = Basis(kNode+3, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
						qn = Basis(kNode+3, k, kDim-1, kMin, DK, DKrecip, 3, R1T2, R1T2);
						P[kNode][kNode+3] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[kNode+3][kNode] += 0.5 * ((pm * pn) + eq * (qm * qn));						
					}
				}
			}
		}
		
		// Cholesky decomp
		for (int i=0;i<kDim;i++) {
			for (int j=i;j<kDim;j++) {
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
		
		for (int k = 0; k < kDim; k++) {
			kL[kDim*4*var + k*4] = p[k];
			//cout << jL[jDim*4*var + j*4] << " ";
			for (int n=1;n<4;n++) {
				if ((k-n) >= 0) {
					kL[kDim*4*var + k*4+n] = P[k][k-n];
				}
				//cout << jL[jDim*4*var + j*4 + n] << " ";
			} //cout << endl;
		} //cout << endl;
	}
	
	// Clear the temporary memory
	for (int k = 0; k < kDim; k++) {
		delete[] P[k];
	}
	delete[] P;
	delete[] p;
	
	
	/*
	//Analytic Test set		
	p = new real[iDim];
	for (int i = 0; i < iDim; i++) {
		p[i] = 0.;
	}
	int var = 0;
	for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
		for (int imu = -1; imu <= 1; imu += 2) {
			real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
			int ii = (int)((i - iMin)*DIrecip);
			for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
				if ((iNode < 0) or (iNode >= iDim)) continue;
				real pm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
				if (fabs(i) < 45) {
					p[iNode] += 0.5*pm*0.005;
				} else {
					p[iNode] = 0;
				}
			}
		}
	}
	real* y = new real[iDim];
	// Solve for A's
	real sum = 0;
	int k;
		 
	//Solve for A's using compact storage
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
			si = (si-y[0])/1e6;
			//real analytic = sin(i*2*Pi/iDim) + sin(2*i*2*Pi/iDim) + sin(4*i*2*Pi/iDim) + sin(8*i*2*Pi/iDim) + sin(16*i*2*Pi/iDim); 
			real analytic = i*i;
			cout << i << "\t" << p[iIndex] << "\t" << y[iIndex] << "\t" << si << "\t" << analytic << "\t" << si - analytic << endl;
		}
	} */
	
	return true;
	
}

bool CostFunction3D::outputAnalysis(const QString& suffix, real* Astate, bool updateMish)
{
	
	real* fieldNodes = new real[12*iDim*jDim*kDim];	
		
	// H --> to Mish for output
	QString samuraiout = "SAMURAI_" + suffix + ".out";
	ofstream samuraistream(samuraiout.toAscii().data());
	samuraistream << "X\tY\tZ\trhoE\tu\tv\tw\tVorticity\tDivergence\tqv\trho\tT\tP\th\n";
	samuraistream.precision(10);
	real CoriolisF = 6e-5;
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		for (int ihalf = 0; ihalf <=1; ihalf++) {
			for (int imu = -ihalf; imu <= ihalf; imu++) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5*ihalf));
				if (i > ((iDim-1)*DI + iMin)) continue;
				
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					for (int jhalf =0; jhalf <=1; jhalf++) {
						for (int jmu = -jhalf; jmu <= jhalf; jmu++) {
							real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5*jhalf));
							if (j > ((jDim-1)*DJ + jMin)) continue;	
							
							for (int kIndex = 0; kIndex < kDim; kIndex++) {
								for (int khalf =0; khalf <=1; khalf++) {
									for (int kmu = -khalf; kmu <= khalf; kmu++) {
										real k = kMin + DK * (kIndex + (0.5*sqrt(1./3.) * kmu + 0.5*khalf));
										if (k > ((kDim-1)*DK + kMin)) continue;	
										
										real rhoBar = rhoBase*exp(-rhoInvScaleHeight*k);
										real qBar = 19.562 - 0.004066*k + 7.8168e-7*k*k;
										real hBar = 3.5e5;
										
										int ii = (int)((i - iMin)*DIrecip);
										int jj = (int)((j - jMin)*DJrecip);
										int kk = (int)((k - kMin)*DKrecip);
										real ibasis = 0.;
										real jbasis = 0.;
										real kbasis = 0.;
										real idbasis = 0.;
										real jdbasis = 0.;
										real rhov = 0.;
										real rhou = 0.;
										real rhow = 0.;
										real chi = 0.;
										real hprime = 0.;
										real qvprime = 0.;
										real rhoprime = 0.;
										real vorticity = 0.;
										real divergence = 0.;
										for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
											for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
												for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
													if ((iNode < 0) or (iNode >= iDim) or 
														(jNode < 0) or (jNode >= jDim) or
														(kNode < 0) or (kNode >= kDim)) continue;
													ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
													jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
													kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
													idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, R1T2, R1T2);
													jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, R1T2, R1T2);
													real basis3x = ibasis*jbasis*kbasis;
													real ucoeff = Astate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode];
													real vcoeff = Astate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 1];
													rhou +=  Astate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 0] * basis3x;
													rhov +=  Astate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 1] * basis3x;

													vorticity += 100 * (vcoeff * idbasis * jbasis - ucoeff * ibasis * jdbasis); //10-5 s-1;
													divergence += 100 * (ucoeff * idbasis * jbasis + vcoeff * ibasis * jdbasis);
													rhow += Astate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 2] * basis3x;
													hprime += Astate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 3] * basis3x;
													qvprime += Astate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 4] * basis3x;
													rhoprime += Astate[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 5] * basis3x;
												}
											}
										}
										// Get it to relevant variables for the flux calculation
										real rhoa = rhoBar + rhoprime / 100;
										real rhoq = (qBar + qvprime) * rhoa / 1000.;
										real rho = rhoa + rhoq;
										real v = rhov / rho;
										real u = rhou / rho;
										real w = rhow / rho;
										real wspd = sqrt(u*u + v*v);
										real KE = 0.5*rho*(v*v + u*u + w*w);
										real rhoE = rho*(hBar + hprime*1.e3) + KE;
										real T = ((hBar + hprime*1.e3) - 2.501e3*(qBar + qvprime) - 9.81*j)/1005.7;
										real press = T*rhoa*287./100.;
										//real e = rhoE/rho;
										
										// Output it
										if (updateMish and imu and jmu and kmu and ihalf and jhalf and khalf) {
											int bgI = iIndex*2 + (imu+1)/2;
											int bgJ = jIndex*2 + (jmu+1)/2;
											int bgK = kIndex*2 + (kmu+1)/2;
											bgFields[varDim*(iDim-1)*2*(jDim-1)*2*bgK + varDim*(iDim-1)*2*bgJ +varDim*bgI + 0] = rhou;
											bgFields[varDim*(iDim-1)*2*(jDim-1)*2*bgK + varDim*(iDim-1)*2*bgJ +varDim*bgI + 1] = rhov;
											bgFields[varDim*(iDim-1)*2*(jDim-1)*2*bgK + varDim*(iDim-1)*2*bgJ +varDim*bgI + 2] = rhow;
											bgFields[varDim*(iDim-1)*2*(jDim-1)*2*bgK + varDim*(iDim-1)*2*bgJ +varDim*bgI + 3] = hprime;
											bgFields[varDim*(iDim-1)*2*(jDim-1)*2*bgK + varDim*(iDim-1)*2*bgJ +varDim*bgI + 4] = qvprime;
											bgFields[varDim*(iDim-1)*2*(jDim-1)*2*bgK + varDim*(iDim-1)*2*bgJ +varDim*bgI + 5] = rhoprime;
										}
										
										if (!ihalf and !jhalf and !khalf){
											// On the node
											samuraistream << scientific << i << "\t" << j << "\t"  << k << "\t" << rhoE
											<< "\t" << u << "\t" << v << "\t" << w << "\t" << vorticity << "\t" << divergence
											<< "\t" << qvprime << "\t" << rhoprime << "\t" << T << "\t" << press <<  "\t" << hprime
											<< "\n";
											
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 0] = wspd;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 1] = chi;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 2] = u;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 3] = v;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 4] = w;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 5] = vorticity;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 6] = qvprime;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 7] = rhoprime;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 8] = T;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 9] = press;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 10] = hprime;
											fieldNodes[12*iDim*jDim*kIndex + 12*iDim*jIndex + 12*iIndex + 11] = divergence;
										}
									}
								}
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
	*os++ = "X";
	*os++ = "Y";
	*os++ = "Z";
	*os++ = "Type";
	*os++ = "Weight 1";
	*os++ = "Weight 2";
	*os++ = "Weight 3";
	*os++ = "Weight 4";
	*os++ = "Weight 5";
	*os++ = "Weight 6";	
	*os++ = "Analysis";
	*os++ = "Background";
	qcstream << endl;
	
	ostream_iterator<double> od(qcstream, "\t ");
	for (int m = 0; m < mObs; m++) {
		int mi = m*12;
		real i = obsVector[mi+2];
		real j = obsVector[mi+3];
		real k = obsVector[mi+4];
		real tempsum = 0;
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		int kk = (int)((k - kMin)*DKrecip);
		real ibasis = 0;
		real jbasis = 0;
		real kbasis = 0;
		real idbasis = 0;
		real jdbasis = 0;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
					if ((iNode < 0) or (iNode >= iDim) or 
						(jNode < 0) or (jNode >= jDim) or
						(kNode < 0) or (kNode >= kDim)) continue;
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
					int cIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
					tempsum += stateC[cIndex + 0] * ibasis * jbasis * kbasis * obsVector[mi+6];
					tempsum += stateC[cIndex + 1] * ibasis * jbasis * kbasis * obsVector[mi+7];
					tempsum += stateC[cIndex + 2] * ibasis * jbasis * kbasis * obsVector[mi+8];
					tempsum += stateC[cIndex + 3] * ibasis * jbasis * kbasis * obsVector[mi+9];
					tempsum += stateC[cIndex + 4] * ibasis * jbasis * kbasis * obsVector[mi+10];
					tempsum += stateC[cIndex + 5] * ibasis * jbasis * kbasis * obsVector[mi+11];
				}
			}
		}
				
		for (int t=0; t<12; t++) {
			//obstream >> temp;
			*od++ = obsVector[mi+t];
		}
		*od++ = tempsum;
		*od++ = obsVector[mi]-innovation[m];
		qcstream << endl;
		
	}
	QString fileName = "samurai_XYZ_" + suffix;
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
	fieldNames << "WS" << "CH" << "U" << "V" << "W" << "VO" << "QV" << "RO" << "T" << "P" << "HP" << "DV";
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
	id[160] = (int)iMin*100;
	id[161] = (int)iMax*100;
	id[162] = (int)iDim;
	id[163] = (int)(DI * 1000);
	id[164] = 1;
	
	// Y Header
	id[165] = (int)jMin*100;
	id[166] = (int)jMax*100;
	id[167] = (int)jDim;
	id[168] = (int)(DJ * 1000);
	id[169] = 2;
	
	// Z Header
	id[170] = (int)kMin*1000;
	id[171] = (int)kMax*1000;
	id[172] = (int)kDim;
	id[173] = int(DK * 1000);
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
	for(int k = 0; k < kDim; k++) {
		out << reset << "level" << qSetFieldWidth(2) << k+1 << endl;
		for(int j = 0; j < jDim; j++) {
			out << reset << "azimuth" << qSetFieldWidth(3) << j+1 << endl;
			for(int n = 0; n < fieldNames.size(); n++) {
				out << reset << left << fieldNames.at(n) << endl;
				int line = 0;
				for (int i = 0; i < iDim;  i++){
					out << reset << qSetRealNumberPrecision(3) << scientific << qSetFieldWidth(10) << 
						fieldNodes[fieldNames.size()*iDim*jDim*k + fieldNames.size()*iDim*j + fieldNames.size()*i + n];
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

void CostFunction3D::obAdjustments() {
	
	// Load the obs locally and weight the nonlinear observation operators by interpolated bg fields
	for (int m = 0; m < mObs; m++) {
		int mi = m*12;
		for (int ob = 0; ob < 12; ob++) {
			obsVector[mi+ob] = rawObs[mi+ob];
		}
		real type = obsVector[mi+5];
		if (type <= 1) continue; 
		
		real i = obsVector[mi+2];
		real j = obsVector[mi+3];
		real k = obsVector[mi+4];
		real rhoprime = 0.;
		real qvprime = 0.;
		
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		int kk = (int)((k - kMin)*DKrecip);
		real ibasis = 0.;
		real jbasis = 0.;
		real kbasis = 0.;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
				if ((iNode < 0) or (iNode >= iDim) or 
					(jNode < 0) or (jNode >= jDim) or
					(kNode < 0) or (kNode >= kDim)) continue;
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, R1T2, R1T2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, R1T2, R1T2);
					kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, R1T2, R1T2);
					qvprime += bgState[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * kbasis;
					rhoprime += bgState[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * kbasis;
				}
			}
		}
		real rhoBar = rhoBase*exp(-rhoInvScaleHeight*k);
		real qBar = 19.562 - 0.004066*k + 7.8168e-7*k*k;
		real rhoa = rhoBar + rhoprime / 100;
		real rhoq = (qBar + qvprime) * rhoa / 1000.;
		real rhoBG = rhoa + rhoq;
		
		// Adjust the relevant fields
		if (type == MetObs::sfmr) {
			obsVector[mi] *= rhoBG;
			// Still need to think carefully about nonlinear SFMR operator in Cartesian space
			//real wsBG = vBG*vBG; // + uBG*uBG;
			//obsVector[mi+2] = 2.*vBG/wsBG;
			//obsVector[mi+3] = 0.; //2.*uBG/wsBG;
		}
		if ((type == MetObs::radar) or (type == MetObs::qscat) or (type == MetObs::ascat)) { 
			obsVector[mi] *= rhoBG;
		}
	}
}	

real CostFunction3D::Basis(int m, real x, int M, real xmin, 
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
				b *= DXrecip * DXrecip;
			}
			break;
		case 3:
			if ((z > 1.0) and (z < 2.0)) {
				b = 1;
			} else if (z < 1.0) {
				b = -3.;
			}
			b *= ((delta > 0) ? -1.0 : 1.0) * DXrecip * DXrecip * DXrecip;
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
				bmod *= DXrecip * DXrecip;
			}
			break;
		case 3:
			if ((z > 1.0) and (z < 2.0)) {
				bmod = 1;
			} else if (z < 1.0) {
				bmod = -3.;
			}
			bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip * DXrecip * DXrecip;
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
				bmod *= DXrecip * DXrecip;
			}
			break;
		case 3:
			if ((z > 1.0) and (z < 2.0)) {
				bmod = 1;
			} else if (z < 1.0) {
				bmod = -3.;
			}
			bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip * DXrecip * DXrecip;
			break;
	}
	
	b += coeffmod * bmod;
	
	return b;
	
}
