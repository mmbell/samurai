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
#include <netcdfcpp.h>

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
	delete[] iL;
	delete[] jL;
	delete[] stateU;
	delete[] stateV;
	delete[] fieldNodes;	
}

void CostFunction3D::initialize(const real& imin, const real& imax, const int& idim,
									const real& jmin, const real& jmax, const int& jdim,
									const real& kmin, const real& kmax, const int& kdim,
									real* bgU, real* obs)
{

	// Initialize number of variables
	varDim = 6;
	
	// Set the constant height for the XY case
	constHeight = 0.;
	
	// Initialize background errors and filter scales
	bgError[0] = 50.;
	bgError[1] = 50.;
	bgError[2] = 50.;
	bgError[3] = 25.;
	bgError[4] = 10.;	
	bgError[5] = 10.;

	// These should be set by the xml file
	iBCL[0] = R1T2; iBCR[0] = R1T2; jBCL[0] = R1T2; jBCR[0] = R1T2; kBCL[0] = R1T2; kBCR[0] = R1T2;
	iBCL[1] = R1T2; iBCR[1] = R1T2; jBCL[1] = R1T2; jBCR[1] = R1T2; kBCL[1] = R1T2; kBCR[1] = R1T2;
	iBCL[2] = R1T2; iBCR[2] = R1T2; jBCL[2] = R1T2; jBCR[2] = R1T2; kBCL[2] = R1T0; kBCR[2] = R1T0;
	iBCL[3] = R1T2; iBCR[3] = R1T2; jBCL[3] = R1T2; jBCR[3] = R1T2; kBCL[3] = R1T2; kBCR[3] = R1T2;
	iBCL[4] = R1T2; iBCR[4] = R1T2; jBCL[4] = R1T2; jBCR[4] = R1T2;	kBCL[4] = R1T2; kBCR[4] = R1T2;
	iBCL[5] = R1T2; iBCR[5] = R1T2; jBCL[5] = R1T2; jBCR[5] = R1T2; kBCL[5] = R1T2; kBCR[5] = R1T2;
	
	// Also should be set in XML
	referenceState = jordan;
	
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

	//	Mass continuity weight
	mcWeight = 100 * (DI * DJ * DK);
	
	// Set up the initial recursive filter
	real iFilterScale = 2.;
	real jFilterScale = 2.;
	real kFilterScale = 2.;	
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

	// These are local to this one
	CTHTd = new real[nState];
	stateU = new real[nState];
	stateV = new real[nState];
	obsVector = new real[mObs*12];
	int nodes = iDim*jDim*kDim;
	HCq = new real[mObs+nodes];
	innovation = new real[mObs+nodes];	
	fieldNodes = new real[24*nodes];	
	bState = iDim*jDim*kDim*varDim;	
	bgState = new real[bState];
	stateB = new real[bState];
	stateA = new real[bState];
	stateC = new real[bState];
	
	// Set up the spline matrices
	setupSplines();

}	

void CostFunction3D::initState(const real& iFilterScale, const real& jFilterScale, const real& kFilterScale)
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
	
	// Change the filter scale
	iFilter->setFilterLengthScale(iFilterScale);
	jFilter->setFilterLengthScale(jFilterScale);
	kFilter->setFilterLengthScale(kFilterScale);
	
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

	cout << "Beginning analysis...\n";
	
	// HTd
	calcHTranspose(innovation, stateC);
	
	SCtranspose(stateC, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtransform(stateA, stateB);
	
	SBtranspose(stateB, CTHTd);
			
}	

double CostFunction3D::funcValue(double* state)
{

	double qIP, obIP, mcIP;
	qIP = 0.;
	obIP = 0.;
	mcIP = 0.;
	
	updateHCq(state);

	// Compute inner product of state and mass residualvectors
	for (int n = 0; n < nState; n++) {
		qIP += state[n]*state[n];
	}
		
	// Subtract d from HCq to yield mObs length vector and compute inner product
	//#pragma omp parallel for reduction(+:obIP)
	for (int m = 0; m < mObs; m++) {
		obIP += (HCq[m]-innovation[m])*(obsVector[m*12+1])*(HCq[m]-innovation[m]);
	}
	
	// Mass continuity on the nodes
	//#pragma omp parallel for reduction(+:mcIP)
	for (int kIndex = 0; kIndex < kDim; kIndex++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				int hIndex = mObs + iDim*jDim*kIndex + iDim*jIndex + iIndex;
				mcIP += (HCq[hIndex]-innovation[hIndex])*mcWeight*(HCq[hIndex]-innovation[hIndex]);
			}
		}
	}
	
	double J = 0.5*(qIP + obIP + mcIP);
	return J;
	
}

void CostFunction3D::funcGradient(double* state, double* gradient)
{
	
	updateHCq(state);
		
	// HTHCq
	calcHTranspose(HCq, stateC);
	
	SCtranspose(stateC, stateA);
	
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
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
					if ((iNode < 0) or (iNode >= iDim) or 
						(jNode < 0) or (jNode >= jDim) or
						(kNode < 0) or (kNode >= kDim)) continue;
					int cIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
					for (int var = 0; var < varDim; var++) {
						if (obsVector[mi+6 + var] == 0) continue;
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
						kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
						real basis3x = ibasis * jbasis * kbasis;
						tempsum += stateC[cIndex + var] *  basis3x * obsVector[mi+6 + var];
					}
				}
			}
		}
		HCq[m] = tempsum;
	}
	
	#pragma omp parallel for
	for (int kIndex = 0; kIndex < kDim; kIndex++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				real i = iIndex*DI + iMin;
				real j = jIndex*DJ + jMin;
				real k = kIndex*DK + kMin;
				
				int ii = iIndex;
				int jj = jIndex;
				int kk = kIndex;
				int hIndex = mObs + iDim*jDim*kIndex + iDim*jIndex + iIndex;

				real ibasis = 0.;
				real jbasis = 0.;
				real kbasis = 0.;
				real idbasis = 0.;
				real jdbasis = 0.;
				real kdbasis = 0.;
				real tempsum = 0.;
				for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
					for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
						for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
							if ((iNode < 0) or (iNode >= iDim) or 
								(jNode < 0) or (jNode >= jDim) or
								(kNode < 0) or (kNode >= kDim)) continue;
							int cIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
							
							idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, iBCL[0], iBCR[0]);
							jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[0], jBCR[0]);
							kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[0], kBCR[0]);							
							tempsum += (stateC[cIndex] * idbasis * jbasis * kbasis);
							
							ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[1], iBCR[1]);
							jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, jBCL[1], jBCR[1]);
							kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[1], kBCR[1]);							
							tempsum += (stateC[cIndex+1] * ibasis * jdbasis * kbasis);
							
							ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[2], iBCR[2]);
							jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[2], jBCR[2]);
							kdbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 1, kBCL[2], kBCR[2]);							
							tempsum += (stateC[cIndex+2] * ibasis * jbasis * kdbasis);
							
						}
					}
				}
				HCq[hIndex] = tempsum;
			}
		}
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
	
	real innovationRMS = 0.;
	real mcbgRMS = 0.;
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
		real ibasis = 0.;
		real jbasis = 0.;
		real kbasis = 0.;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
					if ((iNode < 0) or (iNode >= iDim) or 
						(jNode < 0) or (jNode >= jDim) or
						(kNode < 0) or (kNode >= kDim)) continue;
					for (int var = 0; var < varDim; var++) {
						if (obsVector[mi+6 + var] == 0) continue;
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
						kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
					    int stateIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
						real basis3x = ibasis * jbasis * kbasis;
						tempsum += bgState[stateIndex + var] * basis3x * obsVector[mi+6+var];
					}
				}
			}
		}
		innovation[m] -= tempsum;
		innovationRMS += (innovation[m]*innovation[m]);
	}
	
    #pragma omp parallel for reduction(+:mcbgRMS)
	for (int kIndex = 0; kIndex < kDim; kIndex++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				real i = iIndex*DI + iMin;
				real j = jIndex*DJ + jMin;
				real k = kIndex*DK + kMin;
				
				int ii = iIndex;
				int jj = jIndex;
				int kk = kIndex;
				int hIndex = mObs + iDim*jDim*kIndex + iDim*jIndex + iIndex;
				
				real ibasis = 0.;
				real jbasis = 0.;
				real kbasis = 0.;
				real idbasis = 0.;
				real jdbasis = 0.;
				real kdbasis = 0.;
				real tempsum = 0.;
				for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
					for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
						for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
							if ((iNode < 0) or (iNode >= iDim) or 
								(jNode < 0) or (jNode >= jDim) or
								(kNode < 0) or (kNode >= kDim)) continue;
							int bgIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;

							idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, iBCL[0], iBCR[0]);
							jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[0], jBCR[0]);
							kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[0], kBCR[0]);							
							tempsum += (bgState[bgIndex] * idbasis * jbasis * kbasis);
							
							ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[1], iBCR[1]);
							jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, jBCL[1], jBCR[1]);
							kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[1], kBCR[1]);							
							tempsum += (bgState[bgIndex + 1] * ibasis * jdbasis * kbasis);
										
							ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[2], iBCR[2]);
							jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[2], jBCR[2]);
							kdbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 1, kBCL[2], kBCR[2]);							
							tempsum += (bgState[bgIndex + 2] * ibasis * jbasis * kdbasis);
						}
					}
				}
				innovation[hIndex] = -tempsum;
				mcbgRMS += (innovation[hIndex]*innovation[hIndex]);
			}
		}
	}
		
	if (mObs) innovationRMS /= mObs;
	innovationRMS = sqrt(innovationRMS);
	cout << "Innovation RMS : " << innovationRMS << endl;
	mcbgRMS /= (iDim*jDim*kDim);
	mcbgRMS = sqrt(mcbgRMS);
	cout << "Background Mass Continuity RMS : " << mcbgRMS << endl;

}	

void CostFunction3D::calcHTranspose(const real* yhat, real* Astate) 
{
	
	// Clear the Astate
	for (int b = 0; b < bState; b++) {
		Astate[b] = 0.;
	}
	
	// Calculate H Transpose	
	#pragma omp parallel for
	for (int kIndex = 0; kIndex < kDim; kIndex++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
			
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
					for (int var = 0; var < varDim; var++) {
						if (obsVector[mi+6 + var] == 0) continue;
						ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
						jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
						kbasis = Basis(kIndex, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
						int aIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex;
						real qbasise = qhat* ibasis * jbasis * kbasis *invError;
						#pragma omp atomic
						Astate[aIndex + var] += qbasise * obsVector[mi+6+var];
					}
				}
				
				// Mass continuity transpose				
				int ii = iIndex;
				int jj = jIndex;
				int kk = kIndex;
				for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
					for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
						for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
							if ((iNode < 0) or (iNode >= iDim) or 
								(jNode < 0) or (jNode >= jDim) or
								(kNode < 0) or (kNode >= kDim)) continue;
							
							real i = iNode*DI + iMin;
							real j = jNode*DJ + jMin;
							real k = kNode*DK + kMin;
							real ibasis = 0.;
							real jbasis = 0.;
							real kbasis = 0.;
							real idbasis = 0.;
							real jdbasis = 0.;
							real kdbasis = 0.;
							
							int hIndex = mObs + iDim*jDim*kNode + iDim*jNode + iNode;
							int aIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex;
							
							idbasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 1, iBCL[0], iBCR[0]);
							jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[0], jBCR[0]);
							kbasis = Basis(kIndex, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[0], kBCR[0]);							
							#pragma omp atomic
							Astate[aIndex] += mcWeight * (yhat[hIndex] * idbasis * jbasis * kbasis);
							
							ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[1], iBCR[1]);
							jdbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 1, jBCL[1], jBCR[1]);
							kbasis = Basis(kIndex, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[1], kBCR[1]);							
							#pragma omp atomic
							Astate[aIndex + 1] += mcWeight * (yhat[hIndex] * ibasis * jdbasis * kbasis);
							
							ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[2], iBCR[2]);
							jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[2], jBCR[2]);
							kdbasis = Basis(kIndex, k, kDim-1, kMin, DK, DKrecip, 1, kBCL[2], kBCR[2]);							
							#pragma omp atomic
							Astate[aIndex + 2] += mcWeight * (yhat[hIndex] * ibasis * jbasis * kdbasis);
						}		
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
	for (int kIndex = 0; kIndex < (kDim-1); kIndex++) {
		for (int kmu = -1; kmu <= 1; kmu += 2) {
			real k = kMin + DK * (kIndex + (0.5*sqrt(1./3.) * kmu + 0.5));
			int kk = (int)((k - kMin)*DKrecip);
			for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
				if ((kNode < 0) or (kNode >= kDim)) continue;
		
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
										
										int uI = iIndex*2 + (imu+1)/2;
										int uJ = jIndex*2 + (jmu+1)/2;
										int uK = kIndex*2 + (kmu+1)/2;
										int uIndex = varDim*(iDim-1)*2*(jDim-1)*2*uK +varDim*(iDim-1)*2*uJ +varDim*uI;
										int bIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
										for (int var = 0; var < varDim; var++) {
											int ui = uIndex + var;
											if (Ustate[ui] == 0) continue;
											int bi = bIndex + var;
											real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
											real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
											real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
											#pragma omp atomic
											Bstate[bi] += 0.125 * Ustate[ui] * ibasis * jbasis * kbasis; 
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
	for (int n = 0; n < nState; n++) {
		Ustate[n] = 0;
	}
	
	#pragma omp parallel for
	for (int kIndex = 0; kIndex < (kDim-1); kIndex++) {
		for (int kmu = -1; kmu <= 1; kmu += 2) {
			real k = kMin + DK * (kIndex + (0.5*sqrt(1./3.) * kmu + 0.5));
			int kk = (int)((k - kMin)*DKrecip);
			for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
				if ((kNode < 0) or (kNode >= kDim)) continue;
		
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
										int uI = iIndex*2 + (imu+1)/2;
										int uJ = jIndex*2 + (jmu+1)/2;
										int uK = kIndex*2 + (kmu+1)/2;
										int uIndex = varDim*(iDim-1)*2*(jDim-1)*2*uK +varDim*(iDim-1)*2*uJ +varDim*uI;
										int bIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
										for (int var = 0; var < varDim; var++) {
											int bi = bIndex + var;
											if (Bstate[bi] == 0) continue;
											int ui = uIndex + var;
											real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
											real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
											real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
                                            #pragma omp atomic
											Ustate[ui] += 0.125 * Bstate[bi] * ibasis * jbasis * kbasis; 
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
	#pragma omp parallel for
	for (int var = 0; var < varDim; var++) {
		// These are local for parallelization
		real* iTemp = new real[iDim];
		real* jTemp = new real[jDim];
		real* kTemp = new real[kDim];
		
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
		delete[] iTemp;
		delete[] jTemp;
		delete[] kTemp;
	}
}

void CostFunction3D::SCtranspose(real* Cstate, real* Astate)
{
	
	// Isoptropic Recursive filter for speed, no anisotropic "triad" working yet 
	#pragma omp parallel for
	for (int var = 0; var < varDim; var++) {
		// These are local for parallelization
		real* iTemp = new real[iDim];
		real* jTemp = new real[jDim];
		real* kTemp = new real[kDim];
		
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
		delete[] iTemp;
		delete[] jTemp;
		delete[] kTemp;
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
	
	real cutoff_wl = 2;
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
					real pm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
					real qm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 3, iBCL[var], iBCR[var]);
					real pn, qn;
					P[iNode][iNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((iNode+1) < iDim) {
						pn = Basis(iNode+1, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
						qn = Basis(iNode+1, i, iDim-1, iMin, DI, DIrecip, 3, iBCL[var], iBCR[var]);
						P[iNode][iNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+1][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((iNode+2) < iDim) {
						pn = Basis(iNode+2, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
						qn = Basis(iNode+2, i, iDim-1, iMin, DI, DIrecip, 3, iBCL[var], iBCR[var]);
						P[iNode][iNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+2][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((iNode+3) < iDim) {
						pn = Basis(iNode+3, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
						qn = Basis(iNode+3, i, iDim-1, iMin, DI, DIrecip, 3, iBCL[var], iBCR[var]);
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
	
	cutoff_wl = 2;
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
					real pm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
					real qm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 3, jBCL[var], jBCR[var]);
					real pn, qn;
					P[jNode][jNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((jNode+1) < jDim) {
						pn = Basis(jNode+1, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
						qn = Basis(jNode+1, j, jDim-1, jMin, DJ, DJrecip, 3, jBCL[var], jBCR[var]);
						P[jNode][jNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+1][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((jNode+2) < jDim) {
						pn = Basis(jNode+2, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
						qn = Basis(jNode+2, j, jDim-1, jMin, DJ, DJrecip, 3, jBCL[var], jBCR[var]);
						P[jNode][jNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+2][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));						
					}
					if ((jNode+3) < jDim) {
						pn = Basis(jNode+3, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
						qn = Basis(jNode+3, j, jDim-1, jMin, DJ, DJrecip, 3, jBCL[var], jBCR[var]);
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
	
	cutoff_wl = 2;
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
					real pm = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
					real qm = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 3, kBCL[var], kBCR[var]);
					real pn, qn;
					P[kNode][kNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((kNode+1) < kDim) {
						pn = Basis(kNode+1, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
						qn = Basis(kNode+1, k, kDim-1, kMin, DK, DKrecip, 3, kBCL[var], kBCR[var]);
						P[kNode][kNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[kNode+1][kNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((kNode+2) < kDim) {
						pn = Basis(kNode+2, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
						qn = Basis(kNode+2, k, kDim-1, kMin, DK, DKrecip, 3, kBCL[var], kBCR[var]);
						P[kNode][kNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[kNode+2][kNode] += 0.5 * ((pm * pn) + eq * (qm * qn));						
					}
					if ((kNode+3) < kDim) {
						pn = Basis(kNode+3, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
						qn = Basis(kNode+3, k, kDim-1, kMin, DK, DKrecip, 3, kBCL[var], kBCR[var]);
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
	
	cout << "Outputting " << suffix.toStdString() << "...\n";
	// H --> to Mish for output
	QString samuraiout = "samurai_XYZ_" + suffix + ".out";
	ofstream samuraistream(samuraiout.toAscii().data());
	samuraistream << "X\tY\tZ\trhoE\tu\tv\tw\tVorticity\tDivergence\tqv'\trho'\tT'\tP'\th\t";
	samuraistream << "rhoux\trhouy\trhouz\trhovx\trhovy\trhovz\trhowx\trhowy\trhowz\tMC residual\n";
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
							
							real tpw = 0;
							
							for (int kIndex = 0; kIndex < kDim; kIndex++) {
								for (int khalf =0; khalf <=1; khalf++) {
									for (int kmu = -khalf; kmu <= khalf; kmu++) {
										real k = kMin + DK * (kIndex + (0.5*sqrt(1./3.) * kmu + 0.5*khalf));
										if (k > ((kDim-1)*DK + kMin)) continue;	
										
										real heightm = 1000*k;
										real rhoBar = getReferenceVariable(rhoaref, heightm);
										real qBar = getReferenceVariable(qvbhypref, heightm);
										real hBar = getReferenceVariable(href, heightm);

										int ii = (int)((i - iMin)*DIrecip);
										int jj = (int)((j - jMin)*DJrecip);
										int kk = (int)((k - kMin)*DKrecip);
										real ibasis = 0.;
										real jbasis = 0.;
										real kbasis = 0.;
										real idbasis = 0.;
										real jdbasis = 0.;
										real kdbasis = 0.;
										real rhov = 0.;
										real rhou = 0.;
										real rhow = 0.;
										real rhovdx = 0.;
										real rhoudx = 0.;
										real rhowdx = 0.;
										real rhovdy = 0.;
										real rhoudy = 0.;
										real rhowdy = 0.;
										real rhovdz = 0.;
										real rhoudz = 0.;
										real rhowdz = 0.;
										real hprime = 0.;
										real qvprime = 0.;
										real rhoprime = 0.;
										for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
											for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
												for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
													if ((iNode < 0) or (iNode >= iDim) or 
														(jNode < 0) or (jNode >= jDim) or
														(kNode < 0) or (kNode >= kDim)) continue;
													for (int var = 0; var < varDim; var++) {
														ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
														jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
														kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
														idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, iBCL[var], iBCR[var]);
														jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, jBCL[var], jBCR[var]);
														kdbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 1, kBCL[var], kBCR[var]);
														real basis3x = ibasis*jbasis*kbasis;
														int aIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
														switch (var) {
															case 0:
																rhou +=  Astate[aIndex] * basis3x;
																rhoudx += Astate[aIndex] * idbasis * jbasis * kbasis;
																rhoudy += Astate[aIndex] * ibasis * jdbasis * kbasis;
																rhoudz += Astate[aIndex] * ibasis * jbasis * kdbasis;
																break;
															case 1:
																rhov +=  Astate[aIndex + 1] * basis3x;
																rhovdx += Astate[aIndex + 1] * idbasis * jbasis * kbasis;
																rhovdy += Astate[aIndex + 1] * ibasis * jdbasis * kbasis;
																rhovdz += Astate[aIndex + 1] * ibasis * jbasis * kdbasis;
																break;
															case 2:
																rhow += Astate[aIndex + 2] * basis3x;
																rhowdx += Astate[aIndex + 2] * idbasis * jbasis * kbasis;
																rhowdy += Astate[aIndex + 2] * ibasis * jdbasis * kbasis;
																rhowdz += Astate[aIndex + 2] * ibasis * jbasis * kdbasis;
																break;
															case 3:
																hprime += Astate[aIndex + 3] * basis3x;
																break;
															case 4:
																qvprime += Astate[aIndex + 4] * basis3x;
																break;
															case 5:
																rhoprime += Astate[aIndex + 5] * basis3x;
																break;
														}
													}
												}
											}
										}
										
										// Get it to relevant variables for the flux calculation
										real rhoa = rhoBar + rhoprime / 100;
										real qv = bhypInvTransform(qBar + qvprime);
										real rhoq = qv * rhoa / 1000.;
										real rho = rhoa + rhoq;
										real v = rhov / rho;
										real u = rhou / rho;
										real w = rhow / rho;
										real wspd = sqrt(u*u + v*v);
										real KE = 0.5*rho*(v*v + u*u + w*w);
										real rhoE = rho*(hBar + hprime*1.e3) + KE;
										real h = hBar + hprime*1.e3;
										real temp = (h - 2.501e3*(qv) - 9.81*heightm)/1005.7;
										real airpress = temp*rhoa*287./100.;
										real satvp = 6.1078 * exp(5.0065 * log(273.15/temp)) * exp((5.0065 + 19.83923) * (1 - 273.15/temp));
										real satqv = 1000.0 * 0.622 * satvp / airpress;
										real relhum = -999.;
										if ((satqv != 0) and (satqv >= qv))
											relhum = 100*qv/satqv;
										real vp = temp*rhoq*461./100.;
										real press = airpress + vp;
										
										// Vorticity units are 10-5
										real vorticity = 100*(rhovdx - rhoudy)/rho;
										real divergence = 100*(rhoudx + rhovdy)/rho;
										real s1 = 100*(rhoudx - rhovdy)/rho;
										real s2 = 100*(rhovdx + rhoudy)/rho;
										real strain = sqrt(s1*s1 + s2*s2);
										real okuboweiss = vorticity*vorticity - s1*s1 -s2*s2;
										real mcresidual = rhoudx + rhovdy + rhowdz;
										
										// Output it
										if (updateMish and imu and jmu and kmu and ihalf and jhalf and khalf) {
											int bgI = iIndex*2 + (imu+1)/2;
											int bgJ = jIndex*2 + (jmu+1)/2;
											int bgK = kIndex*2 + (kmu+1)/2;
											int bIndex = varDim*(iDim-1)*2*(jDim-1)*2*bgK + varDim*(iDim-1)*2*bgJ +varDim*bgI;
											bgFields[bIndex + 0] = rhou;
											bgFields[bIndex + 1] = rhov;
											bgFields[bIndex + 2] = rhow;
											bgFields[bIndex + 3] = hprime;
											bgFields[bIndex + 4] = qvprime;
											bgFields[bIndex + 5] = rhoprime;
										}
										
										real pprime = press - getReferenceVariable(pressref, heightm)/100.;
										real tprime = temp - getReferenceVariable(tempref, heightm);
										
										if (!ihalf and !jhalf and !khalf){
											// On the node
											samuraistream << scientific << i << "\t" << j << "\t"  << k << "\t" << rhoE
											<< "\t" << u << "\t" << v << "\t" << w << "\t" << vorticity << "\t" << divergence
											<< "\t" << qvprime*2 << "\t" << rhoprime << "\t" << tprime << "\t" << pprime <<  "\t" << hprime << "\t"
											<< rhoudx << "\t" << rhoudy << "\t" << rhoudz << "\t"
											<< rhovdx << "\t" << rhovdy << "\t" << rhovdz << "\t"
											<< rhowdx << "\t" << rhowdy << "\t" << rhowdz << "\t" << mcresidual << "\n";

											// Sum up the TPW in the vertical, top level is tpw
											tpw += qv * rhoa * DK;
											
											int fIndex = iDim*jDim*kDim; 
											int posIndex = iDim*jDim*kIndex + iDim*jIndex + iIndex;
											fieldNodes[fIndex * 0 + posIndex] = u;
											fieldNodes[fIndex * 1 + posIndex] = v;
											fieldNodes[fIndex * 2 + posIndex] = w;
											fieldNodes[fIndex * 3 + posIndex] = wspd;
											fieldNodes[fIndex * 4 + posIndex] = relhum;
											fieldNodes[fIndex * 5 + posIndex] = hprime;
											fieldNodes[fIndex * 6 + posIndex] = qvprime*2;
											fieldNodes[fIndex * 7 + posIndex] = rhoprime;
											fieldNodes[fIndex * 8 + posIndex] = tprime;
											fieldNodes[fIndex * 9 + posIndex] = pprime;
											fieldNodes[fIndex * 10 + posIndex] = vorticity;
											fieldNodes[fIndex * 11 + posIndex] = divergence;
											fieldNodes[fIndex * 12 + posIndex] = okuboweiss;
											fieldNodes[fIndex * 13 + posIndex] = strain;
											fieldNodes[fIndex * 14 + posIndex] = tpw;
											fieldNodes[fIndex * 15 + posIndex] = rhou;
											fieldNodes[fIndex * 16 + posIndex] = rhov;
											fieldNodes[fIndex * 17 + posIndex] = rhow;
											fieldNodes[fIndex * 18 + posIndex] = rho;
											fieldNodes[fIndex * 19 + posIndex] = press;
											fieldNodes[fIndex * 20 + posIndex] = temp;
											fieldNodes[fIndex * 21 + posIndex] = qv;
											fieldNodes[fIndex * 22 + posIndex] = h;
											fieldNodes[fIndex * 23 + posIndex] = rhowdz;
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

	QString fileName = "samurai_XYZ_" + suffix;
	QString outFileName;
	if(QDir::isAbsolutePath(fileName)) {
		outFileName = fileName;
	}
	else {
		outFileName = QDir::current().filePath(fileName);
	}
		
	QString cdfFileName = outFileName + ".nc";
	if (!writeNetCDF(cdfFileName))
		cout << "Error writing netcdf file " << cdfFileName.toStdString() << endl; 	
	
	// Write the Obs to a summary text file
	QString qcout = "samurai_QC_" + suffix + ".out";
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
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
					if ((iNode < 0) or (iNode >= iDim) or 
						(jNode < 0) or (jNode >= jDim) or
						(kNode < 0) or (kNode >= kDim)) continue;
					int cIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
					for (int var = 0; var < varDim; var++) {
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
						kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
						real basis3x = ibasis * jbasis * kbasis;
						tempsum += stateC[cIndex + var] * basis3x * obsVector[mi+6+var];
					}
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
	
	// Write out the CAPPI to an asi file
	// Initialize header
	int id[511];
	for (int n = 1; n <= 510; n++) {
		id[n]=-999;
	}
	
	// Calculate headers
	QStringList fieldNames;
	fieldNames  << "U" << "V" << "W" << "WS" << "RH"<< "HP" << "QP" << "RP" << "TP" << "PP" << "VO" << "DV" << "OW" << "S" << "PW";
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
	QString asiFileName = outFileName + ".asi";
	QFile asiFile(asiFileName);
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
					fieldNodes[iDim*jDim*kDim*n +iDim*jDim*k + iDim*j + i];
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
	
	return true;
	
}     

bool CostFunction3D::writeNetCDF(const QString& netcdfFile)
{
	NcError err(NcError::verbose_nonfatal);
	int NC_ERR = 0;
	
	// Create the file.
	NcFile dataFile(netcdfFile.toAscii(), NcFile::Replace);
	
	// Check to see if the file was created.
	if(!dataFile.is_valid())
		return NC_ERR;
	
	// Define the dimensions. NetCDF will hand back an ncDim object for
	// each.
	NcDim *lvlDim, *latDim, *lonDim, *timeDim;
	if (!(lvlDim = dataFile.add_dim("altitude", kDim)))
		return NC_ERR;
	if (!(latDim = dataFile.add_dim("latitude", jDim)))
		return NC_ERR;
	if (!(lonDim = dataFile.add_dim("longitude", iDim)))
		return NC_ERR;
	// Add an unlimited dimension...
	if (!(timeDim = dataFile.add_dim("time")))
		return NC_ERR;
	
	// Define the coordinate variables.
	NcVar *latVar, *lonVar, *lvlVar, *timeVar;
	if (!(latVar = dataFile.add_var("latitude", ncFloat, latDim)))
		return NC_ERR;
	if (!(lonVar = dataFile.add_var("longitude", ncFloat, lonDim)))
		return NC_ERR;
	if (!(lvlVar = dataFile.add_var("altitude", ncFloat, lvlDim)))
		return NC_ERR;
	if (!(timeVar = dataFile.add_var("time", ncFloat, timeDim)))
		return NC_ERR;
	
	// Define units attributes for coordinate vars. This attaches a
	// text attribute to each of the coordinate variables, containing
	// the units.
	if (!latVar->add_att("units", "degrees_north"))
		return NC_ERR;
	if (!lonVar->add_att("units", "degrees_east"))
		return NC_ERR;
	if (!lvlVar->add_att("units", "km"))
		return NC_ERR;
	if (!timeVar->add_att("units", "seconds since 1970-01-01 00:00:00 +0000"))
		return NC_ERR;
	
	// Define the netCDF variables 
	NcVar *u, *v, *w, *wspd, *relhum, *hprime, *qvprime, *rhoprime, *tprime, *pprime;
	NcVar *vorticity, *divergence, *okuboweiss, *strain, *tpw, *rhou, *rhov, *rhow;
	NcVar *rho, *press, *temp, *qv, *h, *rhowdz;
		
	if (!(u = dataFile.add_var("U", ncFloat, timeDim, 
									 lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(v = dataFile.add_var("V", ncFloat, timeDim, 
							   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(w = dataFile.add_var("W", ncFloat, timeDim, 
							   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(wspd = dataFile.add_var("WSPD", ncFloat, timeDim, 
							   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(relhum = dataFile.add_var("RH", ncFloat, timeDim, 
							   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(hprime = dataFile.add_var("HP", ncFloat, timeDim, 
							   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(qvprime = dataFile.add_var("QVP", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rhoprime = dataFile.add_var("RHOAP", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(tprime = dataFile.add_var("TP", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(pprime = dataFile.add_var("PP", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(vorticity = dataFile.add_var("VORT", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(divergence = dataFile.add_var("DIV", ncFloat, timeDim, 
									   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(okuboweiss = dataFile.add_var("OW", ncFloat, timeDim, 
									   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(strain = dataFile.add_var("STRAIN", ncFloat, timeDim, 
									 lvlDim, latDim, lonDim)))
		return NC_ERR;       
	if (!(tpw = dataFile.add_var("TPW", ncFloat, timeDim, 
								lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rhou = dataFile.add_var("RHOU", ncFloat, timeDim, 
							   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rhov = dataFile.add_var("RHOV", ncFloat, timeDim, 
							   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rhow = dataFile.add_var("RHOW", ncFloat, timeDim, 
							   lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rho = dataFile.add_var("RHOA", ncFloat, timeDim, 
									  lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(press = dataFile.add_var("P", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(temp = dataFile.add_var("T", ncFloat, timeDim, 
								  lvlDim, latDim, lonDim)))
		return NC_ERR;	
	if (!(qv = dataFile.add_var("QV", ncFloat, timeDim, 
								lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(h = dataFile.add_var("H", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(rhowdz = dataFile.add_var("DRHOWDZ", ncFloat, timeDim, 
								  lvlDim, latDim, lonDim)))
		return NC_ERR;
	
	
	// Define units attributes for data variables.
	if (!u->add_att("units", "m s-1"))
		return NC_ERR;
	if (!v->add_att("units", "m s-1"))
		return NC_ERR;
	if (!w->add_att("units", "m s-1"))
		return NC_ERR;
	if (!wspd->add_att("units", "m s-1"))
		return NC_ERR;
	if (!relhum->add_att("units", "percent"))
		return NC_ERR;
	if (!hprime->add_att("units", "kJ"))
		return NC_ERR;
	if (!qvprime->add_att("units", "g kg-1")) 
		return NC_ERR;
	if (!rhoprime->add_att("units", "kg m-3")) 
		return NC_ERR;
	if (!tprime->add_att("units", "K")) 
		return NC_ERR;
	if (!pprime->add_att("units", "hPa")) 
		return NC_ERR;
	if (!vorticity->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!divergence->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!okuboweiss->add_att("units", "10-10s-1")) 
		return NC_ERR;
	if (!strain->add_att("units", "10-5s-1")) 
		return NC_ERR;       
	if (!tpw->add_att("units", "mm")) 
		return NC_ERR;
	if (!rhou->add_att("units", "kg m-2s-1")) 
		return NC_ERR;
	if (!rhov->add_att("units", "kg m-2s-1")) 
		return NC_ERR;
	if (!rhow->add_att("units", "kg m-2s-1")) 
		return NC_ERR;
	if (!rho->add_att("units", "kg m-2s-1")) 
		return NC_ERR;
	if (!press->add_att("units", "hPa")) 
		return NC_ERR;
	if (!temp->add_att("units", "K")) 
		return NC_ERR;	
	if (!qv->add_att("units", "g kg-1")) 
		return NC_ERR;
	if (!h->add_att("units", "kJ")) 
		return NC_ERR;
	if (!rhowdz->add_att("units", "kg m-3s-1")) 
		return NC_ERR;
	
	// Define long names for data variables.
	if (!u->add_att("long_name", "u wind component"))
		return NC_ERR;
	if (!v->add_att("long_name", "v wind component"))
		return NC_ERR;
	if (!w->add_att("long_name", "w wind component"))
		return NC_ERR;
	if (!wspd->add_att("long_name", "wind speed"))
		return NC_ERR;
	if (!relhum->add_att("long_name", "relative humidity"))
		return NC_ERR;
	if (!hprime->add_att("long_name", "moist static energy perturbation"))
		return NC_ERR;
	if (!qvprime->add_att("long_name", "water vapor mixing ratio perturbation")) 
		return NC_ERR;
	if (!rhoprime->add_att("long_name", "air density perturbation")) 
		return NC_ERR;
	if (!tprime->add_att("long_name", "temperature perturbation")) 
		return NC_ERR;
	if (!pprime->add_att("long_name", "pressure perturbation")) 
		return NC_ERR;
	if (!vorticity->add_att("long_name", "vertical vorticity")) 
		return NC_ERR;
	if (!divergence->add_att("long_name", "horizontal divergence")) 
		return NC_ERR;
	if (!okuboweiss->add_att("long_name", "Okubo-Weiss parameter")) 
		return NC_ERR;
	if (!strain->add_att("long_name", "horizontal strain")) 
		return NC_ERR;       
	if (!tpw->add_att("long_name", "total precipitable water")) 
		return NC_ERR;
	if (!rhou->add_att("long_name", "mass-weighted u wind component")) 
		return NC_ERR;
	if (!rhov->add_att("long_name", "mass-weighted v wind component")) 
		return NC_ERR;
	if (!rhow->add_att("long_name", "mass-weighted w wind component")) 
		return NC_ERR;
	if (!rho->add_att("long_name", "density")) 
		return NC_ERR;
	if (!press->add_att("long_name", "pressure")) 
		return NC_ERR;
	if (!temp->add_att("long_name", "temperature")) 
		return NC_ERR;	
	if (!qv->add_att("long_name", "water vapor mixing ratio")) 
		return NC_ERR;
	if (!h->add_att("long_name", "moist static energy")) 
		return NC_ERR;
	if (!rhowdz->add_att("long_name", "vertical mass flux gradient")) 
		return NC_ERR;	
	
	// Write the coordinate variable data to the file.
	real *lats = new real[iDim];
	real *lons = new real[jDim];
	real *levs = new real[kDim];
	int time[2];
	
	// Hard code this for now, but it needs to be dynamic
	time[0] = 1282046400;
	real latReference = 15.400;
	real lonReference = -68.603;
	double latrad =latReference * 1.745329251994e-02;
	double fac_lat = 111.13209 - 0.56605 * cos(2.0 * latrad)
	+ 0.00012 * cos(4.0 * latrad) - 0.000002 * cos(6.0 * latrad);
	double fac_lon = 111.41513 * cos(latrad)
	- 0.09455 * cos(3.0 * latrad) + 0.00012 * cos(5.0 * latrad);
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		real i = iMin + DI * iIndex;
		lons[iIndex] = lonReference + i/fac_lon;
	}
	if (!lonVar->put(lons, iDim))
		return NC_ERR;       

	for (int jIndex = 0; jIndex < jDim; jIndex++) {
		real j = jMin + DJ * jIndex;
		lats[jIndex] = latReference + j/fac_lat;
	}
	if (!latVar->put(lats, jDim))
		return NC_ERR;       
	
	for (int kIndex = 0; kIndex < kDim; kIndex++) {
		real k = kMin + DK * kIndex;
		levs[kIndex] = k;
	}
	if (!lvlVar->put(levs, kDim))
		return NC_ERR; 
	
	if (!timeVar->put(time, 1))
		return NC_ERR;
	
	// Write the data.
	for (int rec = 0; rec < 1; rec++) 
	{
		if (!u->put_rec(&fieldNodes[0], rec))
			return NC_ERR;
		if (!v->put_rec(&fieldNodes[iDim*jDim*kDim*1], rec))
			return NC_ERR;
		if (!w->put_rec(&fieldNodes[iDim*jDim*kDim*2], rec))
			return NC_ERR;
		if (!wspd->put_rec(&fieldNodes[iDim*jDim*kDim*3], rec))
			return NC_ERR;
		if (!relhum->put_rec(&fieldNodes[iDim*jDim*kDim*4], rec))
			return NC_ERR;
		if (!hprime->put_rec(&fieldNodes[iDim*jDim*kDim*5], rec))
			return NC_ERR;
		if (!qvprime->put_rec(&fieldNodes[iDim*jDim*kDim*6], rec))
			return NC_ERR;
		if (!rhoprime->put_rec(&fieldNodes[iDim*jDim*kDim*7], rec)) 
			return NC_ERR;
		if (!tprime->put_rec(&fieldNodes[iDim*jDim*kDim*8], rec)) 
			return NC_ERR;
		if (!pprime->put_rec(&fieldNodes[iDim*jDim*kDim*9], rec))
			return NC_ERR;
		if (!vorticity->put_rec(&fieldNodes[iDim*jDim*kDim*10], rec))
			return NC_ERR;
		if (!divergence->put_rec(&fieldNodes[iDim*jDim*kDim*11], rec)) 
			return NC_ERR;
		if (!okuboweiss->put_rec(&fieldNodes[iDim*jDim*kDim*12], rec)) 
			return NC_ERR;
		if (!strain->put_rec(&fieldNodes[iDim*jDim*kDim*13], rec)) 
			return NC_ERR;       
		if (!tpw->put_rec(&fieldNodes[iDim*jDim*kDim*14], rec)) 
			return NC_ERR;
		if (!rhou->put_rec(&fieldNodes[iDim*jDim*kDim*15], rec)) 
			return NC_ERR;
		if (!rhov->put_rec(&fieldNodes[iDim*jDim*kDim*16], rec)) 
			return NC_ERR;
		if (!rhow->put_rec(&fieldNodes[iDim*jDim*kDim*17], rec)) 
			return NC_ERR;
		if (!rho->put_rec(&fieldNodes[iDim*jDim*kDim*18], rec)) 
			return NC_ERR;
		if (!press->put_rec(&fieldNodes[iDim*jDim*kDim*19], rec)) 
			return NC_ERR;
		if (!temp->put_rec(&fieldNodes[iDim*jDim*kDim*20], rec))
			return NC_ERR;	
		if (!qv->put_rec(&fieldNodes[iDim*jDim*kDim*21], rec)) 
			return NC_ERR;
		if (!h->put_rec(&fieldNodes[iDim*jDim*kDim*22], rec)) 
			return NC_ERR;
		if (!rhowdz->put_rec(&fieldNodes[iDim*jDim*kDim*23], rec)) 
			return NC_ERR;
	}
	
	// The file is automatically closed by the destructor. This frees
	// up any internal netCDF resources associated with the file, and
	// flushes any buffers.
	delete[] lats;
	delete[] lons;
	delete[] levs;
	
	return 1;
	
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
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[4], iBCR[4]);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[4], jBCR[4]);
					kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[4], kBCR[4]);
					qvprime += bgState[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * kbasis;					
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[5], iBCR[5]);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[5], jBCR[5]);
					kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[5], kBCR[5]);
					rhoprime += bgState[varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode + 5] * ibasis * jbasis * kbasis;
				}
			}
		}
		real heightm = 1000*k;
		real rhoBar = getReferenceVariable(rhoaref, heightm);
		real qBar = getReferenceVariable(qvbhypref, heightm);
		real qv = bhypInvTransform(qBar + qvprime);
		real rhoa = rhoBar + rhoprime / 100;
		real rhoBG = rhoa + rhoa*qv/1000.;
		
		// Adjust the relevant fields
		if (type == MetObs::sfmr) {
			obsVector[mi] *= rhoBG;
			// Still need to think carefully about nonlinear SFMR operator in Cartesian space
			//real wsBG = vBG*vBG; // + uBG*uBG;
			//obsVector[mi+2] = 2.*vBG/wsBG;
			//obsVector[mi+3] = 0.; //2.*uBG/wsBG;
		}
		if ((type == MetObs::radar) or (type == MetObs::qscat) 
			or (type == MetObs::ascat) or (type == MetObs::AMV)) { 
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
				z = 2.0 - z;
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

real CostFunction3D::getReferenceVariable(int refVariable, real heightm)
{
	real qvbhypcoeff[5];
	real rhoacoeff[5];
	real dpdzcoeff[5];
	
	if (referenceState == jordan) {
		qvbhypcoeff[0] = 9.5108;
		qvbhypcoeff[1] = -0.002761;
		qvbhypcoeff[2] = 3.0362e-7;
		qvbhypcoeff[3] = -1.476e-11;
		qvbhypcoeff[4] = 2.6515e-16;
		
		rhoacoeff[0] = 1.1451;
		rhoacoeff[1] = -0.00010098;
		rhoacoeff[2] = 3.1546e-09;
		rhoacoeff[3] = -2.6251e-14;
		rhoacoeff[4] = -5.5556e-19;
		
		dpdzcoeff[0] = -11.432;
		dpdzcoeff[1] = 0.0010267;
		dpdzcoeff[2] = -2.7622e-08;
		dpdzcoeff[3] = -5.2937e-13;
		dpdzcoeff[4] = 3.4713e-17;
		
	}
	
	if (refVariable == qvbhypref) {
		real qvbhyp = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i); 
			qvbhyp += qvbhypcoeff[i] * power;
		}
		if (qvbhyp < 0.) qvbhyp = 0.;
		return qvbhyp;
	} else if (refVariable == rhoaref) {
		real rhoa = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i); 
			rhoa += rhoacoeff[i] * power;
		}
		return rhoa;
	} else if (refVariable == rhoref) {
		real rho = 0.;
		real qvbhyp = 0.;
		real rhoa = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i); 
			rhoa += rhoacoeff[i] * power;
			qvbhyp += qvbhypcoeff[i] * power;
		}
		if (qvbhyp < 0.) qvbhyp = 0.;
		real qv = bhypInvTransform(qvbhyp);
		rho = rhoa*qv/1000. + rhoa;
		return rho;
	} else if ((refVariable == href) or (refVariable = tempref) or (refVariable == pressref)) {
		// Integrate hydrostatic equation to get pressure and/or solve for T or h
		real press = 0.;
		real temp = 0.;
		real rho = 0.;
		real qvbhyp = 0.;
		real rhoa = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i);
			real power1 = pow(heightm, i+1);
			press += dpdzcoeff[i] * power1 / (i+1);
			rhoa += rhoacoeff[i] * power;
			qvbhyp += qvbhypcoeff[i] * power;
		}
		if (qvbhyp < 0.) qvbhyp = 0.;
		real qv = bhypInvTransform(qvbhyp);
		rho = rhoa*qv/1000. + rhoa;
		press += 101510.0;
		temp = press/(286.9*rhoa + 461.5*rhoa*qv/1000.);
		real h = 1005.7*temp + 9.81*heightm + 2.5e3*qv;
		switch (refVariable) {
			case href:
				return h;
			case tempref:
				return temp;
			case pressref:
				return press;
			default:
				break;
		}
	}
	
	return 0;
}

real CostFunction3D::bhypTransform(real qv)
{
	
	real qvbhyp = 0.5*((qv + 1.e-7) - 1.e-14/(qv + 1.e-7));
	return qvbhyp;
	
}

real CostFunction3D::bhypInvTransform(real qvbhyp)
{
	real qv = 0.;
	if (qvbhyp > 0) {
		qv = sqrt(qvbhyp*qvbhyp + 1.e-14) + qvbhyp - 1.e-7;
	}
	return qv;
}

