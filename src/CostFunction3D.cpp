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
#include <GeographicLib/TransverseMercatorExact.hpp>

CostFunction3D::CostFunction3D(const int& numObs, const int& stateSize)
	: CostFunction(numObs, stateSize)
{
	// Set up the boundary condition hash
	bcHash["R0"] = R0;
	bcHash["R1T0"] = R1T0;
	bcHash["R1T1"] = R1T1;
	bcHash["R1T2"] = R1T2;
	bcHash["R1T10"] = R1T10;
	bcHash["R2T10"] = R2T10;
	bcHash["R2T20"] = R2T20;
	bcHash["R3"] = R3;
	bcHash["PERIODIC"] = PERIODIC;	
    
    // Set the derivative array
    derivative[0][0] = 0;
    derivative[0][1] = 0;
    derivative[0][2] = 0;
    
    derivative[1][0] = 1;
    derivative[1][1] = 0;
    derivative[1][2] = 0;

    derivative[2][0] = 0;
    derivative[2][1] = 1;
    derivative[2][2] = 0;

    derivative[3][0] = 0;
    derivative[3][1] = 0;
    derivative[3][2] = 1;
    
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
	delete[] tempState;
	delete[] tempGradient;
	delete[] xt;
	delete[] df;
	delete[] CTHTd;
	delete[] stateU;
	delete[] obsVector;
	delete[] HCq;
	delete[] innovation;
	delete[] finalAnalysis;
	delete[] bgState;
	delete[] bgStdDev;
	delete[] stateA;
	delete[] stateB;
	delete[] stateC;
	delete[] iL;
	delete[] jL;
	delete[] kL;
	delete[] basis0;
	delete[] basis1;
	
}

void CostFunction3D::initialize(const QHash<QString, QString>* config, real* bgU, real* obs, ReferenceState* ref)
{

	// Initialize number of variables
	varDim = 7;
    derivDim = 4;
	configHash = config;

	// Horizontal boundary conditions
	int ibc = bcHash.value(configHash->value("i_bc"));	
	iBCL[0] = ibc; iBCR[0] = ibc; 
	iBCL[1] = ibc; iBCR[1] = ibc; 
	iBCL[2] = ibc; iBCR[2] = ibc; 
	iBCL[3] = ibc; iBCR[3] = ibc; 
	iBCL[4] = ibc; iBCR[4] = ibc; 
	iBCL[5] = ibc; iBCR[5] = ibc; 
	iBCL[6] = ibc; iBCR[6] = ibc; 
	
	int jbc = bcHash.value(configHash->value("j_bc"));	
	jBCL[0] = jbc; jBCR[0] = jbc; 
	jBCL[1] = jbc; jBCR[1] = jbc;
	jBCL[2] = jbc; jBCR[2] = jbc;
	jBCL[3] = jbc; jBCR[3] = jbc;
	jBCL[4] = jbc; jBCR[4] = jbc;
	jBCL[5] = jbc; jBCR[5] = jbc;
	jBCL[6] = jbc; jBCR[6] = jbc;
	
	int kbc = bcHash.value(configHash->value("k_bc"));
	kBCL[0] = kbc; kBCR[0] = kbc;
	kBCL[1] = kbc; kBCR[1] = kbc;
	kBCL[3] = kbc; kBCR[3] = kbc;
	kBCL[4] = kbc; kBCR[4] = kbc;
	kBCL[5] = kbc; kBCR[5] = kbc;
	kBCL[6] = kbc; kBCR[6] = kbc;
	
	// Hard-code vertical boundary conditions on W 
	kBCL[2] = R1T0; kBCR[2] = R1T0;

	// Define the Reference state
    refstate = ref;
	
	// Assign local object pointers
	bgFields = bgU;
	rawObs = obs;
	iMin = configHash->value("i_min").toFloat();
	iMax = configHash->value("i_max").toFloat();
	DI = configHash->value("i_incr").toFloat();
	iDim = (int)((iMax - iMin)/DI) + 1;
	jMin = configHash->value("j_min").toFloat();
	jMax = configHash->value("j_max").toFloat();
	DJ = configHash->value("j_incr").toFloat();
	jDim = (int)((jMax - jMin)/DJ) + 1;
	kMin = configHash->value("k_min").toFloat();
	kMax = configHash->value("k_max").toFloat();
	DK = configHash->value("k_incr").toFloat();
	kDim = (int)((kMax - kMin)/DK) + 1;

	DIrecip = 1./DI;
	DJrecip = 1./DJ;
	DKrecip = 1./DK;
	
	// Allocate memory for the final gridded analysis
	int nodes = iDim*jDim*kDim;
	finalAnalysis = new real[nodes*45];
	
	// Adjust the internal, variable domain to exclude boundaries in R0, R2, and R3 cases
	adjustInternalDomain(1);
	
	// Redefine nodes with new domain
	nodes = iDim*jDim*kDim;
	
	// Set up the initial recursive filter
	iFilter = new RecursiveFilter(4);
	jFilter = new RecursiveFilter(4);
	kFilter = new RecursiveFilter(4);
	
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
	obsVector = new real[mObs*(7+varDim*derivDim)];
	HCq = new real[mObs+nodes];
	innovation = new real[mObs+nodes];	
	bgState = new real[nState];
	bgStdDev = new real[nState];
	stateA = new real[nState];
	stateB = new real[nState];
	stateC = new real[nState];
	iL = new real[varDim*iDim*4];
	jL = new real[varDim*jDim*4];
	kL = new real[varDim*kDim*4];

	// Precalculate the basis functions for lookup table option
	basis0 = new real[2000000];
	basis1 = new real[2000000];
	fillBasisLookup();
	
}	

void CostFunction3D::initState(const int iteration)
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
		stateA[n] = 0.0;
		stateB[n] = 0.0;
		stateC[n] = 0.0;
	}
	
	// Initialize background errors and filter scales
	bgError[0] = configHash->value("bg_rhou_error").toFloat();
	bgError[1] = configHash->value("bg_rhov_error").toFloat();
	bgError[2] = configHash->value("bg_rhow_error").toFloat();
	bgError[3] = configHash->value("bg_tempk_error").toFloat();
	bgError[4] = configHash->value("bg_qv_error").toFloat();	
	bgError[5] = configHash->value("bg_rhoa_error").toFloat();
	bgError[6] = configHash->value("bg_qr_error").toFloat();	
	
	// Set up the recursive filter
	iFilterScale = configHash->value("i_filter_length").toFloat();
	jFilterScale = configHash->value("j_filter_length").toFloat();
	kFilterScale = configHash->value("k_filter_length").toFloat();		
	iFilter->setFilterLengthScale(iFilterScale);
	jFilter->setFilterLengthScale(jFilterScale);
	kFilter->setFilterLengthScale(kFilterScale);
	
	// Set up the spline matrices
	setupSplines();

	// Flag whether or not to print the subgrid information
	outputMish = configHash->value("output_mish").toInt();
	
	// Mass continuity weight
	mcWeight = configHash->value("mc_weight").toFloat();
	cout << "Mass continuity weight set to " << mcWeight << endl;
	
	if (iteration == 1) {
		cout << "Initializing background..." << endl;
		// Set up the background state
		for (int n = 0; n < nState; n++) {
			bgState[n] = 0.0;
			bgStdDev[n] = 0.0;
		}
		
		// SB Transform on the original bg fields
		SBtransform(bgFields, stateB);
		
		// SA transform = bg B's -> bg A's
		SAtransform(stateB, bgState);
	}
	
	// Using a constant bg error variance for now, but this could be variable across the nodes
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					int bIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var;
					bgStdDev[bIndex] = bgError[var];
				}
			}
		}
	}
	
	// Compute and display the variable BG errors and RMS of values
	for (int var = 0; var < varDim; var++) {
		real varScale = 0;
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
			real errPct = 100*bgError[var]/varScale;
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
	outputAnalysis("background", bgState);

	cout << "Beginning analysis...\n";
		
	// HTd
	calcHTranspose(innovation, stateC);
	
	SCtranspose(stateC, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtransform(stateA, CTHTd);
				
}	

real CostFunction3D::funcValue(real* state)
{

	real qIP, obIP, mcIP;
	qIP = 0.;
	obIP = 0.;
	mcIP = 0.;
	
	updateHCq(state);

	// Compute inner product of state vector
	#pragma omp parallel for reduction(+:qIP)
	for (int n = 0; n < nState; n++) {
		qIP += state[n]*state[n];
	}
		
	// Subtract d from HCq to yield mObs length vector and compute inner product
	#pragma omp parallel for reduction(+:obIP)
	for (int m = 0; m < mObs; m++) {
        int obIndex = m*(7+varDim*derivDim) + 1; 
		obIP += (HCq[m]-innovation[m])*(obsVector[obIndex])*(HCq[m]-innovation[m]);
	}
		
	real J = 0.5*(qIP + obIP);
	return J;
	
}

void CostFunction3D::funcGradient(real* state, real* gradient)
{
	
	updateHCq(state);
		
	// HTHCq
	calcHTranspose(HCq, stateC);
	
	SCtranspose(stateC, stateA);
	
	SAtransform(stateA, stateB);
		
	for (int n = 0; n < nState; n++) {
		gradient[n] = state[n] + stateB[n] - CTHTd[n];
	}
	
	
}

void CostFunction3D::updateHCq(real* state)
{
	
	// S (SA transform) yield A's, put it in the tempState
	SAtransform(state, stateA);
	
	SCtransform(stateA, stateC);
	
	// H
	#pragma omp parallel for
	for (int m = 0; m < mObs; m++) {
		int mi = m*(7+varDim*derivDim);
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
        for (int var = 0; var < varDim; var++) {
            for (int d = 0; d < derivDim; d++) {
                int wgt_index = mi + (7*(d+1)) + var;
                if (!obsVector[wgt_index]) continue;
                for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                    int iNode = iiNode;
                    if ((iBCL[var] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                    if ((iBCR[var] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                    if ((iNode < 0) or (iNode >= iDim)) continue;
                    ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, derivative[d][0], iBCL[var], iBCR[var]);
                    
                    for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                        int jNode = jjNode;
                        if ((jBCL[var] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                        if ((jBCR[var] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                        if ((jNode < 0) or (jNode >= jDim)) continue;
                        jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, derivative[d][1], jBCL[var], jBCR[var]);
                        
                        for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
                            kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, derivative[d][2], kBCL[var], kBCR[var]);
                            int cIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
                            real basis = ibasis * jbasis * kbasis;
                            tempsum += stateC[cIndex + var] * basis * obsVector[wgt_index];
                        }
                    }
                }
            }
		}
		HCq[m] = tempsum;
	}

}

void CostFunction3D::updateBG()
{

	// S (SA transform) yield A's
	SAtransform(currState, stateA);
	
	SCtransform(stateA, stateC);
	
	outputAnalysis("increment", stateC);	
	
	// In BG update we are directly summing C + A
	ofstream cstream("samurai_Coefficients.out");
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
		
	outputAnalysis("analysis", bgState);
	
}

void CostFunction3D::calcInnovation()
{
	// Initialize and fill the innovation vector
	cout << "Initializing innovation vector..." << endl;
	for (int m = 0; m < mObs; m++) {
		HCq[m] = 0.0;
		innovation[m] = obsVector[m*(7+varDim*derivDim)];
	}
	
	real innovationRMS = 0.;
	#pragma omp parallel for reduction(+:innovationRMS)
	for (int m = 0; m < mObs; m++) {
		int mi = m*(7+varDim*derivDim);
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
        for (int var = 0; var < varDim; var++) {
            for (int d = 0; d < derivDim; d++) {
                int wgt_index = mi + (7*(d+1)) + var;
                if (!obsVector[wgt_index]) continue;
                for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
                    kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, derivative[d][2], kBCL[var], kBCR[var]);
                    for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                        int jNode = jjNode;
                        if ((jBCL[var] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                        if ((jBCR[var] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                        if ((jNode < 0) or (jNode >= jDim)) continue;
                        jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, derivative[d][1], jBCL[var], jBCR[var]);
                        
                        for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                            int iNode = iiNode;
                            if ((iBCL[var] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                            if ((iBCR[var] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                            if ((iNode < 0) or (iNode >= iDim)) continue;
                            
                            int stateIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
                            ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, derivative[d][0], iBCL[var], iBCR[var]);
                            tempsum += bgState[stateIndex + var] * ibasis * jbasis * kbasis * obsVector[wgt_index];
                        }
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
	for (int n = 0; n < nState; n++) {
		Astate[n] = 0.;
	}
	
	// Calculate H Transpose	
	//#pragma omp parallel for
	for (int m = 0; m < mObs; m++) {
		// Sum over obs this time
		// Multiply state by H weights
		int mi = m*(7+varDim*derivDim);
		real i = obsVector[mi+2];
		int ii = (int)((i - iMin)*DIrecip);
		
		real j = obsVector[mi+3];
		int jj = (int)((j - jMin)*DJrecip);
		
		real k = obsVector[mi+4];
		int kk = (int)((k - kMin)*DKrecip);
        for (int var = 0; var < varDim; var++) {
            for (int d = 0; d < derivDim; d++) {
                int wgt_index = mi + (7*(d+1)) + var;
                if (!obsVector[wgt_index]) continue;
                for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); kNode++) {
                    real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, derivative[d][2], kBCL[var], kBCR[var]);
                    for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                        int jNode = jjNode;
                        if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                        if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                        if ((jNode < 0) or (jNode >= jDim)) continue;
                        real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, derivative[d][1], jBCL[var], jBCR[var]);
                        
                        for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                            int iNode = iiNode;
                            if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                            if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                            if ((iNode < 0) or (iNode >= iDim)) continue;
                            
                            int aIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;			                            
                            real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, derivative[d][0], iBCL[var], iBCR[var]);
                            real invError = obsVector[mi+1];
                            real qbasise = yhat[m] * ibasis * jbasis * kbasis * invError;
                            //#pragma omp atomic
                            Astate[aIndex + var] += qbasise * obsVector[wgt_index];
                        }
                    }
                }
            }
        }
	}	
	
}

bool CostFunction3D::SAtransform(const real* Bstate, real* Astate)
{

	//#pragma omp parallel for
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

void CostFunction3D::SBtransform(const real* Ustate, real* Bstate)
{
	// Clear the Bstate
	for (int n = 0; n < nState; n++) {
		Bstate[n] = 0.;
	}
	real gausspoint = 0.5*sqrt(1./3.);
	
//#pragma omp parallel for
    for (int var = 0; var < varDim; var++) {
        for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
            for (int imu = -1; imu <= 1; imu += 2) {
                real i = iMin + DI * (iIndex + (gausspoint * imu + 0.5));
                int ii = (int)((i - iMin)*DIrecip);
                for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                    int iNode = iiNode;
                    if ((iBCL[var] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                    if ((iBCR[var] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                    if ((iNode < 0) or (iNode >= iDim)) continue;
                    
                    real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
                    int uI = iIndex*2 + (imu+1)/2;
                    
                    for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
                        for (int jmu = -1; jmu <= 1; jmu += 2) {
                            real j = jMin + DJ * (jIndex + (gausspoint * jmu + 0.5));
                            int jj = (int)((j - jMin)*DJrecip);
                            for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {                            
                                int jNode = jjNode;
                                if ((jBCL[var] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                                if ((jBCR[var] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                                if ((jNode < 0) or (jNode >= jDim)) continue;
                                
                                real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
                                int uJ = jIndex*2 + (jmu+1)/2;
                                real ijbasis = ibasis * jbasis;
                                for (int kIndex = 0; kIndex < (kDim-1); kIndex++) {
                                    for (int kmu = -1; kmu <= 1; kmu += 2) {
                                        real k = kMin + DK * (kIndex + (gausspoint * kmu + 0.5));
                                        int kk = (int)((k - kMin)*DKrecip);
                                        for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
                                            if ((kNode < 0) or (kNode >= kDim)) continue;
                                            real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
                                            real ijkbasis = 0.125 * ijbasis * kbasis;
                                            int uK = kIndex*2 + (kmu+1)/2;
                                            int uIndex = varDim*(iDim-1)*2*(jDim-1)*2*uK +varDim*(iDim-1)*2*uJ +varDim*uI;
                                            int bIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
                                            
											int ui = uIndex + var;
											if (Ustate[ui] == 0) continue;
											
											int bi = bIndex + var;
											//#pragma omp atomic
											Bstate[bi] += Ustate[ui] * ijkbasis; 
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


void CostFunction3D::SBtranspose(const real* Bstate, real* Ustate)
{
		
	// Clear the Ustate
	for (int n = 0; n < nState; n++) {
		Ustate[n] = 0;
	}
	real gausspoint = 0.5*sqrt(1./3.);
	
	//#pragma omp parallel for
    for (int var = 0; var < varDim; var++) {
        for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
            for (int imu = -1; imu <= 1; imu += 2) {
                real i = iMin + DI * (iIndex + (gausspoint * imu + 0.5));
                int ii = (int)((i - iMin)*DIrecip);
                for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                    int iNode = iiNode;
                    if ((iBCL[var] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                    if ((iBCR[var] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                    if ((iNode < 0) or (iNode >= iDim)) continue;
                    
                    real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[0], iBCR[0]);
                    int uI = iIndex*2 + (imu+1)/2;
                    
                    for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
                        for (int jmu = -1; jmu <= 1; jmu += 2) {
                            real j = jMin + DJ * (jIndex + (gausspoint * jmu + 0.5));
                            int jj = (int)((j - jMin)*DJrecip);
                            for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                                int jNode = jjNode;
                                if ((jBCL[var] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                                if ((jBCR[var] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                                if ((jNode < 0) or (jNode >= jDim)) continue;
                                
                                real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[0], jBCR[0]);
                                int uJ = jIndex*2 + (jmu+1)/2;
                                real ijbasis = ibasis * jbasis;
                                for (int kIndex = 0; kIndex < (kDim-1); kIndex++) {
                                    for (int kmu = -1; kmu <= 1; kmu += 2) {
                                        real k = kMin + DK * (kIndex + (gausspoint * kmu + 0.5));
                                        int kk = (int)((k - kMin)*DKrecip);
                                        for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
                                            if ((kNode < 0) or (kNode >= kDim)) continue;
                                            real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
                                            real ijkbasis = 0.125 * ijbasis * kbasis;
                                            int bIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
                                            int uK = kIndex*2 + (kmu+1)/2;
                                            int uIndex = varDim*(iDim-1)*2*(jDim-1)*2*uK +varDim*(iDim-1)*2*uJ +varDim*uI;										
                                            
											int bi = bIndex + var;
											if (Bstate[bi] == 0) continue;
											
											int ui = uIndex + var;
											//#pragma omp atomic
											Ustate[ui] += Bstate[bi] * ijkbasis;
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

void CostFunction3D::SCtransform(const real* Astate, real* Cstate)
{
	// Disable recursive filter if less than 1
	if ((iFilterScale < 0) and (jFilterScale < 0) and (kFilterScale < 0)) {
		#pragma omp parallel for
		for (int n = 0; n < nState; n++) {
			Cstate[n]= Astate[n] * bgStdDev[n];
		}
	} else {
		// Isotropic Recursive filter, no anisotropic "triad" working yet
              //#pragma omp parallel for
		for (int var = 0; var < varDim; var++) {
			
			//FK
			// These are local for parallelization
			real* iTemp = new real[iDim];
			real* jTemp = new real[jDim];
			real* kTemp = new real[kDim];
			
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					for (int kIndex = 0; kIndex < kDim; kIndex++) {
						kTemp[kIndex] = Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
					}
					if (kFilterScale > 0) kFilter->filterArray(kTemp, kDim);
					for (int kIndex = 0; kIndex < kDim; kIndex++) {
						Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = kTemp[kIndex];
					}
				}
			}
			
			//FJ
			for (int kIndex = 0; kIndex < kDim; kIndex++) {
				for (int iIndex = 0; iIndex < iDim; iIndex++) {
					for (int jIndex = 0; jIndex < jDim; jIndex++) {
						jTemp[jIndex] = Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
					}
					if (jFilterScale > 0) jFilter->filterArray(jTemp, jDim);
					for (int jIndex = 0; jIndex < jDim; jIndex++) {
						Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
					}
				}
			}
			//FI
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					for (int iIndex = 0; iIndex < iDim; iIndex++) {
						iTemp[iIndex] = Cstate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
					}
					if (iFilterScale > 0) iFilter->filterArray(iTemp, iDim);
					for (int iIndex = 0; iIndex < iDim; iIndex++) {
						// D
						int cIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var;
						Cstate[cIndex] = iTemp[iIndex] * bgStdDev[cIndex]; 
					}
				}
			}
			delete[] iTemp;
			delete[] jTemp;
			delete[] kTemp;
		} 
	}
}

void CostFunction3D::SCtranspose(const real* Cstate, real* Astate)
{
	if ((iFilterScale < 0) and (jFilterScale < 0) and (kFilterScale < 0)) {
		#pragma omp parallel for
		for (int n = 0; n < nState; n++) {
			Astate[n]= Cstate[n] * bgStdDev[n];
		}
	} else {
		// Isotropic Recursive filter, no anisotropic "triad" working yet 
              //_#pragma omp parallel for
		for (int var = 0; var < varDim; var++) {
			
			// These are local for parallelization
			real* iTemp = new real[iDim];
			real* jTemp = new real[jDim];
			real* kTemp = new real[kDim];
			
			//FI & D
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				for (int kIndex = 0; kIndex < kDim; kIndex++) {
					for (int iIndex = 0; iIndex < iDim; iIndex++) {
						int cIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var;
						iTemp[iIndex] = Cstate[cIndex] * bgStdDev[cIndex];
					}
					if (iFilterScale > 0) iFilter->filterArray(iTemp, iDim);
					for (int iIndex = 0; iIndex < iDim; iIndex++) {
						Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex]; 
					}
				}
			}
			//FJ
			for (int kIndex = 0; kIndex < kDim; kIndex++) {
				for (int iIndex = 0; iIndex < iDim; iIndex++) {
					for (int jIndex = 0; jIndex < jDim; jIndex++) {
						jTemp[jIndex] = Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
					}
					if (jFilterScale > 0) jFilter->filterArray(jTemp, jDim);
					for (int jIndex = 0; jIndex < jDim; jIndex++) {
						Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
					}
				}
			}
			//FK
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					for (int kIndex = 0; kIndex < kDim; kIndex++) {
						kTemp[kIndex] = Astate[varDim*iDim*jDim*kIndex + varDim*iDim*jIndex + varDim*iIndex + var];
					}
					if (kFilterScale > 0) kFilter->filterArray(kTemp, kDim);
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
}

bool CostFunction3D::setupSplines()
{
	
	// Do the spline via a Cholesky decomposition
	// and manipulate the DC Filter
	real Pi = acos(-1.);
	real** P = new real*[iDim];
	real* p = new real[iDim];
	for (int i = 0; i < iDim; i++) {
		P[i] = new real[iDim];
		p[i] = 0.;
	}
		
	for (int i = 0; i < varDim*iDim*4; i++) {
		iL[i] = 0;
	}
	
	real cutoff_wl = configHash->value("i_spline_cutoff").toFloat();
	cout << "i Spline cutoff set to " << cutoff_wl << endl;
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
				for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                    int iNode = iiNode;
                    if ((iBCL[var] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                    if ((iBCR[var] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
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
				real sum=P[i][j];
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
	for (int j = 0; j < jDim; j++) {
		P[j] = new real[jDim];
		p[j] = 0.;
	}
	
	for (int j = 0; j < varDim*jDim*4; j++) {
		jL[j] = 0;
	}
	
	cutoff_wl = configHash->value("j_spline_cutoff").toFloat();
	cout << "j Spline cutoff set to " << cutoff_wl << endl;
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
				for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                    int jNode = jjNode;
                    if ((jBCL[var] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                    if ((jBCR[var] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
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
				real sum=P[i][j];
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
	for (int k = 0; k < kDim; k++) {
		P[k] = new real[kDim];
		p[k] = 0.;
	}
	
	for (int k = 0; k < varDim*kDim*4; k++) {
		kL[k] = 0;
	}
	
	cutoff_wl = configHash->value("k_spline_cutoff").toFloat();
	cout << "k Spline cutoff set to " << cutoff_wl << endl;
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
				for (int kNode = (kk-1); kNode <= (kk+2); ++kNode) {
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
				real sum=P[i][j];
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

	return true;
	
}


void CostFunction3D::obAdjustments() {
	
	// Load the obs locally and weight the nonlinear observation operators by interpolated bg fields
	for (int m = 0; m < mObs; m++) {
		int mi = m*(7+varDim*derivDim);
		for (int ob = 0; ob < (7+varDim*derivDim); ob++) {
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
		
		for (int iiNode = ii-1; iiNode <= ii+2; ++iiNode) {
            int iNode = iiNode;
            if ((iBCL[4] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
            if ((iBCL[4] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
            if ((iNode < 0) or (iNode >= iDim)) continue;

			for (int jjNode = jj-1; jjNode <= jj+2; ++jjNode) {
                int jNode = jjNode;
                if ((jBCL[4] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                if ((jBCL[4] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                if ((jNode < 0) or (jNode >= jDim)) continue;

				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
                    if ((kNode < 0) or (kNode >= kDim)) continue;

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
		real rhoBar = refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm);
		real qBar = refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm);
		real qv = refstate->bhypInvTransform(qBar + qvprime);
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
			or (type == MetObs::ascat) or (type == MetObs::AMV)
			or (type == MetObs::lidar)) { 
			obsVector[mi] *= rhoBG;
		}
	}
}	

void CostFunction3D::fillBasisLookup()
{

	real ONESIXTH = 1./6.;
	for (int i=0; i < 2000000; i++) {
		real z = 2.0 - real(i)/1000000.;
		real b = (z*z*z) * ONESIXTH;
		z -= 1.0;
		if (z > 0)
			b -= (z*z*z) * 4 * ONESIXTH;
		basis0[i] = b;
		
		z = 2.0 - real(i)/1000000.;
		b = (z*z) * ONESIXTH;
		z -= 1.0;
		if (z > 0)
			b -= (z*z) * 4 * ONESIXTH;
		basis1[i] = b;
	}
	
}

real CostFunction3D::Basis(const int& m, const real& x, const int& M,const real& xmin, 
							const real& DX, const real& DXrecip, const int& derivative,
							const int& BL, const int& BR, const real& lambda)
{
	real b = 0;
	real xm = xmin + (m * DX);
	real delta = (x - xm) * DXrecip;
	real z = fabs(delta);
	// real ONESIXTH = 1./6.; real FOURSIXTH = 4./6.;
	if (z < 2.0) {
		real zi = z*1000000.;
		int z1 = int(zi);
		switch (derivative) {
			case 0:
				// Cheapest approximation
				b = basis0[z1];
				
				// Slightly more expensive
				//b = basis0[z1] + (basis0[z1+1]-basis0[z1])*(zi - z1);

				// Unapproximated
				/* z = 2.0 - z;
				b = (z*z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					b -= (z*z*z) * FOURSIXTH; */
				
				break;
			case 1:
				b = basis1[z1];
				
				//b = basis1[z1] + (basis1[z1+1]-basis1[z1])*(zi - z1);
				
				/* z = 2.0 - z;
				b = (z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					b -= (z*z) * FOURSIXTH; */

				b *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
				break;
			case 2:
				z = 2.0 - z;
				b = z;
				z -= 1.0;
				if (z > 0)
					b -= z * 4;
				b *= DXrecip * DXrecip;
				break;
			case 3:
				if (z > 1.0) {
					b = 1;
				} else if (z < 1.0) {
					b = -3.;
				}
				b *= ((delta > 0) ? -1.0 : 1.0) * DXrecip * DXrecip * DXrecip;
				break;
		}
	}
	if ((m > 1) and (m < M-1)) return b;
	
	// Add the boundary conditions if we get this far
	real bc = BasisBC(b, m, x, M, xmin, DX, DXrecip, derivative, BL, BR, lambda);
	return bc;
}

real CostFunction3D::BasisBC(real b, const int& m, const real& x, const int& M, const real& xmin, 
							  const real& DX, const real& DXrecip, const int& derivative,
							  const int& BL, const int& BR, const real& lambda)
{
	
	real bmod = 0;
	int node = -2;
	real coeffmod = 0.;
	// real ONESIXTH = 1./6.; real FOURSIXTH = 4./6.;
	if (m == 0) {
		// Left BC
		switch (BL) {
			case -1:
				// No boundary condition, but buffered so use R1T2
				node = -1;
				coeffmod = 2.;
				break;
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
				// For R2 BCs, the 0 node is recast as 1, with BCs applied to -1 and -2 
				node = -2;
				coeffmod = 1.;
				break;
			case 5:
				node = -2;
				coeffmod = -1.;
				break;
			case 6:
				return b;	
			case 7:
				node = M+1;
				coeffmod = 1.;
				break;
		}
	} else if (m == 1) {
		// Left BC
		switch (BL) {
			case -1:
				// No boundary condition, but buffered so use R1T0
				node = -1;
				coeffmod = -1.;
				break;			
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
				return b;
			case 5:
				return b;
			case 6:
				return b;
			case 7:
				node = M+2;
				coeffmod = 1.;
				break;
		}
	} else if (m == M) {
		// Right BC
		switch (BR) {
			case -1:
				// No boundary condition, but buffered so use R1T0
				node = M+1;
				coeffmod = 2.;
				break;					
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
				node = M+2;
				coeffmod = 1.;
				break;					
			case 5:
				node = M+2;
				coeffmod = -1.;
				break;					
			case 6:
				return b;
			case 7:
                node = -1;
                coeffmod = 1.;
                break;
		} 
	} else if (m == (M-1)) {
		// Right BC
		switch (BR) {
			case -1:
				// No boundary condition, but buffered so use R1T0
				node = M+1;
				coeffmod = -1.;
				break;
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
				return b;
			case 5:
				return b;
			case 6:
				return b;
			case 7:
                return b;
		}
	}	
	
	real xm = xmin + (node * DX);
	real delta = (x - xm) * DXrecip;
	real z = fabs(delta);
	if (z < 2.0) {
		real zi = z*1000000.;
		int z1 = int(zi);
		switch (derivative) {
			case 0:
				// Cheapest approximation
				bmod = basis0[z1];
				
				// Slightly more expensive
				//bmod = basis0[z1] + (basis0[z1+1]-basis0[z1])*(zi - z1);
				
				// Unapproximated
				/* z = 2.0 - z;
				 bmod = (z*z*z) * ONESIXTH;
				 z -= 1.0;
				 if (z > 0)
				 bmod -= (z*z*z) * FOURSIXTH; */
				
				break;
			case 1:
				bmod = basis1[z1];
				
				//bmod = basis1[z1] + (basis1[z1+1]-basis1[z1])*(zi - z1);
				
				/* z = 2.0 - z;
				 bmod = (z*z) * ONESIXTH;
				 z -= 1.0;
				 if (z > 0)
				 bmod -= (z*z) * FOURSIXTH; */
				
				bmod *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
				break;
			case 2:
				z = 2.0 - z;
				bmod = z;
				z -= 1.0;
				if (z > 0)
					bmod -= z * 4;
				bmod *= DXrecip * DXrecip;
				break;
			case 3:
				if (z > 1.0) {
					bmod = 1;
				} else if (z < 1.0) {
					bmod = -3.;
				}
				bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip * DXrecip * DXrecip;
				break;
		}
	}	
	b += coeffmod * bmod;
	
	// R2 needs one more addition
	if ((BL == 4) and (m == 0)) {
		node = -1;
		coeffmod = -0.5;
	} else if ((BR == 4) and (m == M)) {
		node = M+1;
		coeffmod = -0.5;
	} else {
		return b;
	}
	
	xm = xmin + (node * DX);
	delta = (x - xm) * DXrecip;
	z = fabs(delta);
	if (z < 2.0) {
		real zi = z*1000000.;
		int z1 = int(zi);
		switch (derivative) {
			case 0:
				bmod = basis0[z1];
				//bmod = basis0[z1] + (basis0[z1+1]-basis0[z1])*(zi - z1);	
				break;
			case 1:
				bmod = basis1[z1];
				//bmod = basis1[z1] + (basis1[z1+1]-basis1[z1])*(zi - z1);
				bmod *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
				break;
			case 2:
				z = 2.0 - z;
				bmod = z;
				z -= 1.0;
				if (z > 0)
					bmod -= z * 4;
				bmod *= DXrecip * DXrecip;
				break;
			case 3:
				if (z > 1.0) {
					bmod = 1;
				} else if (z < 1.0) {
					bmod = -3.;
				}
				bmod *= ((delta > 0) ? -1.0 : 1.0) * DXrecip * DXrecip * DXrecip;
				break;
		}
	}
	b += coeffmod * bmod;
	
	return b;
	
}

void CostFunction3D::adjustInternalDomain(int increment)
{
	if (configHash->value("i_bc") == "R0") {
		// Increase the "internal" size of the grid for the R0 condition
		iMin -= DI*increment;
		iMax += DI*increment;
		iDim += 2*increment;
	} else if ((configHash->value("i_bc") == "R2T10") or
			   (configHash->value("i_bc") == "R2T10")) {
		// Decrease the "internal" size of the grid for the R2 condition
		iMin += DI*increment;
		iMax -= DI*increment;
		iDim -= 2*increment;
	} else if (configHash->value("i_bc") == "R3") {
		// Decrease the "internal" size of the grid for the R3 conditions
		iMin += DI*2*increment;
		iMax -= DI*2*increment;
		iDim -= 4*increment;
	}
	
	if (configHash->value("j_bc") == "R0") {
		jMin -= DJ*increment;
		jMax += DJ*increment;
		jDim += 2*increment;
	} else if ((configHash->value("j_bc") == "R2T10") or
			   (configHash->value("j_bc") == "R2T20")) {
		jMin += DJ*increment;
		jMax -= DJ*increment;
		jDim -= 2*increment;
	} else if (configHash->value("j_bc") == "R3") {
		jMin += DJ*2*increment;
		jMax -= DJ*2*increment;
		jDim -= 4*increment;
	}
	
	if (configHash->value("k_bc") == "R0") {
		kMin -= DK*increment;
		kMax += DK*increment;
		kDim += 2*increment;
	} else if ((configHash->value("k_bc") == "R2T20") or
			   (configHash->value("k_bc") == "R2T20")) {
		kMin += DK*increment;
		kMax -= DK*increment;
		kDim -= 2*increment;
	} else if (configHash->value("k_bc") == "R3") {
		kMin += DK*2*increment;
		kMax -= DK*2*increment;
		kDim -= 4*increment;
	}
}	

