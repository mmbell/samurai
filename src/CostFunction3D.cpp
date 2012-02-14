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
	delete[] fieldNodes;
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

void CostFunction3D::initialize(const QHash<QString, QString>* config, real* bgU, real* obs)
{

	// Initialize number of variables
	varDim = 7;
	configHash = config;

	// Horizontal boundary conditions
	int ibc = bcHash.value(configHash->value("ibc"));	
	iBCL[0] = ibc; iBCR[0] = ibc; 
	iBCL[1] = ibc; iBCR[1] = ibc; 
	iBCL[2] = ibc; iBCR[2] = ibc; 
	iBCL[3] = ibc; iBCR[3] = ibc; 
	iBCL[4] = ibc; iBCR[4] = ibc; 
	iBCL[5] = ibc; iBCR[5] = ibc; 
	iBCL[6] = ibc; iBCR[6] = ibc; 
	
	int jbc = bcHash.value(configHash->value("jbc"));	
	jBCL[0] = jbc; jBCR[0] = jbc; 
	jBCL[1] = jbc; jBCR[1] = jbc;
	jBCL[2] = jbc; jBCR[2] = jbc;
	jBCL[3] = jbc; jBCR[3] = jbc;
	jBCL[4] = jbc; jBCR[4] = jbc;
	jBCL[5] = jbc; jBCR[5] = jbc;
	jBCL[6] = jbc; jBCR[6] = jbc;
	
	int kbc = bcHash.value(configHash->value("kbc"));
	kBCL[0] = kbc; kBCR[0] = kbc;
	kBCL[1] = kbc; kBCR[1] = kbc;
	kBCL[3] = kbc; kBCR[3] = kbc;
	kBCL[4] = kbc; kBCR[4] = kbc;
	kBCL[5] = kbc; kBCR[5] = kbc;
	kBCL[6] = kbc; kBCR[6] = kbc;
	
	// Hard-code vertical boundary conditions on W 
	kBCL[2] = R1T0; kBCR[2] = R1T0;

	// Define the Reference state
	if (configHash->value("refstate") == "dunion_mt") {
		referenceState = dunion_mt;
	}
	
	// Assign local object pointers
	bgFields = bgU;
	rawObs = obs;
	iMin = configHash->value("imin").toFloat();
	iMax = configHash->value("imax").toFloat();
	DI = configHash->value("iincr").toFloat();
	iDim = (int)((iMax - iMin)/DI) + 1;
	jMin = configHash->value("jmin").toFloat();
	jMax = configHash->value("jmax").toFloat();
	DJ = configHash->value("jincr").toFloat();
	jDim = (int)((jMax - jMin)/DJ) + 1;
	kMin = configHash->value("kmin").toFloat();
	kMax = configHash->value("kmax").toFloat();
	DK = configHash->value("kincr").toFloat();
	kDim = (int)((kMax - kMin)/DK) + 1;

	DIrecip = 1./DI;
	DJrecip = 1./DJ;
	DKrecip = 1./DK;

	// Increase the "internal" size of the grid for the zero BC condition
	if (configHash->value("ibc") == "R0") {
		iMin -= DI;
		iMax += DI;
		iDim += 2;
	}
	if (configHash->value("jbc") == "R0") {
		jMin -= DJ;
		jMax += DJ;
		jDim += 2;
	}
	if (configHash->value("kbc") == "R0") {
		kMin -= DK;
		kMax += DK;
		kDim += 2;
		// Keep the vertical ones the same, since they are outside the domain anyway?
		//kBCL[2] = R0; kBCR[2] = R0;
	}
	
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
	obsVector = new real[mObs*14];
	int nodes = iDim*jDim*kDim;
	HCq = new real[mObs+nodes];
	innovation = new real[mObs+nodes];	
	fieldNodes = new real[nodes*33];	
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
	bgError[0] = configHash->value("uerror").toFloat();
	bgError[1] = configHash->value("verror").toFloat();
	bgError[2] = configHash->value("werror").toFloat();
	bgError[3] = configHash->value("terror").toFloat();
	bgError[4] = configHash->value("qverror").toFloat();	
	bgError[5] = configHash->value("rhoerror").toFloat();
	bgError[6] = configHash->value("qrerror").toFloat();	
	
	// Set up the recursive filter
	iFilterScale = configHash->value("ifilter").toFloat();
	jFilterScale = configHash->value("jfilter").toFloat();
	kFilterScale = configHash->value("kfilter").toFloat();		
	iFilter->setFilterLengthScale(iFilterScale);
	jFilter->setFilterLengthScale(jFilterScale);
	kFilter->setFilterLengthScale(kFilterScale);
	
	// Set up the spline matrices
	setupSplines();

	// Flag whether or not to print the subgrid information
	outputMish = configHash->value("output_mish").toInt();
	
	// Mass continuity weight
	mcWeight = configHash->value("mcweight").toFloat();
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
		obIP += (HCq[m]-innovation[m])*(obsVector[m*14+1])*(HCq[m]-innovation[m]);
	}
	
	// Mass continuity on the nodes
	#pragma omp parallel for reduction(+:mcIP)
	for (int kIndex = 0; kIndex < kDim; kIndex++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				int hIndex = mObs + iDim*jDim*kIndex + iDim*jIndex + iIndex;
				mcIP += (HCq[hIndex]-innovation[hIndex])*mcWeight*(HCq[hIndex]-innovation[hIndex]);
			}
		}
	}
	
	real J = 0.5*(qIP + obIP + mcIP);
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
		int mi = m*14;
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
		for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
            int iNode = iiNode;
            if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
            if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
            if ((iNode < 0) or (iNode >= iDim)) continue;
			ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[0], iBCR[0]);
			
			for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                int jNode = jjNode;
                if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                if ((jNode < 0) or (jNode >= jDim)) continue;
				jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[0], jBCR[0]);
				
				for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
					kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[0], kBCR[0]);
					int cIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
					real basis = ibasis * jbasis * kbasis;
					for (int var = 0; var < varDim; var++) {
						if (!obsVector[mi+7 + var]) continue;
						if ((kBCL[var] != kBCL[0]) or (kBCR[var] != kBCR[0])) {
							kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
							basis = ibasis * jbasis * kbasis;
						} 
						tempsum += stateC[cIndex + var] * basis * obsVector[mi+7 + var];
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
				for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
					for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
						for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                            int iNode = iiNode;
                            if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                            if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                            if ((iNode < 0) or (iNode >= iDim)) continue;

                            int jNode = jjNode;
                            if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                            if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                            if ((jNode < 0) or (jNode >= jDim)) continue;

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
		innovation[m] = obsVector[m*14];
	}
	
	real innovationRMS = 0.;
	real mcbgRMS = 0.;
	#pragma omp parallel for reduction(+:innovationRMS)
	for (int m = 0; m < mObs; m++) {
		int mi = m*14;
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
        for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
            for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                    int iNode = iiNode;
                    if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                    if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                    if ((iNode < 0) or (iNode >= iDim)) continue;
                    
                    int jNode = jjNode;
                    if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                    if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                    if ((jNode < 0) or (jNode >= jDim)) continue;

					int stateIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
					for (int var = 0; var < varDim; var++) {
						if (obsVector[mi+7 + var] == 0) continue;
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
						kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
						tempsum += bgState[stateIndex + var] * ibasis * jbasis * kbasis * obsVector[mi+7+var];
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
				real tempsum = 0;
				for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
                    for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                        for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                            int iNode = iiNode;
                            if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                            if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                            if ((iNode < 0) or (iNode >= iDim)) continue;
                            
                            int jNode = jjNode;
                            if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                            if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                            if ((jNode < 0) or (jNode >= jDim)) continue;

							int bgIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;

							real idbasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 1, iBCL[0], iBCR[0]);
							real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[0], jBCR[0]);
							real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[0], kBCR[0]);							
							tempsum += (bgState[bgIndex] * idbasis * jbasis * kbasis);
							
							real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[1], iBCR[1]);
							real jdbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 1, jBCL[1], jBCR[1]);
							kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[1], kBCR[1]);							
							tempsum += (bgState[bgIndex + 1] * ibasis * jdbasis * kbasis);
										
							ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[2], iBCR[2]);
							jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[2], jBCR[2]);
							real kdbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 1, kBCL[2], kBCR[2]);							
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
	for (int n = 0; n < nState; n++) {
		Astate[n] = 0.;
	}
	
	// Calculate H Transpose	
	//#pragma omp parallel for
	for (int m = 0; m < mObs; m++) {
		// Sum over obs this time
		// Multiply state by H weights
		int mi = m*14;
		real i = obsVector[mi+2];
		int ii = (int)((i - iMin)*DIrecip);
		
		real j = obsVector[mi+3];
		int jj = (int)((j - jMin)*DJrecip);
		
		real k = obsVector[mi+4];
		int kk = (int)((k - kMin)*DKrecip);
		for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); kNode++) {
            for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                    int iNode = iiNode;
                    if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                    if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                    if ((iNode < 0) or (iNode >= iDim)) continue;
                    
                    int jNode = jjNode;
                    if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                    if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                    if ((jNode < 0) or (jNode >= jDim)) continue;

					int aIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;			

					real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[0], iBCR[0]);
					real jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[0], jBCR[0]);
					real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[0], kBCR[0]);
					real invError = obsVector[mi+1];
					real qbasise = yhat[m] * ibasis * jbasis * kbasis *invError;
					for (int var = 0; var < varDim; var++) {
						if (!obsVector[mi+7 + var]) continue;
						if ((kBCL[var] != iBCL[0]) or (kBCR[var] != kBCR[0])) {
							kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
							qbasise = yhat[m] * ibasis * jbasis * kbasis *invError;
						}
						//#pragma omp atomic
						Astate[aIndex + var] += qbasise * obsVector[mi+7+var];
					 }
				}
			}
		}
	}
	
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int kIndex = 0; kIndex < kDim; kIndex++) {
				
				int aIndex = varDim*iDim*jDim*kIndex + varDim*iDim*jIndex +varDim*iIndex;	
				
				// Mass continuity transpose				
				int ii = iIndex;
				int jj = jIndex;
				int kk = kIndex;
				for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
                    for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                        for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                            int iNode = iiNode;
                            if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                            if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                            if ((iNode < 0) or (iNode >= iDim)) continue;
                            
                            int jNode = jjNode;
                            if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                            if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                            if ((jNode < 0) or (jNode >= jDim)) continue;
														
							real i = iNode*DI + iMin;
							real j = jNode*DJ + jMin;
							real k = kNode*DK + kMin;
							
							int hIndex = mObs + iDim*jDim*kNode + iDim*jNode + iNode;
							if (yhat[hIndex] == 0) continue;
							
							real idbasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 1, iBCL[0], iBCR[0]);
							real jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[0], jBCR[0]);
							real kbasis = Basis(kIndex, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[0], kBCR[0]);							
							//#pragma omp atomic
							Astate[aIndex] += mcWeight * (yhat[hIndex] * idbasis * jbasis * kbasis);
							
							real ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[1], iBCR[1]);
							real jdbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 1, jBCL[1], jBCR[1]);
							// Extra calcs unless horizontal BC are different
							//kbasis = Basis(kIndex, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[1], kBCR[1]);							
							//#pragma omp atomic
							Astate[aIndex + 1] += mcWeight * (yhat[hIndex] * ibasis * jdbasis * kbasis);
							
							//ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[2], iBCR[2]);
							//jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[2], jBCR[2]);
							real kdbasis = Basis(kIndex, k, kDim-1, kMin, DK, DKrecip, 1, kBCL[2], kBCR[2]);							
							//#pragma omp atomic
							Astate[aIndex + 2] += mcWeight * (yhat[hIndex] * ibasis * jbasis * kdbasis);
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
	for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
		for (int imu = -1; imu <= 1; imu += 2) {
			real i = iMin + DI * (iIndex + (gausspoint * imu + 0.5));
			int ii = (int)((i - iMin)*DIrecip);
			for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                int iNode = iiNode;
                if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                if ((iNode < 0) or (iNode >= iDim)) continue;
                
				real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[0], iBCR[0]);
				int uI = iIndex*2 + (imu+1)/2;
				
				for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
					for (int jmu = -1; jmu <= 1; jmu += 2) {
						real j = jMin + DJ * (jIndex + (gausspoint * jmu + 0.5));
						int jj = (int)((j - jMin)*DJrecip);
						for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {                            
                            int jNode = jjNode;
                            if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                            if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
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
										real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[0], kBCR[0]);
										real ijkbasis = 0.125 * ijbasis * kbasis;
										int uK = kIndex*2 + (kmu+1)/2;
										int uIndex = varDim*(iDim-1)*2*(jDim-1)*2*uK +varDim*(iDim-1)*2*uJ +varDim*uI;
										int bIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;

										for (int var = 0; var < varDim; var++) {
											int ui = uIndex + var;
											if (Ustate[ui] == 0) continue;
											
											int bi = bIndex + var;
											if ((kBCL[var] != kBCL[0]) or (kBCR[var] != kBCL[0])) {
												kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
												ijkbasis = 0.125 * ibasis * jbasis * kbasis;
											}
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
	for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
		for (int imu = -1; imu <= 1; imu += 2) {
			real i = iMin + DI * (iIndex + (gausspoint * imu + 0.5));
			int ii = (int)((i - iMin)*DIrecip);
			for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                int iNode = iiNode;
                if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                if ((iNode < 0) or (iNode >= iDim)) continue;
                
				real ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[0], iBCR[0]);
				int uI = iIndex*2 + (imu+1)/2;
				
				for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
					for (int jmu = -1; jmu <= 1; jmu += 2) {
						real j = jMin + DJ * (jIndex + (gausspoint * jmu + 0.5));
						int jj = (int)((j - jMin)*DJrecip);
						for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                            int jNode = jjNode;
                            if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                            if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
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
										real kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[0], kBCR[0]);
										real ijkbasis = 0.125 * ijbasis * kbasis;
										int bIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
										int uK = kIndex*2 + (kmu+1)/2;
										int uIndex = varDim*(iDim-1)*2*(jDim-1)*2*uK +varDim*(iDim-1)*2*uJ +varDim*uI;										

										for (int var = 0; var < varDim; var++) {
											int bi = bIndex + var;
											if (Bstate[bi] == 0) continue;
											
											int ui = uIndex + var;
											if ((kBCL[var] != kBCL[0]) or (kBCR[var] != kBCL[0])) {
												kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
												ijkbasis = 0.125 * ibasis * jbasis * kbasis;
											}
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
                    if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                    if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
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
                    if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                    if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
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

bool CostFunction3D::outputAnalysis(const QString& suffix, real* Astate)
{
	
	cout << "Outputting " << suffix.toStdString() << "...\n";
	// H --> to Mish for output
	QString samuraiout = "samurai_XYZ_" + suffix + ".out";
	ofstream samuraistream(samuraiout.toAscii().data());
	samuraistream << "X\tY\tZ\trhoE\tu\tv\tw\tVorticity\tDivergence\tqv'\trho'\tT'\tP'\th\t";
	samuraistream << "udx\tudy\tudz\tvdx\tvdy\tvdz\twdx\twdy\twdz\trhowdz\tMC residual\tdBZ\n";
	samuraistream.precision(10);
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		for (int ihalf = 0; ihalf <= outputMish; ihalf++) {
			for (int imu = -ihalf; imu <= ihalf; imu++) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5*ihalf));
				if (i > ((iDim-1)*DI + iMin)) continue;
				
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					for (int jhalf =0; jhalf <=outputMish; jhalf++) {
						for (int jmu = -jhalf; jmu <= jhalf; jmu++) {
							real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5*jhalf));
							if (j > ((jDim-1)*DJ + jMin)) continue;	
							
							real tpw = 0;
							
							for (int kIndex = 0; kIndex < kDim; kIndex++) {
								for (int khalf =0; khalf <=outputMish; khalf++) {
									for (int kmu = -khalf; kmu <= khalf; kmu++) {
										real k = kMin + DK * (kIndex + (0.5*sqrt(1./3.) * kmu + 0.5*khalf));
										if (k > ((kDim-1)*DK + kMin)) continue;	
										
										real heightm = 1000*k;
										real rhoBar = getReferenceVariable(rhoaref, heightm);
										real qBar = getReferenceVariable(qvbhypref, heightm);
										real tBar = getReferenceVariable(tempref, heightm);

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
										real rhovdx = 0.; real rhoudx = 0.; real rhowdx = 0.;
										real rhovdy = 0.; real rhoudy = 0.; real rhowdy = 0.;
										real rhovdz = 0.; real rhoudz = 0.; real rhowdz = 0.;
										real tprime = 0.;
										real rhoadx = 0.; real rhoady = 0.; real rhoadz = 0.;
										real qvdx = 0.; real qvdy = 0.; real qvdz = 0.;
										real qvprime = 0.;
										real rhoprime = 0.;
										real qrprime = 0.;
										for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
                                            for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                                                int iNode = iiNode;
                                                if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                                                if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                                                if ((iNode < 0) or (iNode >= iDim)) continue;
                                                
                                                for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                                                    int jNode = jjNode;
                                                    if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                                                    if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                                                    if ((jNode < 0) or (jNode >= jDim)) continue;
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
																tprime += Astate[aIndex + 3] * basis3x;
																break;
															case 4:
																qvprime += Astate[aIndex + 4] * basis3x;
																qvdx += Astate[aIndex + 4] * idbasis * jbasis * kbasis;
																qvdy += Astate[aIndex + 4] * ibasis * jdbasis * kbasis;
																qvdz += Astate[aIndex + 4] * ibasis * jbasis * kdbasis; 
																break;
															case 5:
																rhoprime += Astate[aIndex + 5] * basis3x;
																rhoadx += Astate[aIndex + 5] * idbasis * jbasis * kbasis;
																rhoady += Astate[aIndex + 5] * ibasis * jdbasis * kbasis;
																rhoadz += Astate[aIndex + 5] * ibasis * jbasis * kdbasis;
																break;
															case 6:
																qrprime += Astate[aIndex + 6] * basis3x;
																break;
														}
													}
												}
											}
										}
																				
										// Output it
										
										// Skip the border nodes if the boundary conditions are zero
										int fIndex = iDim*jDim*kDim; 
										int posIndex = iDim*jDim*kIndex + iDim*jIndex + iIndex;
										if (configHash->value("ibc") == "R0") {
											if ((iIndex == 0) or (iIndex == iDim-1) or
												(jIndex == 0) or (jIndex == jDim-1)) {
												for (int n = 0; n < 33; ++n) {
													fieldNodes[fIndex * n + posIndex] = -999.0;
												}
												continue;
											}
										}
										
										if (configHash->value("ibc") == "R0") {
											if ((kIndex == 0) or (kIndex == kDim-1)) {
												for (int n = 0; n < 33; ++n) {
													fieldNodes[fIndex * n + posIndex] = -999.0;
												}
												continue;
											}
										}
										
										real rhoa = rhoBar + rhoprime / 100;
										real qv = bhypInvTransform(qBar + qvprime);
										real qr; 
										QString gridref = configHash->value("qrvariable");
										if (gridref == "dbz") {
											qr = qrprime*10. - 35.;
											if (qr < -35.) {
												qr = -999.;
											}
										} else {
											qr = bhypInvTransform(qrprime);
										}
										real rhoq = qv * rhoa / 1000.;
										real rho = rhoa + rhoq;
										real v = rhov / rho;
										real u = rhou / rho;
										real w = rhow / rho;
										real wspd = sqrt(u*u + v*v);
										real KE = 0.5*rho*(v*v + u*u + w*w);
										real temp = tBar + tprime;
										real h = 1005.7*temp + 2.501e3*qv + 9.81*heightm;
										real rhoE = rho*h + KE;
										real airpress = temp*rhoa*287./100.;
										real tempc = temp - 273.15;
										real satvp = 6.112 * exp((17.62 * tempc)/(243.12 + tempc));
										real vp = temp*rhoq*461./100.;
										real relhum = -999.;
										if (satvp != 0)
											relhum = 100*vp/satvp;
										real press = airpress + vp;
										
										real pprime = press - getReferenceVariable(pressref, heightm)/100.;
										real hprime = h - getReferenceVariable(href, heightm);
										
										// Calculate the kinematic derivatives
										// rhoa derivatives divided by 100
										// qv derivatives multipled by 2 to account for hyperbolic transform, not exact but close enough
										real rhodx = rhoadx * (1. + qv/1000.) / 100. + rhoa * qvdx/500.;
										real rhody = rhoady * (1. + qv/1000.) / 100. + rhoa * qvdy/500.;
										real rhodz = rhoadz * (1. + qv/1000.) / 100. + rhoa * qvdz/500.;
										real rhobardz = 1000 * getReferenceVariable(rhoref, heightm, 1);
										rhodz += rhobardz;
										// Units 10-5
										real udx = 100. * (rhoudx - u*rhodx) / rho;
										real udy = 100. * (rhoudy - u*rhody) / rho;
										real udz = 100. * (rhoudz - u*rhodz) / rho;
										
										real vdx = 100. * (rhovdx - v*rhodx) / rho;
										real vdy = 100. * (rhovdy - v*rhody) / rho;
										real vdz = 100. * (rhovdz - v*rhodz) / rho;
										
										real wdx = 100. * (rhowdx - w*rhodx) / rho;
										real wdy = 100. * (rhowdy - w*rhody) / rho;
										real wdz = 100. * (rhowdz - w*rhodz) / rho;
										
										// Vorticity units are 10-5
										real vorticity = (vdx - udy);
										real divergence = (udx + vdy);
										real s1 = (udx - vdy);
										real s2 = (vdx + udy);
										real strain = sqrt(s1*s1 + s2*s2);
										real okuboweiss = vorticity*vorticity - s1*s1 -s2*s2;
										real mcresidual = rhoudx + rhovdy + rhowdz;

										QString refmask = configHash->value("mask_reflectivity");
										if (refmask != "None") {
											real refthreshold = refmask.toFloat();
											if (qr < refthreshold) {
												u = -999.;
												v = -999.;
												w = -999.;
												wspd = -999.;
												relhum = -999.;
												hprime = -999.;
												qvprime = -999.;
												rhoprime = -999.;
												tprime = -999.;
												pprime = -999.;
												vorticity = -999.;
												divergence = -999.;
												okuboweiss = -999.;
												strain = -999.;
												tpw = -999.;
												rhou = -999.;
												rhov = -999.;
												rhow = -999.;
												rho = -999.;
												press = -999.;
												temp = -999.;
												qv = -999.;
												h = -999.;
												qr = -999.;
												udx = -999.;
												vdx = -999.;
												wdx = -999.;
												udy = -999.;
												vdy = -999.;
												wdy = -999.;
												udz = -999.;
												vdz = -999.;
												wdz = -999.;
												rhoE = -999.;
											}
										}

										samuraistream << scientific << i << "\t" << j << "\t"  << k << "\t" << rhoE
										<< "\t" << u << "\t" << v << "\t" << w << "\t" << vorticity << "\t" << divergence
										<< "\t" << qvprime*2 << "\t" << rhoprime << "\t" << tprime << "\t" << pprime <<  "\t" << hprime << "\t"
										<< udx << "\t" << udy << "\t" << udz << "\t"
										<< vdx << "\t" << vdy << "\t" << vdz << "\t"
										<< wdx << "\t" << wdy << "\t" << wdz << "\t" 
										<< rhowdz * 100. << "\t" << mcresidual << "\t" << qr << "\n";
										
										// Sum up the TPW in the vertical, top level is tpw
										tpw += qv * rhoa * DK;
										
										// On the nodes	
										if (!ihalf and !jhalf and !khalf){
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
											fieldNodes[fIndex * 23 + posIndex] = qr;
											fieldNodes[fIndex * 24 + posIndex] = udx;
											fieldNodes[fIndex * 25 + posIndex] = vdx;
											fieldNodes[fIndex * 26 + posIndex] = wdx;
											fieldNodes[fIndex * 27 + posIndex] = udy;
											fieldNodes[fIndex * 28 + posIndex] = vdy;
											fieldNodes[fIndex * 29 + posIndex] = wdy;
											fieldNodes[fIndex * 30 + posIndex] = udz;
											fieldNodes[fIndex * 31 + posIndex] = vdz;
											fieldNodes[fIndex * 32 + posIndex] = wdz;
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
	// Write out to a netCDF file	
	QString cdfFileName = outFileName + ".nc";
	if (!writeNetCDF(cdfFileName))
		cout << "Error writing netcdf file " << cdfFileName.toStdString() << endl; 	
	
	// Write out to an asi file
	QString asiFileName = outFileName + ".asi";
	if (!writeAsi(asiFileName))
		cout << "Error writing asi file " << asiFileName.toStdString() << endl; 	
		
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
	*os++ = "Time";
	*os++ = "rhou";
	*os++ = "rhov";
	*os++ = "rhow";
	*os++ = "T'";
	*os++ = "qv'";
	*os++ = "rhoa'";
	*os++ = "qr";
	*os++ = "Analysis";
	*os++ = "Background";
	qcstream << endl;
	qcstream.precision(10);

	ostream_iterator<real> od(qcstream, "\t ");
	for (int m = 0; m < mObs; m++) {
		int mi = m*14;
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
		for (int kNode = max(kk-1,0); kNode <= min(kk+2,kDim-1); ++kNode) {
            for (int iiNode = (ii-1); iiNode <= (ii+2); ++iiNode) {
                int iNode = iiNode;
                if ((iBCL[0] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                if ((iBCL[0] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                if ((iNode < 0) or (iNode >= iDim)) continue;
                
                for (int jjNode = (jj-1); jjNode <= (jj+2); ++jjNode) {
                    int jNode = jjNode;
                    if ((jBCL[0] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                    if ((jBCL[0] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                    if ((jNode < 0) or (jNode >= jDim)) continue;

					int aIndex = varDim*iDim*jDim*kNode + varDim*iDim*jNode +varDim*iNode;
					for (int var = 0; var < varDim; var++) {
						ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 0, iBCL[var], iBCR[var]);
						jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 0, jBCL[var], jBCR[var]);
						kbasis = Basis(kNode, k, kDim-1, kMin, DK, DKrecip, 0, kBCL[var], kBCR[var]);
						tempsum += Astate[aIndex + var] * ibasis * jbasis * kbasis * obsVector[mi+7+var];
					}
				}
			}
		}
				
		for (int t=0; t < 6; t++) {
			*od++ = obsVector[mi+t];
		}
		int unixtime = (int)obsVector[mi+6];
		QDateTime obtime;
		obtime.setTime_t(unixtime);
		obtime.setTimeSpec(Qt::UTC);
		QString timestring = obtime.toString("hh:mm:ss.zzz");
		qcstream << timestring.toStdString() << "\t";
		
		// Multiply the weight by the ob -- Observations.in has individual weights already
		for (int t=7; t<14; t++) {
			*od++ = obsVector[mi+t] * obsVector[mi];
		}
		
		*od++ = tempsum;
		*od++ = obsVector[mi]-innovation[m];
		qcstream << endl;
		
	}
	
	return true;
	
}     

bool CostFunction3D::writeNetCDF(const QString& netcdfFileName)
{
	NcError err(NcError::verbose_nonfatal);
	int NC_ERR = 0;
	
	// Create the file.
	NcFile dataFile(netcdfFileName.toAscii(), NcFile::Replace);
	
	// Check to see if the file was created.
	if(!dataFile.is_valid())
		return NC_ERR;
	
	// Define the dimensions. NetCDF will hand back an ncDim object for
	// each.
	NcDim *lvlDim, *latDim, *lonDim, *timeDim;
	if (!(lonDim = dataFile.add_dim("longitude", iDim)))
		return NC_ERR;
	if (!(latDim = dataFile.add_dim("latitude", jDim)))
		return NC_ERR;
	if (!(lvlDim = dataFile.add_dim("altitude", kDim)))
		return NC_ERR;
	// Add an unlimited dimension...
	if (!(timeDim = dataFile.add_dim("time")))
		return NC_ERR;
	
	// Define the coordinate variables.
	NcVar *latVar, *lonVar, *lvlVar, *timeVar;
	if (!(lonVar = dataFile.add_var("longitude", ncFloat, lonDim)))
		return NC_ERR;
	if (!(latVar = dataFile.add_var("latitude", ncFloat, latDim)))
		return NC_ERR;
	if (!(lvlVar = dataFile.add_var("altitude", ncFloat, lvlDim)))
		return NC_ERR;
	if (!(timeVar = dataFile.add_var("time", ncInt, timeDim)))
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
	NcVar *rho, *press, *temp, *qv, *h, *qr;
	NcVar *dudx, *dvdx, *dwdx, *dudy, *dvdy, *dwdy, *dudz, *dvdz, *dwdz;
		
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
	if (!(qr = dataFile.add_var("QR", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dudx = dataFile.add_var("DUDX", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dvdx = dataFile.add_var("DVDX", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dwdx = dataFile.add_var("DWDX", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dudy = dataFile.add_var("DUDY", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dvdy = dataFile.add_var("DVDY", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dwdy = dataFile.add_var("DWDY", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dudz = dataFile.add_var("DUDZ", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dvdz = dataFile.add_var("DVDZ", ncFloat, timeDim, 
									lvlDim, latDim, lonDim)))
		return NC_ERR;
	if (!(dwdz = dataFile.add_var("DWDZ", ncFloat, timeDim, 
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
	if (!rho->add_att("units", "kg m-3")) 
		return NC_ERR;
	if (!press->add_att("units", "hPa")) 
		return NC_ERR;
	if (!temp->add_att("units", "K")) 
		return NC_ERR;	
	if (!qv->add_att("units", "g kg-1")) 
		return NC_ERR;
	if (!h->add_att("units", "kJ")) 
		return NC_ERR;
	if (!qr->add_att("units", "g kg-1")) 
		return NC_ERR;
	if (!dudx->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dvdx->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dwdx->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dudy->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dvdy->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dwdy->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dudz->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dvdz->add_att("units", "10-5s-1")) 
		return NC_ERR;
	if (!dwdz->add_att("units", "10-5s-1")) 
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
	if (!qr->add_att("long_name", "precipitation mixing ratio")) 
		return NC_ERR;
	if (!dudx->add_att("long_name", "wind gradient")) 
		return NC_ERR;	
	if (!dvdx->add_att("long_name", "wind gradient")) 
		return NC_ERR;	
	if (!dwdx->add_att("long_name", "wind gradient")) 
		return NC_ERR;
	if (!dudy->add_att("long_name", "wind gradient")) 
		return NC_ERR;	
	if (!dvdy->add_att("long_name", "wind gradient")) 
		return NC_ERR;	
	if (!dwdy->add_att("long_name", "wind gradient")) 
		return NC_ERR;
	if (!dudz->add_att("long_name", "wind gradient")) 
		return NC_ERR;	
	if (!dvdz->add_att("long_name", "wind gradient")) 
		return NC_ERR;	
	if (!dwdz->add_att("long_name", "wind gradient")) 
		return NC_ERR;
	
	// Define missing data
	if (!u->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!v->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!w->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!wspd->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!relhum->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!hprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!qvprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rhoprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!tprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!pprime->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!vorticity->add_att("missing_value", -999.f)) 
		return NC_ERR;
	if (!divergence->add_att("missing_value", -999.f)) 
		return NC_ERR;
	if (!okuboweiss->add_att("missing_value", -999.f)) 
		return NC_ERR;
	if (!strain->add_att("missing_value", -999.f))
		return NC_ERR;       
	if (!tpw->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rhou->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rhov->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rhow->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!rho->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!press->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!temp->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!qv->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!h->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!qr->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dudx->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dvdx->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dwdx->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dudy->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dvdy->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dwdy->add_att("missing_value", -999.f))
		return NC_ERR;
	if (!dudz->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dvdz->add_att("missing_value", -999.f))
		return NC_ERR;	
	if (!dwdz->add_att("missing_value", -999.f))
		return NC_ERR;
	
	// Define _Fill_Value for NCL users
	if (!u->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!v->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!w->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!wspd->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!relhum->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!hprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!qvprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rhoprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!tprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!pprime->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!vorticity->add_att("_FillValue", -999.f)) 
		return NC_ERR;
	if (!divergence->add_att("_FillValue", -999.f)) 
		return NC_ERR;
	if (!okuboweiss->add_att("_FillValue", -999.f)) 
		return NC_ERR;
	if (!strain->add_att("_FillValue", -999.f))
		return NC_ERR;       
	if (!tpw->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rhou->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rhov->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rhow->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!rho->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!press->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!temp->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!qv->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!h->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!qr->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dudx->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dvdx->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dwdx->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dudy->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dvdy->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dwdy->add_att("_FillValue", -999.f))
		return NC_ERR;
	if (!dudz->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dvdz->add_att("_FillValue", -999.f))
		return NC_ERR;	
	if (!dwdz->add_att("_FillValue", -999.f))
		return NC_ERR;
	
	// Write the coordinate variable data to the file.
	real *lons = new real[iDim];
	real *lats = new real[jDim];
	real *levs = new real[kDim];
	int time[2];
	
	// Reference time and position from center file 
	time[0] = configHash->value("reftime").toInt();
	real latReference = configHash->value("reflat").toFloat();	
	real lonReference = configHash->value("reflon").toFloat();
	real refX, refY;
	
	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM;
	tm.Forward(lonReference, latReference, lonReference, refX, refY);
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		real i = (iMin + DI * iIndex)*1000;
		real j = (jMin + DJ * (jDim/2))*1000;
		real latnull = 0;
		tm.Reverse(lonReference,refX + i, refY + j, latnull, lons[iIndex]);
	}
	
	for (int jIndex = 0; jIndex < jDim; jIndex++) {
		real i = (iMin + DI * (iDim/2))*1000;
		real j = (jMin + DJ * jIndex)*1000;
		real lonnull = 0;
		tm.Reverse(lonReference,refX + i, refY + j, lats[jIndex], lonnull);		
	}
	
	if (!lonVar->put(lons, iDim))
		return NC_ERR;       

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
		if (!qr->put_rec(&fieldNodes[iDim*jDim*kDim*23], rec)) 
			return NC_ERR;
		if (!dudx->put_rec(&fieldNodes[iDim*jDim*kDim*24], rec)) 
			return NC_ERR;
		if (!dvdx->put_rec(&fieldNodes[iDim*jDim*kDim*25], rec)) 
			return NC_ERR;
		if (!dwdx->put_rec(&fieldNodes[iDim*jDim*kDim*26], rec)) 
			return NC_ERR;
		if (!dudy->put_rec(&fieldNodes[iDim*jDim*kDim*27], rec)) 
			return NC_ERR;
		if (!dvdy->put_rec(&fieldNodes[iDim*jDim*kDim*28], rec)) 
			return NC_ERR;
		if (!dwdy->put_rec(&fieldNodes[iDim*jDim*kDim*29], rec)) 
			return NC_ERR;
		if (!dudz->put_rec(&fieldNodes[iDim*jDim*kDim*30], rec)) 
			return NC_ERR;
		if (!dvdz->put_rec(&fieldNodes[iDim*jDim*kDim*31], rec)) 
			return NC_ERR;
		if (!dwdz->put_rec(&fieldNodes[iDim*jDim*kDim*32], rec)) 
			return NC_ERR;

	}
	
	// The file is automatically closed by the destructor. This frees
	// up any internal netCDF resources associated with the file, and
	// flushes any buffers.
	delete[] lats;
	delete[] lons;
	delete[] levs;
	
	return true;
	
}

bool CostFunction3D::writeAsi(const QString& asiFileName)
{
	// Initialize header
	int id[511];
	for (int n = 1; n <= 510; n++) {
		id[n]=-999;
	}
	
	// Calculate headers
	QStringList fieldNames;
	fieldNames  << "U" << "V" << "W" << "WS" << "RH"<< "HP" << "QP" << "RP" << "TP" << "PP" << "VO" << "DV" << "OW" << "S" << "PW"
	<< "MU" << "MV" << "MW" << "RO" << "PS" << "TK" << "QV" << "HH" << "DZ";
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

void CostFunction3D::obAdjustments() {
	
	// Load the obs locally and weight the nonlinear observation operators by interpolated bg fields
	for (int m = 0; m < mObs; m++) {
		int mi = m*14;
		for (int ob = 0; ob < 14; ob++) {
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
			for (int jjNode = jj-1; jjNode <= jj+2; ++jjNode) {
				for (int kNode = kk-1; kNode <= kk+2; ++kNode) {
                    int iNode = iiNode;
                    if ((iBCL[4] == PERIODIC) and (iNode < 0)) iNode = iDim-1;
                    if ((iBCL[4] == PERIODIC) and (iNode > (iDim-1))) iNode = iiNode - iDim;
                    if ((iNode < 0) or (iNode >= iDim)) continue;
                    
                    int jNode = jjNode;
                    if ((jBCL[4] == PERIODIC) and (jNode < 0)) jNode = jDim-1;
                    if ((jBCL[4] == PERIODIC) and (jNode > (jDim-1))) jNode = jjNode - jDim;
                    if ((jNode < 0) or (jNode >= jDim)) continue;
   
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
				// There is no contribution from this node 
				return 0;
			case 5:
				// There is no contribution from this node 
				return 0;
			case 6:
				// There is no contribution from this node 
				return 0;	
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
				node = -1;
				coeffmod = 1.;
				break;
			case 5:
				node = -1;
				coeffmod = -1.;
				break;				
			case 6:
				// There is no contribution from this node 
				return 0;	
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
				// There is no contribution from this node 
				return 0.;
			case 5:
				// There is no contribution from this node 
				return 0.;
			case 6:
				// There is no contribution from this node 
				return 0.;
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
	if ((BL == 4) and (m == 1)) {
		node = 0;
		coeffmod = -0.5;
	} else if ((BR == 4) and (m == M-1)) {
		node = M;
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


real CostFunction3D::getReferenceVariable(const int& refVariable, const real& heightm, const int& dz)
{
	real qvbhypcoeff[5];
	real rhoacoeff[5];
	real dpdzcoeff[5];

	if (referenceState == dunion_mt) {
		qvbhypcoeff[0] = 9.4826;
		qvbhypcoeff[1] = -0.0026721;
		qvbhypcoeff[2] = 2.8312e-07;
		qvbhypcoeff[3] = -1.3217e-11;
		qvbhypcoeff[4] = 2.2749e-16;
		
		rhoacoeff[0] = 1.1439;
		rhoacoeff[1] = -0.00010117;
		rhoacoeff[2] = 3.2486e-09;
		rhoacoeff[3] = -3.4898e-14;
		rhoacoeff[4] = -2.6925e-19;
		
		dpdzcoeff[0] = -11.432;
		dpdzcoeff[1] = 0.0010633;
		dpdzcoeff[2] = -4.0545e-08;
		dpdzcoeff[3] =  7.9634e-13;
		dpdzcoeff[4] = -5.8778e-18;
		
	}
	
	if (refVariable == qvbhypref) {
		real qvbhyp = 0.;
		for (int i = dz; i < 5; i++) {
			real power = pow(heightm, i-dz);
			if (dz) {
				qvbhyp += qvbhypcoeff[i] * power * i;
			} else {
				qvbhyp += qvbhypcoeff[i] * power;
			}
		}
		if (!dz) {
			if (qvbhyp < 0.) qvbhyp = 0.;
		}
		return qvbhyp;
	} else if (refVariable == rhoaref) {
		real rhoa = 0.;
		for (int i = dz; i < 5; i++) {
			real power = pow(heightm, i-dz); 
			if (dz) {
				rhoa += rhoacoeff[i] * power * i;
			} else {
				rhoa += rhoacoeff[i] * power;
			}
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
		if (dz) {
			real rhoadz = 0.;
			real qvdz = 0.;
			for (int i = dz; i < 5; i++) {
				real power = pow(heightm, i-dz); 
				rhoadz += rhoacoeff[i] * power * i;
				qvdz += qvbhypcoeff[i] * power * i;
			}
			rho = rhoadz * (1 + qv/1000.) + rhoa * (qvdz/500.);
			return rho;
		}
		return rho;
	} else if ((refVariable == href) or (refVariable == tempref) or (refVariable == pressref)) {
		// Integrate hydrostatic equation to get pressure and/or solve for T or h
		real press = 0.;
		real temp = 0.;
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
		if (dz) {
			real rhoadz = 0.;
			real qvdz = 0.;
			real dpdz = 0.;
			for (int i = dz; i < 5; i++) {
				real power = pow(heightm, i-dz); 
				rhoadz += rhoacoeff[i] * power * i;
				qvdz += qvbhypcoeff[i] * power * i;
				dpdz += dpdzcoeff[i] * power;
			}
			real alphadz = 1/(rhoadz * (286.9 + 461.5*qv/1000.) + 461.5 * rhoa * (qvdz/500.));
			real dtdz = press*alphadz + dpdz/(286.9*rhoa + 461.5*rhoa*qv/1000.);
			real dhdz = 1005.7*dtdz + 9.81 + 2.5e3*qvdz;
			switch (refVariable) {
				case href:
					return dhdz;
				case tempref:
					return dtdz;
				case pressref:
					return dpdz;
				default:
					break;
			}
			
		}
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

