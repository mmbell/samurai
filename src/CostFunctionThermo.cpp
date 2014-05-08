/*
 *  CostFunctionThermo.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionThermo.h"
#include <cmath>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QFile>
#include <QDir>
#include <QDateTime>
#include <netcdfcpp.h>
#include <GeographicLib/TransverseMercatorExact.hpp>

CostFunctionThermo::CostFunctionThermo(const int& numObs, const int& stateSize)
	: CostFunction3D(numObs, stateSize)
{   
}

CostFunctionThermo::~CostFunctionThermo()
{
}

void CostFunctionThermo::finalize()
{	    
}

void CostFunctionThermo::initialize(const QHash<QString, QString>* config, real* bgU, real* obs, ReferenceState* ref)
{
	// Initialize number of variables
	configHash = config;

	/* Set the output path */
	outputPath.setPath(configHash->value("output_directory"));
	
	// Horizontal boundary conditions
	iBCL[0] = bcHash.value(configHash->value("i_pip_bcL"));
    iBCR[0] = bcHash.value(configHash->value("i_pip_bcR"));
	iBCL[1] = bcHash.value(configHash->value("i_thetarhop_bcL"));
    iBCR[1] = bcHash.value(configHash->value("i_thetarhop_bcR"));
	iBCL[2] = bcHash.value(configHash->value("i_ftheta_bcL"));
    iBCR[2] = bcHash.value(configHash->value("i_ftheta_bcR"));

	jBCL[0] = bcHash.value(configHash->value("j_pip_bcL"));
    jBCR[0] = bcHash.value(configHash->value("j_pip_bcR"));
	jBCL[1] = bcHash.value(configHash->value("j_thetarhop_bcL"));
    jBCR[1] = bcHash.value(configHash->value("j_thetarhop_bcR"));
	jBCL[2] = bcHash.value(configHash->value("j_ftheta_bcL"));
    jBCR[2] = bcHash.value(configHash->value("j_ftheta_bcR"));

	kBCL[0] = bcHash.value(configHash->value("k_pip_bcL"));
    kBCR[0] = bcHash.value(configHash->value("k_pip_bcR"));
	kBCL[1] = bcHash.value(configHash->value("k_thetarhop_bcL"));
    kBCR[1] = bcHash.value(configHash->value("k_thetarhop_bcR"));
	kBCL[2] = bcHash.value(configHash->value("k_ftheta_bcL"));
    kBCR[2] = bcHash.value(configHash->value("k_ftheta_bcR"));

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
		
	// Adjust the internal, variable domain to include boundaries
	adjustInternalDomain(1);
    
	// Define nodes with internal domain
	int nodes = iDim*jDim*kDim;
    
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
	obsVector = new real[mObs*(obMetaSize+varDim*derivDim)];
	HCq = new real[mObs+nodes];
	innovation = new real[mObs+nodes];	
	bgState = new real[nState];
	bgStdDev = new real[nState];
	stateA = new real[nState];
	stateB = new real[nState];
	stateC = new real[nState];
    
    if (iBCL[0] == PERIODIC) {
        iLDim = iDim-2;   
    } else {
        iLDim = 4;
    }
    if (jBCL[0] == PERIODIC) {
        jLDim = jDim-2;   
    } else {
        jLDim = 4;
    }
    kLDim = 4;
    for (int var = 0; var < varDim; ++var) {
        iRank[var] = iDim - rankHash[iBCL[var]] - rankHash[iBCR[var]];
        if (iBCL[var] == PERIODIC) iRank[var]--;
        jRank[var] = jDim - rankHash[jBCL[var]] - rankHash[jBCR[var]];
        if (jBCL[var] == PERIODIC) jRank[var]--;
        kRank[var] = kDim - rankHash[kBCL[var]] - rankHash[kBCR[var]];
        // Need to verify that Rank is sufficient (>2?)
        iL[var] = new real[iRank[var]*iLDim];
        jL[var] = new real[jRank[var]*jLDim];
        kL[var] = new real[kRank[var]*kLDim];
        iGamma[var] = new real[iRank[var]*iDim];
        jGamma[var] = new real[jRank[var]*jDim];
        kGamma[var] = new real[kRank[var]*kDim];

    }
	// Precalculate the basis functions for lookup table option
	basisappx = configHash->value("spline_approximation").toInt();
	if (basisappx > 0) {
		basis0 = new real[2000000];
		basis1 = new real[2000000];
		fillBasisLookup();
	}
    
    // Initialize the Fourier transforms
    iFFTin = (double*) fftw_malloc(sizeof(double) * iDim);
    iFFTout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * iDim);
    iForward = fftw_plan_dft_r2c_1d(iDim, iFFTin, iFFTout, FFTW_MEASURE);
    iBackward = fftw_plan_dft_c2r_1d(iDim, iFFTout, iFFTin, FFTW_MEASURE);
    jFFTin = (double*) fftw_malloc(sizeof(double) * jDim);
    jFFTout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * jDim);
    jForward = fftw_plan_dft_r2c_1d(jDim, jFFTin, jFFTout, FFTW_MEASURE);
    jBackward = fftw_plan_dft_c2r_1d(jDim, jFFTout, jFFTin, FFTW_MEASURE);
    iMaxWavenumber = configHash->value("i_max_wavenumber").toFloat();
    jMaxWavenumber = configHash->value("j_max_wavenumber").toFloat();
}

void CostFunctionThermo::initState(const int iteration)
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
	bgError[0] = configHash->value("bg_pip_error").toFloat();
	bgError[1] = configHash->value("bg_thetarhop_error").toFloat();
	bgError[2] = configHash->value("bg_ftheta_error").toFloat();
	bgError[3] = 0;
	bgError[4] = 0;
	bgError[5] = 0;
	bgError[6] = 0;
	
	// Initialize filter scales
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
	//obAdjustments();
	obsVector= rawObs;

	// d = y - HXb
	calcInnovation();
	
	// Output the original background field
	outputAnalysis("background", bgState);

	cout << "Beginning analysis...\n";
		
	// HTd
	calcHTranspose(innovation, stateC);
	
	SCtranspose(stateC, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtranspose(stateA, CTHTd);
				
}	
