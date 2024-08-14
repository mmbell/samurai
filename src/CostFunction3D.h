/*
 *  CostFunction3D.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNC3D_H
#define COSTFUNC3D_H

#include "CostFunction.h"
#include "BSpline.h"
#include "RecursiveFilter.h"
#include "Observation.h"
#include "ReferenceState.h"
#include "VarDriver.h" // added
#include "ErrorData.h"
#include "HashMap.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <string>
#include <unordered_map>
//#include <complex.h>

// Whether we use FFTW or cuFFT's FFTW interface depends on whether we have the 'USE_CUFFTW' define set (Cmake):
#ifdef USE_CUFFTW
  #include <cufftw.h>
#else
  #include <fftw3.h>
#endif


class CostFunction3D : public CostFunction
{

public:

CostFunction3D(const Projection& proj, const int& numObs = 0, const int& stateSize = 0);
	virtual ~CostFunction3D();
    void initialize(HashMap* config, real* bgU, real* obs, ReferenceState* ref);
	void finalize();
	void updateBG();
	void initState(const int iteration);
	bool copyResults(int iDim, int jDim, int kDim,
			 float *u, float *v, float *w, float *th, float *p);

protected:
	double funcValue(double* state);
	void funcGradient(double* state, double* gradient);
	void updateHCq(double* state);
	double funcValueAndGradient(double *state, double *gradient);
	void funcHessian(double *x, double *hessian);
	void updateHCq(double* state, double* HCq);
	real Basis(const int& m, const real& x, const int& M,const real& xmin,
			   const real& DX, const real& DXrecip, const int& derivative,
			   const int& BL, const int& BR, const real& lambda = 0);
	real BasisBC(real b, const int& m, const real& x, const int& M,const real& xmin,
			   const real& DX, const real& DXrecip, const int& derivative,
			   const int& BL, const int& BR, const real& lambda = 0);
	void fillBasisLookup();
	bool filterArray(real* array, const int& arrLength);
	bool setupSplines();
	void obAdjustments();
	void solveBC(real* A, real* B);
	bool SAtransform(const real* Bstate, real* Astate);
	bool SAtranspose(const real* Astate, real* Bstate);
	void calcInnovation();
	void calcHTranspose(const real* yhat, real* Astate);
	virtual bool outputAnalysis(const std::string& suffix, real* Astate) = 0;
	void SBtransform(const real* Ustate, real* Bstate);
	void SBtranspose(const real* Bstate, real* Ustate);
	void SCtransform(const real* Astate, real* Cstate);
	void SCtranspose(const real* Cstate, real* Astate);
	void FFtransform(const real* Astate, real* Cstate);

	bool writeAsi(const std::string& asiFileName);
	bool writeNetCDF(const std::string& netcdfFileName);

	void adjustInternalDomain(int increment);
	void calcSplineCoefficients(const int& Dim, const real& eq, const int* BCL, const int* BCR,
                                const real& xmin, const real& DX, const real& DXrecip, const int& LDim,
                                real* L[7], real* gamma[7]);
	bool copy3DArray(real *src, float *dest, int iDim, int jDim, int kDim);
	void calcHmatrix();
	void Htransform(const real* Cstate, real* Hstate);

	// A couple of utilities functions to help query config values
  bool isTrue(const char *flag_in) {
    std::string flag = flag_in;
    if ((*configHash)[flag] == "true") {
      return true;
    }
    return false;
  }

	bool isEqual(const char *flag_in, std::string value) {
    std::string flag = flag_in;
    if ((*configHash)[flag] == value) {
      return true;
    }
    return false;
	}

	void initBkgdErrors();

	bool mishFlag;
	int iDim, jDim, kDim;
	int iLDim, jLDim, kLDim;
	int iRank[7], jRank[7], kRank[7];
  int kRankMax;
	real iMin, iMax, DI, DIrecip;
	real jMin, jMax, DJ, DJrecip;
	real kMin, kMax, DK, DKrecip;
	real* bgFields;
	real* bgState;
	real* bgStdDev;
	real* obsVector;
        real* obsData;  // This only contains the necessary data for the calcHTranspose subroutine
	real* rawObs;
	real* stateA;
	real* stateB;
	real* stateC;
	real* stateU;
	real* CTHTd;
	real* HCq;
	real* innovation;
	real* iL[7];
	real* jL[7];
	real* kL[7];
	real* iGamma[7];
	real* jGamma[7];
	real* kGamma[7];
        real* kGammaL;
	real* kLL;
	real* finalAnalysis;

	uint64_t varDim; // NCAR: promoted to 64-bit, since it should auto-promote calculations with it to 64-bit
	int derivDim;
	int obMetaSize;

	real bgError[7];
	int iBCL[7], iBCR[7], jBCL[7], jBCR[7], kBCL[7], kBCR[7];
	int derivative[4][3];
	real constHeight;
	real mcWeight;
	real latReference, lonReference;
	int iMaxWavenumber[7], jMaxWavenumber[7], kMaxWavenumber[7];
	fftw_plan iForward, jForward, iBackward, jBackward, kForward, kBackward;
	double *iFFTin, *jFFTin, *kFFTin;
	fftw_complex *iFFTout, *jFFTout, *kFFTout;
        bool UseFFT;

	// explicitly store the H matrix in CSR format
	real *H;
	uint64_t *IH; // Array with extent:(nState+1) can take on values [0 to nonzeros]
        uint32_t *JH; // Array with extent:(nonzeros) can take on values [0 to nState-1]

	// explicity store the H^t matrix in CSR format
        real *Ht;
	uint64_t *IHt; // Array with extent(mObs+1_ can take on values [0 to nonzeros]
	uint32_t *JHt; // Array with extent(nonzeros) can take on values [0 to mObs-1]

	int basisappx;
	real* basis0;
	real* basis1;
	HashMap* configHash;
	std::unordered_map<std::string, int> bcHash;
	std::unordered_map<int, int> rankHash;

	enum BoundaryConditionTypes {
		RX = -1,
		R0 = 0,
		R1T0 = 1,
		R1T1 = 2,
		R1T2 = 3,
		R1T10 = 4,
		R2T10 = 5,
		R2T20 = 6,
		R3 = 7,
		PERIODIC = 8
	};

	enum BasisApproximation {
		NONE = 0,
		PARTIAL = 1,
		FULL = 2
	};

	real iFilterScale,jFilterScale, kFilterScale;
	RecursiveFilter* iFilter;
	RecursiveFilter* jFilter;
	RecursiveFilter* kFilter;

	ReferenceState* refstate;
	std::string dataPath, outputPath;

	ErrorData variance;
	#pragma acc declare copyin(iDim,jDim,kDim,varDim,kLDim)
	#pragma acc declare copyin(iFilterScale,jFilterScale,kFilterScale)
	#pragma acc declare copyin(kRankMax)
};

#endif
