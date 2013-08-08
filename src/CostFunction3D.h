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
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <QString>
#include <QHash>
#include <QDir>
//#include <complex.h>
#include <fftw3.h>

class CostFunction3D: public CostFunction
{
	
public:
	CostFunction3D(const int& numObs = 0, const int& stateSize = 0);
	virtual ~CostFunction3D();
    virtual void initialize(const QHash<QString, QString>* config, real* bgU, real* obs, ReferenceState* ref); 
	virtual void finalize();
	void updateBG();
	virtual void initState(const int iteration);
	
protected:
	static const int varDim = 7;
    static const int derivDim = 4;
	virtual double funcValue(double* state);
	virtual void funcGradient(double* state, double* gradient);
	void updateHCq(double* state);
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
	virtual bool outputAnalysis(const QString& suffix, real* Astate) = 0;
	void SBtransform(const real* Ustate, real* Bstate);
	void SBtranspose(const real* Bstate, real* Ustate);
	void SCtransform(const real* Astate, real* Cstate);
	void SCtranspose(const real* Cstate, real* Astate);
	
	bool writeAsi(const QString& asiFileName);
	bool writeNetCDF(const QString& netcdfFileName);
	void adjustInternalDomain(int increment);
	void calcSplineCoefficients(const int& Dim, const real& eq, const int* BCL, const int* BCR,
                                const real& xmin, const real& DX, const real& DXrecip, const int& LDim,
                                real* L[varDim], real* gamma[varDim]);
	bool outputMish;
	int iDim, jDim, kDim;
    int iLDim, jLDim, kLDim;
    int iRank[varDim], jRank[varDim], kRank[varDim];
	real iMin, iMax, DI, DIrecip;
	real jMin, jMax, DJ, DJrecip;
	real kMin, kMax, DK, DKrecip;
	real* bgFields;
	real* bgState;
	real* bgStdDev;
	real* obsVector;
	real* rawObs;
	real* stateA;
	real* stateB;
	real* stateC;
	real* stateU;
	real* CTHTd;
	real* HCq;
	real* innovation;
	real* iL[varDim];
	real* jL[varDim];
	real* kL[varDim];
    real* iGamma[varDim];
    real* jGamma[varDim];
    real* kGamma[varDim];
	real* finalAnalysis;
	real bgError[varDim];
	int iBCL[varDim], iBCR[varDim], jBCL[varDim], jBCR[varDim], kBCL[varDim], kBCR[varDim];
    int derivative[4][3];
	real constHeight;
	real mcWeight;
    int iMaxWavenumber, jMaxWavenumber;
    double *iFFTin, *jFFTin;
    fftw_complex *iFFTout, *jFFTout;
    fftw_plan iForward, jForward, iBackward, jBackward;
    
	int basisappx;
	real* basis0;
	real* basis1;
	const QHash<QString, QString>* configHash;
	QHash<QString, int> bcHash;
    QHash<int, int> rankHash;
    
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
    QDir outputPath;
	
};

#endif
