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

class CostFunction3D: public CostFunction
{
	
public:
	CostFunction3D(const int& numObs = 0, const int& stateSize = 0);
	virtual ~CostFunction3D();
    void initialize(const QHash<QString, QString>* config, real* bgU, real* obs, ReferenceState* ref); 
	void finalize();
	void updateBG();
	void initState(const int iteration);
	
protected:
	double funcValue(double* state);
	void funcGradient(double* state, double* gradient);
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
                                real* L[7], real* gamma[7]);
	bool outputMish;
	int iDim, jDim, kDim;
    int iLDim, jLDim, kLDim;
    int iRank[7], jRank[7], kRank[7];
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
	real* iL[7];
	real* jL[7];
	real* kL[7];
    real* iGamma[7];
    real* jGamma[7];
    real* kGamma[7];
	real* finalAnalysis;
	int varDim;
    int derivDim;
	real bgError[7];
	int iBCL[7], iBCR[7], jBCL[7], jBCR[7], kBCL[7], kBCR[7];
    int derivative[4][3];
	real constHeight;
	real mcWeight;

	real* basis0;
	real* basis1;
	const QHash<QString, QString>* configHash;
	QHash<QString, int> bcHash;
    QHash<int, int> rankHash;
    
	enum BoundaryConditionTypes {
		R0 = -1,
		R1T0 = 0,
		R1T1 = 1,
		R1T2 = 2,
		R1T10 = 3,
		R2T10 = 4,
		R2T20 = 5,
		R3 = 6,
		PERIODIC = 7
	};

	real iFilterScale,jFilterScale, kFilterScale;
	RecursiveFilter* iFilter;
	RecursiveFilter* jFilter;
	RecursiveFilter* kFilter;
	
    ReferenceState* refstate;
    
	
};

#endif
