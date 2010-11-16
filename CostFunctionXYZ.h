/*
 *  CostFunctionXYZ.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCXYZ_H
#define COSTFUNCXYZ_H

#include "CostFunction.h"
#include "BSpline.h"
#include "RecursiveFilter.h"
#include "Observation.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <QString>
#include <QHash>

class CostFunctionXYZ: public CostFunction
{
	
public:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	CostFunctionXYZ(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionXYZ();
    void initialize(const QHash<QString, QString>& config, real* bgU, real* obs); 
	void finalize();
	void updateBG();
	void initState();
	
private:
	double funcValue(double* state);
	void funcGradient(double* state, double* gradient);
	void updateHCq(double* state);
	real FullBasis(int m, real x, int M, real xmin, 
				real DX, real DXrecip, int derivative,
				int BL, int BR, real lambda = 0);
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
	bool SAtransform_ori(real* Bstate, real* Astate);
	void calcInnovation();
	void calcHTranspose(const real* yhat, real* Astate);
	bool outputAnalysis(const QString& suffix, real* Astate, bool updateMish);
	void SBtransform(const real* Ustate, real* Bstate);
	void SBtranspose(const real* Bstate, real* Ustate);
	void SCtransform(const real* Astate, real* Cstate);
	void SCtranspose(const real* Cstate, real* Astate);
	
	real getReferenceVariable(const int& refVariable, const real& heightm, const int& dz = 0);
	real bhypTransform(real qv);
	real bhypInvTransform(real qvbhyp);
	void writeAsi();
	bool writeNetCDF(const QString& netcdfFile);
	
	int iDim, jDim, kDim;
	real iMin, iMax, DI, DIrecip;
	real jMin, jMax, DJ, DJrecip;
	real kMin, kMax, DK, DKrecip;
	real* bgFields;
	real* bgState;
	real* obsVector;
	real* rawObs;
	real* stateA;
	real* stateB;
	real* stateC;
	real* stateU;
	real* CTHTd;
	real* HCq;
	real* innovation;
	real* iL;
	real* jL;
	real* kL;
	real* kLw;
	real* fieldNodes;
	int varDim;
	real bgError[7];
	int iBCL[7], iBCR[7], jBCL[7], jBCR[7], kBCL[7], kBCR[7];
	real bgErrorScale;
	real constHeight;
	real mcWeight;
	int referenceState;
	real basis0[200000];
	real basis1[200000];
	QHash<QString, QString> configHash;

	float BoundaryConditions[9][4];
	enum BoundaryConditionTypes {
		R1T0 = 0,
		R1T1 = 1,
		R1T2 = 2,
		R1T10 = 3,
		R2T10 = 4,
		R2T20 = 5,
		R3 = 6,
		R3X = 7
	};
	
	RecursiveFilter* iFilter;
	RecursiveFilter* jFilter;
	RecursiveFilter* kFilter;

	real rhoBase;
	real rhoInvScaleHeight;
	real qBase;
	real qInvScaleHeight;
	
	enum referenceStates {
		jordan
	};
	
	enum referenceVariables {
		qvbhypref,
		rhoaref,
		rhoref,
		href,
		tempref,
		pressref
	};	
	
};

#endif
