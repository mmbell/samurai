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
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	CostFunction3D(const int& numObs = 0, const int& stateSize = 0);
	~CostFunction3D();
    void initialize(const QHash<QString, QString>& config, real* bgU, real* obs); 
	void finalize();
	void updateBG();
	void initState();
	
private:
	double funcValue(double* state);
	void funcGradient(double* state, double* gradient);
	void updateHCq(double* state);
	real Basis(int m, real x, int M, real xmin, 
				real DX, real DXrecip, int derivative,
				int BL, int BR, real lambda = 0);
	void fillBasisLookup();
	bool filterArray(real* array, const int& arrLength);
	bool setupSplines();
	void obAdjustments();
	void solveBC(real* A, real* B);
	bool SAtransform(real* Bstate, real* Astate);
	bool SAtransform_ori(real* Bstate, real* Astate);
	void calcInnovation();
	void calcHTranspose(const real* yhat, real* Astate);
	bool outputAnalysis(const QString& suffix, real* Astate, bool updateMish);
	void SBtransform(real* Ustate, real* Bstate);
	void SBtranspose(real* Bstate, real* Ustate);
	void SCtransform(real* Astate, real* Cstate);
	void SCtranspose(real* Cstate, real* Astate);
	
	real getReferenceVariable(int refVariable, real heightm);
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
	real* stateV;
	real* CTHTd;
	real* HCq;
	real* innovation;
	real* iL;
	real* jL;
	real* kL;
	real* kLw;
	real* fieldNodes;
	int varDim;
	int bState;
	real bgError[7];
	int iBCL[7], iBCR[7], jBCL[7], jBCR[7], kBCL[7], kBCR[7];
	real bgErrorScale;
	real constHeight;
	real mcWeight;
	int referenceState;
	real basis0[2000];
	real basis1[2000];
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
