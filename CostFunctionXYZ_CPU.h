/*
 *  CostFunctionXYZ_CPU.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCXYZCPU_H
#define COSTFUNCXYZCPU_H

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

class CostFunctionXYZ_CPU: public CostFunction
{
	
public:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	CostFunctionXYZ_CPU(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionXYZ_CPU();
    void initialize(const real& imin, const real& imax, const int& idim,
					const real& imin, const real& imax, const int& jdim,
					const real* iA, const real* jA, real* bgU, real* obs, 
					const unsigned int* IXdim = 0, const vector<real>* IX = 0,
					const unsigned int* JXdim = 0, const vector<real>* JX = 0);
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
	float BasisOri(int m, float x, int M, float xmin, 
				float DX, float DXrecip, int C);
	float DBasisOri(int m, float x, int M, float xmin, 
				 float DX, float DXrecip, int C);
	float DDBasisOri(int m, float x, int M, float xmin, 
				 float DX, float DXrecip, int C);
	float DDDBasisOri(int m, float x, int M, float xmin, 
				 float DX, float DXrecip, int C);
	bool filterArray(real* array, const int& arrLength);
	bool setupSplines();
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
	int iDim;
	const unsigned int* IXDim;
	real iMin, iMax, DI, DIrecip;
	int jDim;
	const unsigned int* JXDim;
	real jMin, jMax, DJ, DJrecip;
	const vector<real>* IXform;
	const vector<real>* JXform;
	real* bgFields;
	real* bgState;
	real* obsVector;
	real* rawObs;
	real* stateA;
	real* stateB;
	real* stateC;
	real* stateU;
	real* Uprime;
	real* CTHTd;
	real* HCq;
	real* innovation;
	const real* iA;
	const real* jA;
	real* iTemp;
	real* jTemp;
	real* iL;
	real* jL;
	int varDim;
	int bState;
	real bgError[5];
	real bgErrorScale;
	real LI, LJ;
	real LF[5];
	int bcLeft[5], bcRight[5], bcTop[5], bcBottom[5];
	
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
	real iFilterCoeff[5];
	real iFilterBC[5][5];
	real jFilterCoeff[5];
	real jFilterBC[5][5];

	real rhoBase;
	real rhoInvScaleHeight;
	real qBase;
	real qInvScaleHeight;
	
};

#endif
