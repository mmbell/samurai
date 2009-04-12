/*
 *  CostFunctionRZ_CPU.h
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCRZCPU_H
#define COSTFUNCRZCPU_H

#include "CostFunction.h"
#include "BSpline.h"
#include "RecursiveFilter.h"
#include "Observation.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>

class CostFunctionRZ_CPU: public CostFunction
{
	
public:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	CostFunctionRZ_CPU(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionRZ_CPU();
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
	float Basis(int m, float x, int M, float xmin, 
				float DX, float DXrecip, int C);
	float DBasis(int m, float x, int M, float xmin, 
				 float DX, float DXrecip, int C);
	float DDBasis(int m, float x, int M, float xmin, 
				 float DX, float DXrecip, int C);
	float DDDBasis(int m, float x, int M, float xmin, 
				 float DX, float DXrecip, int C);
	bool filterArray(real* array, const int& arrLength);
	bool setupSplines();
	void solveBC(real* A, real* B);
	bool SAtransform(real* Bstate, real* Astate);
	bool SAtransform_ori(real* Bstate, real* Astate);

	void calcInnovation();
	void calcHTranspose(const real* yhat, real* Astate);

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
	real* stateA;
	real* stateB;
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
	double varScale[6];
	real bgError[6];

	float BoundaryConditions[9][4];

	RecursiveFilter* iFilter;
	RecursiveFilter* jFilter;
	real iFilterCoeff[5];
	real iFilterBC[5][5];
	real jFilterCoeff[5];
	real jFilterBC[5][5];

	real rhoBase;
	real rhoInvScaleHeight;
};

#endif
