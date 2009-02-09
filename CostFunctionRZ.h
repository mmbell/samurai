/*
 *  CostFunctionRZ.h
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCRZ_H
#define COSTFUNCRZ_H

#include "CostFunction.h"
#include "BSpline.h"
#include "ParametricVortex.h"
#include "RecursiveFilter.h"
#include "Observation.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>

// Define this to run on GPU
//#define CSTGPU

class CostFunctionRZ: public CostFunction
{
	
public:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	CostFunctionRZ(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionRZ();
    void initialize(SplineD* vecs, SplineD* scalars, SplineD* ctrls, SplineD* zs, SplineD* zpsi, 
					vector<real>* bgr, vector<real>* bgz, vector<real>** bgf, 
					unsigned int* ctrlR, vector<real>* RX, Observation* obs);
	void finalize();
	// Should this return an object, not a pointer?
	void getCq(double* Cq);
	void initState();
	
private:
	double funcValue(double* state);
	void funcGradient(double* state, double* gradient);
	void updateHCq(double* state);
	SplineD* vecSpline;
	SplineD* scalarSpline;
	SplineD* ctrlSpline;
	SplineD* CqSpline;
	SplineD* zSpline;
	SplineD* zSplinePsi;
	vector<real>* bgradii;
	vector<real>* bgheights;
	vector<real>** bgFields;
	unsigned int* RState;
	unsigned int maxRadius;
	unsigned int pState;
	unsigned int zState;
	vector<real>* RXform;
	vector<real>* rXform;
	Observation* obsVector;
	RecursiveFilter* filterR;
	RecursiveFilter* filterZ;

	real* HT;
	double* zHT;
	double* stateMod;
	double* innovation;
	double* rCTHT;
	double* zCTHT;
	double** rzCTHT;
	double** rzSTHT;
	double*** CTHTd;
	double* fieldR;
	double* fieldZ;
	double** field;
	double* zBuffer;
	unsigned int numVars;
	double varScale[6];
	double bgError[6];

	// GPU Variables
	float* coeffHost;
	float* obs_h;
	float BoundaryConditions[9][4];
	void updateHCq_GPU(double* state);
	void updateHCq_parallel(double* state);
	float Basis(int m, float x, int M, float xmin, 
				float DX, float DXrecip, float ONESIXTH, int C);
	float DBasis(int m, float x, int M, float xmin, 
				 float DX, float DXrecip, float ONESIXTH, int C);

#ifdef CSTGPU
	float* HCq;
#else
	float* HCq;
	//double* HCq;
#endif
	
};

#endif
