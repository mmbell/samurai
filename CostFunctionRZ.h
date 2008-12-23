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
	typedef BSplineBase<double> SplineBase;
	typedef BSpline<double> SplineD;
	CostFunctionRZ(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionRZ();
    void initialize(SplineD* vecs, SplineD* scalars, SplineD* ctrls, SplineD* zs, SplineD* zpsi, 
					vector<double>* bgr, vector<double>* bgz, vector<double>** bgf, 
					unsigned int* ctrlR, vector<double>* RX, Observation* obs);
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
	vector<double>* bgradii;
	vector<double>* bgheights;
	vector<double>** bgFields;
	unsigned int* RState;
	unsigned int maxRadius;
	unsigned int pState;
	unsigned int zState;
	vector<double>* RXform;
	vector<double>* rXform;
	Observation* obsVector;
	RecursiveFilter* filterR;
	RecursiveFilter* filterZ;

	double* HT;
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
	void updateHCq_GPU(double* state);
#ifdef CSTGPU
	float* HCq;
#else
	float* HCq;
	//double* HCq;
#endif
	
};

#endif
