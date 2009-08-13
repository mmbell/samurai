/*
 *  CostFunctionR.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCR_H
#define COSTFUNCR_H

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

class CostFunctionR : public CostFunction
{
	
public:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	CostFunctionR(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionR();
    void initialize(SplineD* bgs, vector<real>* bgr, vector<real>* bgf,
					SplineD* ctrl, vector<real>* ctrlR, 
					vector<real>* RXform, vector<real>* rXform,  Observation* ob);
	void finalize();
	// Should this return an object, not a pointer?
	void getCq(double* Cq);
	void initState();
	
private:
	double funcValue(double* state);
	void funcGradient(double* state, double* gradient);
	void updateHCq(double* state);
	SplineD* bgSpline;
	SplineD* ctrlSpline;
	SplineD* CqSpline;
	vector<real>* bgradii;
	vector<real>* bgFields;
	vector<real>* ctrlRadii;
	vector<real>* RXform;
	vector<real>* rXform;
	Observation* obsVector;
	RecursiveFilter* filter1d;
	real* HCq;
	real* HTHCq;
	real* HTd;
	double* stateMod;
	double* innovation;
	unsigned int numVars;
	float bgError[6];
};

#endif
