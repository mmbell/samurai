/*
 *  CostFunctionR.h
 *  tcvar
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
	typedef BSplineBase<double> SplineBase;
	typedef BSpline<double> SplineD;
	CostFunctionR(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionR();
    void initialize(SplineD* bgs, vector<double>* bgr, vector<double>* bgf,
					SplineD* ctrl, vector<double>* ctrlR, 
					vector<double>* RXform, vector<double>* rXform,  Observation* ob);
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
	vector<double>* bgradii;
	vector<double>* bgFields;
	vector<double>* ctrlRadii;
	vector<double>* RXform;
	vector<double>* rXform;
	Observation* obsVector;
	RecursiveFilter* filter1d;
	double* HCq;
	double* HTHCq;
	double* HTd;
	double* stateMod;
	double* innovation;
	unsigned int numVars;
	float bgError[6];
};

#endif
