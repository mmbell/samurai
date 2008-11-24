/*
 *  CostFunctionAnalytic.cpp
 *  tcvar
 *
 *  Copyright 2008  Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionAnalytic.h"
#include <iostream>

CostFunctionAnalytic::CostFunctionAnalytic()
	: CostFunction()
{
}

CostFunctionAnalytic::~CostFunctionAnalytic()
{
}

void CostFunctionAnalytic::initialize()
{
	nState = 3;
	// Allocate and fill the current state vector
	currState = new double[nState];
	currGradient = new double[nState];
	tempState = new double[nState];
	tempGradient = new double[nState];
	for (int n = 0; n < nState; n++) {
		currState[n] = 0.0;
		currGradient[n] = 0.0;
		tempState[n] = 0.0;
		tempGradient[n] = 0.0;
	}
	currState[0] = 5;
	currState[1] = 7;
	currState[2] = 10;
	mObs = 0;
}

void CostFunctionAnalytic::finalize()
{
	cout << "Final State: x = " << currState[0] 
	<< " , y = " << currState[1] << " , z = " << currState[2] << endl;
}

double CostFunctionAnalytic::funcValue(double* state)
{
	
	double x = state[0];
	double y = state[1];
	double z = state[2];
	double parab = 10*(x-1)*(x-1) + 20*(y-2)*(y-2) + 30*(z-3)*(z-3) + 40;
	return parab;
	
}

void CostFunctionAnalytic::funcGradient(double* state, double* gradient)
{
	double x = state[0];
	double y = state[1];
	double z = state[2];
	gradient[0] = 2*10*(x-1);
	gradient[1] = 2*20*(y-2);
	gradient[2] = 2*30*(z-3);
}
