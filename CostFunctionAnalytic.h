/*
 *  CostFunctionAnalytic.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */


#ifndef COSTFUNCANALYTIC_H
#define COSTFUNCANALYTIC_H

#include "CostFunction.h"


class CostFunctionAnalytic : public CostFunction
{
		
public:
	CostFunctionAnalytic();
	~CostFunctionAnalytic();
	void initialize();
	void finalize();

private:
	double funcValue(double* state);
	void funcGradient(double* state, double* gradient);

};

#endif
