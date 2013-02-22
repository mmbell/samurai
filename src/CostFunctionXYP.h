/*
 *  CostFunctionXYP.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCXYP_H
#define COSTFUNCXYP_H

#include "CostFunctionXYZ.h"
#include <QString>

class CostFunctionXYP: public CostFunction3D
{
	
public:
	CostFunctionXYP(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionXYP();
	
private:
	bool outputAnalysis(const QString& suffix, real* Astate);
	bool writeAsi(const QString& asiFileName);
	bool writeNetCDF(const QString& netcdfFileName);
	real pMin, pMax, DP;
	int pDim;
};

#endif
