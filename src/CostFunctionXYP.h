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
#include <string>

class CostFunctionXYP: public CostFunction3D
{
	
public:
  CostFunctionXYP(const Projection& proj, const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionXYP();
	
private:
	bool outputAnalysis(const std::string& suffix, real* Astate);
	bool outputAnalysis_thermo(const std::string& suffix, real* Astate);
	bool writeAsi(const std::string& asiFileName);
	bool writeNetCDF(const std::string& netcdfFileName);
	real pMin, pMax, DP;
	int pDim;
};

#endif
