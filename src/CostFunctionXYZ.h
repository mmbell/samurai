/*
 *  CostFunctionXYZ.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCXYZ_H
#define COSTFUNCXYZ_H

#include "CostFunction3D.h"
#include "MetObs.h"
#include "VarDriver.h"
#include "HashMap.h"
#include <string>


class CostFunctionXYZ: public CostFunction3D
{

public:
  CostFunctionXYZ(const Projection& proj, const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionXYZ();

	// VarDriver *vardriver;

private:
	bool outputAnalysis(const std::string& suffix, real* Astate);
	bool outputAnalysis_thermo(const std::string& suffix, real* Astate);
	bool writeAsi(const std::string& asiFileName);
	bool writeNetCDF(const std::string& netcdfFileName);
	bool SItransform(size_t numVars, double *finalAnalysis, double *mishData, real *Astate, ofstream *outStream);

	bool fractl_mode;
	// HashMap configHash;
};

#endif
