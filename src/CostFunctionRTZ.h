/*
 *  CostFunctionRTZ.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCRTZ_H
#define COSTFUNCRTZ_H

#include "CostFunction3D.h"
#include <QString>

class CostFunctionRTZ: public CostFunction3D
{
	
public:
  CostFunctionRTZ(const Projection& proj, const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionRTZ();
	
private:
	bool outputAnalysis(const QString& suffix, real* Astate);
	bool writeAsi(const QString& asiFileName);
	bool writeNetCDF(const QString& netcdfFileName);
};

#endif
