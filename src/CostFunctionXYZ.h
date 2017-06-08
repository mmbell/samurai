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
#include <QString>

class CostFunctionXYZ: public CostFunction3D
{
	
public:
  CostFunctionXYZ(const Projection& proj, const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionXYZ();
	
private:
	bool outputAnalysis(const QString& suffix, real* Astate);
	bool writeAsi(const QString& asiFileName);
	bool writeNetCDF(const QString& netcdfFileName);
};

#endif
