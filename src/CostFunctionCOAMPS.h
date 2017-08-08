/*
 *  CostFunctionCOAMPS.h
 *  samurai
 *
 *  Copyright 2015 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCCOAMPS_H
#define COSTFUNCCOAMPS_H

#include "CostFunction3D.h"
#include <QString>

class CostFunctionCOAMPS: public CostFunction3D
{

public:
	CostFunctionCOAMPS(const Projection& proj, const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionCOAMPS();

private:
	bool outputAnalysis(const QString& suffix, real* Astate);
	bool writeAsi(const QString& asiFileName);
	bool writeNetCDF(const QString& netcdfFileName);
	bool writeFlatfile(const QString& flatFileName, const int var);
	bool copyResults(float *u, float *v, float *w, float *th, float *p);
	int sDim;
};

#endif
