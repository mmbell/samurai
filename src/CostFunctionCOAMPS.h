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
	CostFunctionCOAMPS(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionCOAMPS();

private:
	bool outputAnalysis(const QString& suffix, real* Astate);
	bool writeAsi(const QString& asiFileName);
	bool writeNetCDF(const QString& netcdfFileName);
	bool writeFlatfile(const QString& flatFileName, const int var);
	int sDim;
	const real sigma[42] = { 31885.0, 29385.0, 25085.0, 21985.0, 19685.0, 17885.0,
		16385.0, 15072.5, 13910.0, 12860.0, 11885.0, 10955.0, 10065.0,
		9215.0, 8405.0, 7635.0, 6905.0, 6215.0, 5565.0, 4955.0, 4385.0,
		3855.0, 3365.0, 2915.0, 2505.0, 2135.0, 1805.0, 1515.0, 1265.0,
		1050.0, 860.0, 690.0, 540.0, 410.0, 300.0, 210.0, 140.0, 90.0,
		55.0, 30.0, 10.0, 0.0 };
};

#endif
