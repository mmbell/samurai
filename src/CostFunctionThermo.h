/*
 *  CostFunctionThermo.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNCTHERMO_H
#define COSTFUNCTHERMO_H

#include "CostFunction3D.h"
#include "CostFunctionRTZ.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <QString>
#include <QHash>
#include <QDir>
//#include <complex.h>
#include <fftw3.h>

class CostFunctionThermo: public CostFunction3D
{
	
public:
	CostFunctionThermo(const int& numObs = 0, const int& stateSize = 0);
	~CostFunctionThermo();
    void initialize(const QHash<QString, QString>* config, real* bgU, real* obs, ReferenceState* ref); 
	void finalize();
	void initState(const int iteration);
	
protected:
	bool outputAnalysis(const QString& suffix, real* Astate);
	bool writeAsi(const QString& asiFileName);
	bool writeNetCDF(const QString& netcdfFileName);

};

#endif
