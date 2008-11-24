/*
 *  VarDriver.h
 *  tcvar
 *
 *  Created by Michael Bell on 4/12/08.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVER_H
#define VARDRIVER_H

#include "BSpline.h"
#include "Observation.h"
#include "CostFunctionR.h"
#include "CostFunctionAnalytic.h"
#include "MetObs.h"
#include "TCcenter.h"
#include <iostream>
#include <vector>
#include <QHash>
#include <QDir>
#include <QList>
#include <QString>

using namespace std;

class VarDriver
{

public:
	VarDriver();
	virtual ~VarDriver();
	virtual bool run() = 0;

protected:
	
	float CoriolisF;
	float Pi;
	unsigned int numVars;
	unsigned int numHeights;
	vector<TCcenter> tcVector;

	// Data Processing
	QHash<QString, int> dataSuffix;
	enum dataFormats {
		unknown,
		frd,
		cls,
		sec,
		ten,
		swp
	};
	bool read_frd(QFile& metFile, QList<MetObs>* metObVector);
	bool read_cls(QFile& metFile, QList<MetObs>* metObVector);
	bool read_sec(QFile& metFile, QList<MetObs>* metObVector);
	bool read_ten(QFile& metFile, QList<MetObs>*metObVector);
	bool read_dorade(QFile& metFile, QList<MetObs>* metObVector);
	bool readTCcenters();
	
};

#endif
