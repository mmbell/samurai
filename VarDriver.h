/*
 *  VarDriver.h
 *  samurai
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
#include <QDomDocument>
#include "precision.h"

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
	float rhoBase;
	float rhoInvScaleHeight;
	unsigned int numVars;
	unsigned int numHeights;
	unsigned int maxHeights;
	unsigned int maxJdim;
	unsigned int maxKdim;
	vector<TCcenter> tcVector;
	QDomDocument domDoc;
	QHash<QString, QString> configHash;
	int referenceState;
	
	// Data Processing
	QHash<QString, int> dataSuffix;
	enum dataFormats {
		unknown,
		cen,
		frd,
		cls,
		sec,
		ten,
		swp,
		sfmr,
		wwind,
		qscat,
		ascat,
		nopp,
		eol,
		cimss
	};
	
	enum referenceStates {
		jordan
	};
	
	enum referenceVariables {
		qvbhypref,
		rhoaref,
		rhoref,
		href,
		tempref,
		pressref,
	};
	
	bool read_frd(QFile& metFile, QList<MetObs>* metObVector);
	bool read_cls(QFile& metFile, QList<MetObs>* metObVector);
	bool read_wwind(QFile& metFile, QList<MetObs>* metObVector);
	bool read_eol(QFile& metFile, QList<MetObs>* metObVector);
	bool read_sec(QFile& metFile, QList<MetObs>* metObVector);
	bool read_ten(QFile& metFile, QList<MetObs>*metObVector);
	bool read_dorade(QFile& metFile, QList<MetObs>* metObVector);
	bool read_sfmr(QFile& metFile, QList<MetObs>* metObVector);
	bool read_qscat(QFile& metFile, QList<MetObs>* metObVector);
	bool read_ascat(QFile& metFile, QList<MetObs>* metObVector);
	bool read_nopp(QFile& metFile, QList<MetObs>* metObVector);
	bool read_cimss(QFile& metFile, QList<MetObs>* metObVector);
	bool readTCcenters();
	bool readXMLconfig(const QString& xmlfile);
	real getReferenceVariable(const int& refVariable, const real& heightm, const int& dz = 0);
	real bhypTransform(real qv);
	real bhypInvTransform(real qvbhyp);
	
};

#endif
