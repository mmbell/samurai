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
#include "MetObs.h"
#include "FrameCenter.h"
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
	virtual bool initialize(const QDomElement& configuration) = 0;
	virtual bool run() = 0;
	virtual bool finalize() = 0;

protected:
	
	real CoriolisF;
	real Pi;
	real rhoBase;
	real rhoInvScaleHeight;
	unsigned int numVars;
	unsigned int numHeights;
	unsigned int maxHeights;
	unsigned int maxJdim;
	unsigned int maxKdim;
	vector<FrameCenter> frameVector;
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
		cimss,
		dwl
	};
	
	enum referenceStates {
		dunion_mt
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
	bool read_dwl(QFile& metFile, QList<MetObs>* metObVector);
	bool readFrameCenters();
	bool parseXMLconfig(const QDomElement& config);
	real getReferenceVariable(const int& refVariable, const real& heightm, const int& dz = 0);
	real bhypTransform(const real& qv);
	real bhypInvTransform(const real& qvbhyp);
	
};

#endif
