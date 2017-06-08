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
#include "ReferenceState.h"
#include <iostream>
#include <vector>
#include <QHash>
#include <QDir>
#include <QList>
#include <QString>
#include <QDomDocument>
#include "precision.h"
#include "Projection.h"
#include "samurai.h"

using namespace std;

class VarDriver
{

public:

	VarDriver();
	virtual ~VarDriver();
	virtual bool initialize(const QDomElement& configuration) = 0;
	virtual bool initialize(const samurai_config &configSam) = 0;

	virtual bool run() = 0;
	virtual bool run(int nx, int ny, int nsigma,
			 float dx, float dy,
			 float *sigmas,
			 float *latitude, // 2D arrays
			 float *longitude,
			 float *u1,	// 3D array (nx, ny, nsigma)  
			 float *v1,
			 float *w1,
			 float *th1,
			 float *p1,
			 // These are output values
			 float *usam,	// 3D array
			 float *vsam,
			 float *wsam,
			 float *thsam,
			 float *psam) = 0;

	virtual bool finalize() = 0;

 protected:
	
	real CoriolisF;
	real Pi;
	unsigned int numVars;
	unsigned int numHeights;
	unsigned int maxHeights;
	unsigned int numDerivatives;
	unsigned int maxJdim;
	unsigned int maxKdim;
	vector<FrameCenter> frameVector;
	QDomDocument domDoc;
	QHash<QString, QString> configHash;
	ReferenceState* refstate;
	// Data Processing
	QHash<QString, int> dataSuffix;
	QDir dataPath;
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
		dwl,
		insitu,
		mtp,
		mesonet,
		classnc,
		qcf,
		aeri,
		rad
	};
	Projection projection;
	
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
	bool read_insitu(QFile& metFile, QList<MetObs>* metObVector);
	bool read_mtp(QFile& metFile, QList<MetObs>* metObVector);
	bool read_mesonet(QFile& metFile, QList<MetObs>* metObVector);
	bool read_classnc(QFile& metFile, QList<MetObs>* metObVector);
	bool read_qcf(QFile& metFile, QList<MetObs>* metObVector);
	bool read_aeri(QFile& metFile, QList<MetObs>* metObVector);
	bool read_rad(QFile& metFile, QList<MetObs>* metObVector);
	bool readFrameCenters();
	bool parseXMLconfig(const QDomElement& config);
	bool parseSamuraiConfig(const samurai_config &config);
	Projection::ProjectionType projectionFromConfig();
};

#endif
