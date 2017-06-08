/*
 *  VarDriver3D.h
 *  samurai
 *
 *  Created by Michael Bell on 4/12/08.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVER3D_H
#define VARDRIVER3D_H

#include "VarDriver.h"
#include "BSpline.h"
#include "Observation.h"
#include "CostFunctionXYZ.h"
#include "CostFunctionXYP.h"
#include "CostFunctionRTZ.h"
#include "CostFunctionCOAMPS.h"
#include "MetObs.h"
#include "FrameCenter.h"
#include "BkgdAdapter.h"
#include <iostream>
#include <vector>
#include <QHash>
#include <QDir>
#include <QList>
#include <QString>

using namespace std;

class VarDriver3D : public VarDriver
{

public:
	VarDriver3D();
	~VarDriver3D();

	// ESMF type calls
	bool initialize(const QDomElement& configuration);
	bool initialize(const samurai_config &configSam);
	
	bool run();
	bool run(int nx, int ny, int nsigma,
		 float dx, float dy,
		 float *sigmas,
		 float *latitude,	// 2D arrays
		 float *longitude,
		 float *u1,		// 3D arrays (nx, ny, nsigma)  
		 float *v1,
		 float *w1,
		 float *th1,
		 float *p1,
		 // These are output values
		 float *usam,		// 3D array
		 float *vsam,
		 float *wsam,
		 float *thsam,
		 float *psam);
	bool finalize();

private:

	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;   // This is also declared in ReferenceState.h
	
	// Common methods
	bool validateDriver();
	bool preProcessMetObs();
	bool loadMetObs();
	bool loadBGfromFile();
	
	int loadBackgroundObs(const char *background_fname);
	int loadBackgroundObs(int nx, int ny, int nsigma,
			      float dx, float dy,
			      float *sigmas,
			      float *latitude,
			      float *longitude,
			      float *u1,
			      float *v1,
			      float *w1,
			      float *th1,
			      float *p1);
	int loadBackgroundObs();
	
	int loadBackgroundCoeffs();
	bool adjustBackground(const int& bStateSize);
	void updateAnalysisParams(const int& iteration);
	bool validateXMLconfig();
	
	QList<real> bgIn;
	QList<Observation> obVector;
	int maxIter;

	// Cost Functions
	CostFunction3D* obCost3D;
	CostFunction3D* bgCost3D;

	// Variables passed to Cost function
	real* bgB;
	real* bgU;
	real* bgWeights;
	real* obs;
	real* bgObs;
	real imin, imax, jmin, jmax, kmin, kmax;
	real iincr;
	real jincr;
	real kincr;
	int idim;
	int jdim;
	int kdim;
	int runMode;

	BkgdAdapter *bkgdAdapter;
	
    enum RunModes {
        XYZ,
        RTZ
    };
};

#endif
