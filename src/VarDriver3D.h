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
#include "Xml.h"
#include <iostream>
#include <vector>

using namespace std;

class VarDriver3D : public VarDriver
{

public:
	VarDriver3D();
	virtual ~VarDriver3D();

	// ESMF type calls
	bool initialize(const XMLNode& configuration);
	bool initialize(const samurai_config &configSam);
	bool initialize();

	void dumpBgu();
	void dumpBgIn();

	bool run();
	bool run(int nx, int ny, int nsigma,
		 // ----- new -----
		 char cdtg[10],	// "12Z oct 4 2015 -> "2015100412"
		 int delta,	// delta * iter1 past cdtg
		 int iter1,
		 float imin, float imax, float iincr, // used to come from config
		 float jmin, float jmax, float jincr,
		 // ----- new -----

		 float *sigmas,
		 float *latitude,	// 2D arrays
		 float *longitude,
		 float *terrain_Height,
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
	// A fixed grid is defined in the configuration.
	// A non-fixed grid is given as part of the run(... nx, ny, nz, ...) call

	void setGridFlag(bool flag) { fixedGrid = flag; };
	bool isFixedGrid() { return fixedGrid == true; };

private:

	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;   // This is also declared in ReferenceState.h

	// Common methods
	bool validateDriver();
	bool validateConfig();	// was validateXMLconfig();

	bool validateFixedGrid();					// grid comes from config
	bool validateRunGrid(float nx, float ny, float nz,		// grid comes from run call
			     float imin, float imax, float iincr,
			     float jmin, float jmax, float jincr,
			     float *sigmas);
	bool validateFractlGrid();				// grid comes from Fractl netcdf file
	bool validateGrid();						// called by both validate*Grid above

	void fillRunCenters(char cdtg[10], int delta, int iter, float lat, float lon);
	bool initObCost3D();
	bool gridDependentInit();
	bool preProcessMetObs();
	bool loadMetObs();
	bool loadPreProcessMetObs();
	bool loadBGfromFile();
	bool loadBackgroundCoeffs();
	int loadBackgroundObs(const char *background_fname);
	int loadBackgroundObs(int nx, int ny, int nsigma,
			      char *ctdg, int delta, int iter, // time elements
			      float *sigmas,
			      float *latitude,
			      float *longitude,
						float *terrain_Height,
			      float *u1,
			      float *v1,
			      float *w1,
			      float *th1,
			      float *p1);
	int loadBackgroundObs();
	bool adjustBackground();
	// bool adjustBackground(const int& bStateSize);
	void updateAnalysisParams(const int& iteration);

	bool findReferenceCenter();

	std::vector<real> bgIn;
	std::vector<Observation> obVector;
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

	// These 2 control sizes of data structure. Pulling them here
	// since they were computed the same way in 2 different places.

	int64_t bStateSize;
	int64_t uStateSize;

	BkgdAdapter *bkgdAdapter;

	// Some cost functions (COAMPS) might need access to the sigmas
	float *sigmaTable;

    enum RunModes {
        XYZ,
        RTZ
    };
};

#endif
