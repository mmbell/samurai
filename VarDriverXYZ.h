/*
 *  VarDriverXYZ.h
 *  samurai
 *
 *  Created by Michael Bell on 4/12/08.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVERXYZ_H
#define VARDRIVERXYZ_H

#include "VarDriver.h"
#include "BSpline.h"
#include "Observation.h"
#include "CostFunctionXYZ.h"
#include "MetObs.h"
#include "FrameCenter.h"
#include <iostream>
#include <vector>
#include <QHash>
#include <QDir>
#include <QList>
#include <QString>

using namespace std;

class VarDriverXYZ : public VarDriver
{

public:
	VarDriverXYZ();
	~VarDriverXYZ();
	
	// ESMF type calls
	bool initialize(const QDomElement& configuration);
	bool run();
	bool finalize();
	
private:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	
	// Common methods
	void preProcessMetObs();
	bool loadMetObs();
	bool loadBGfromFile();
	int loadBackgroundObs();
	void updateAnalysisParams();

	vector<Observation> obVector;
	int maxIter;

	// Cost Functions
	CostFunctionXYZ* obCostXYZ;
	CostFunctionXYZ* bgCostXYZ;
	
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
};

#endif
