/*
 *  VarDriver2d.h
 *  tcvar
 *
 *  Created by Michael Bell on 4/12/08.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVER2D_H
#define VARDRIVER2D_H

#include "VarDriver.h"
#include "BSpline.h"
#include "Observation.h"
#include "CostFunctionRZ.h"
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

class VarDriver2d : public VarDriver
{

public:
	VarDriver2d();
	~VarDriver2d();
	bool run();

private:
	typedef BSplineBase<double> SplineBase;
	typedef BSpline<double> SplineD;

	// Different cases
	bool initVortexBG();
	void runVortexBG();
	bool finalizeVortexBG();
	
	bool initParaboloid();
	void runParaboloid();
	bool finalizeParaboloid();
	
	// Common methods
	void EvalSpline (SplineD* spline, ostream &out);
	void bgOut();
	void processMetObs();
	double updateXforms();
	bool writeAsi(const QString& fileName, vector<double>** fields, SplineD* scalar, SplineD* vector);	
	
	SplineD* scalarSpline;
	SplineD* vecSpline;
	SplineD* zSpline;
	SplineD* zSplinePsi;
	SplineD* ctrlSpline;
	SplineD* rXformSpline;
	vector<double> r;
	vector<double> R;
	vector<double> z;
	vector<double> initCtrl;
	vector<double>* RXform;
	vector<double>* rXform;
	unsigned int* RnumGridpts;
	vector<Observation> obVector;
	int bc;
	
	vector<double>** BG;
	vector<double>** BGsave;
		
	// Cost Functions
	CostFunctionRZ* cost2d;
	CostFunctionAnalytic costAnalytic;
};

#endif
