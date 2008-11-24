/*
 *  VarDriver1d.h
 *  tcvar
 *
 *  Created by Michael Bell on 4/12/08.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVER1D_H
#define VARDRIVER1D_H

#include "BSpline.h"
#include "Observation.h"
#include "CostFunctionR.h"
#include "CostFunctionAnalytic.h"
#include "VarDriver.h"
#include "MetObs.h"
#include "TCcenter.h"
#include <iostream>
#include <vector>
#include <QHash>
#include <QDir>
#include <QList>
#include <QString>

using namespace std;

class VarDriver1d : public VarDriver
{

public:
	VarDriver1d();
	~VarDriver1d();
	bool run();

private:
	typedef BSplineBase<double> SplineBase;
	typedef BSpline<double> SplineD;

	// Different cases
	bool initParametricBG();
	void runParametricBG();
	bool finalizeParametricBG();
	
	bool initParaboloid();
	void runParaboloid();
	bool finalizeParaboloid();
	
	// Common methods
	void EvalSpline (SplineD* spline, ostream &out);
	void bgOut();
	void processMetObs();
	double updateXforms();
	
	SplineD* bgSpline;
	SplineD* ctrlSpline;
	SplineD* RXformSpline;
	vector<double> r;
	vector<double> R;
	vector<double> initCtrl;
	vector<double> RXform;
	vector<double> rXform;
	vector<Observation> obVector;
	int bc;
		
	vector<double>* BG;
	vector<double>* BGsave;
		
	// Cost Functions
	CostFunctionR* cost1d;
	CostFunctionAnalytic costAnalytic;
};

#endif
