/*
 *  VarDriverRZ.h
 *  tcvar
 *
 *  Created by Michael Bell on 4/12/08.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVERRZ_H
#define VARDRIVERRZ_H

#include "VarDriver.h"
#include "BSpline.h"
#include "Observation.h"
#include "CostFunctionRZ_CPU.h"
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

class VarDriverRZ : public VarDriver
{

public:
	VarDriverRZ();
	~VarDriverRZ();
	// ESMF type calls
	bool initialize();
	bool run();
	bool finalize();
	
private:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	
	// Common methods
	void EvalSpline (SplineD* spline, ostream &out);
	void bgOut();
	void preProcessMetObs();
	bool loadMetObs();
	bool loadBGfromFile();
	bool bilinearMish();
	real bilinearField(real radius, real height, int var);
	bool setupMishAndRXform();
	double updateXforms();
	bool writeAsi(const QString& fileName, vector<real>** fields);
	
	SplineD* scalarSpline;
	SplineD* vecSpline;
	SplineD* zSpline;
	SplineD* zSplinePsi;
	SplineD* ctrlSpline;
	SplineD* rXformSpline;
	vector<real> r;
	vector<real> R;
	vector<real> z;
	vector<real> initCtrl;
	vector<real>* RXform;
	vector<real>* rXform;
	unsigned int* RnumGridpts;
	vector<Observation> obVector;
	int bc;
	double zincr;
	double rincr;
	double CQTOL;
	int maxIter;
	
	vector<real>** BG;
	vector<real>** BGsave;
	
	// Passable variables
	real* bgB;
	real* bgU;
	real* obs;
	const real* iaScalar;
	const real* iaVector;
	const real* jaScalar;
	const real* jaVector;
	real* ia;
	real* ja;
	real imin, imax, jmin, jmax;
	int idim;
	int jdim;
	
	// Cost Functions
	CostFunctionRZ_CPU* costRZ;

};

#endif
