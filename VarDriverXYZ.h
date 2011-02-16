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
#include "TCcenter.h"
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
	bool bilinearMish();
	real bilinearField(real xPos, real yPos, int var);
	
	vector<real> x;
	vector<real> y;
	unsigned int* RnumGridpts;
	vector<Observation> obVector;
	int bc;
	double iincr;
	double jincr;
	double kincr;
	double CQTOL;
	int maxIter;
	real zLevel;
	
	vector<real>** BG;
	vector<real>** BGsave;
	
	// Passable variables
	real* bgB;
	real* bgU;
	real* bgWeights;
	real* obs;
	real* bgObs;
	const real* iaScalar;
	const real* iaVector;
	const real* jaScalar;
	const real* jaVector;
	real* ia;
	real* ja;
	real imin, imax, jmin, jmax, kmin, kmax;
	int idim;
	int jdim;
	int kdim;
	// Cost Functions
	CostFunctionXYZ* obCostXYZ;
	CostFunctionXYZ* bgCostXYZ;
};

#endif
