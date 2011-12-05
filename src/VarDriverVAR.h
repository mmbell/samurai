/*
 *  VarDriverVAR.h
 *  samurai
 *
 *  Created by Michael Bell on 4/12/08.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVERVAR_H
#define VARDRIVERVAR_H

#include "VarDriver.h"
#include "BSpline.h"
#include "Observation.h"
#include "CostFunctionVAR.h"
#include "MetObs.h"
#include "FrameCenter.h"
#include <iostream>
#include <vector>
#include <QHash>
#include <QDir>
#include <QList>
#include <QString>

using namespace std;

class VarDriverVAR : public VarDriver
{

public:
	VarDriverVAR();
	~VarDriverVAR();
	// ESMF type calls
	bool initialize();
	bool run();
	bool finalize();
	
private:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	
	// Common methods
	void preProcessMetObs();
	bool loadMetObs();
	bool loadBGfromFile();
	bool bilinearMish();
	real bilinearField(real xPos, real yPos, int var);
	bool setupMishAndRXform();
	
	vector<real> x;
	vector<real> y;
	unsigned int* RnumGridpts;
	vector<Observation> obVector;
	int bc;
	double yincr;
	double xincr;
	double CQTOL;
	int maxIter;
	real zLevel;
	
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
	CostFunctionVAR* costVAR;

};

#endif
