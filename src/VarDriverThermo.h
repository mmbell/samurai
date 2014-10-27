/*
 *  VarDriverThermo.h
 *  samurai
 *
 *  Created by Annette Foerster on 8/7/13.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVERTHERMO_H
#define VARDRIVERTHERMO_H

#include "VarDriver3D.h"
#include "VarDriver.h"
#include "BSpline.h"
#include "Observation.h"
#include "CostFunctionThermoXYZ.h"
#include "CostFunctionThermoRTZ.h"
#include "CostFunctionThermo.h"
#include "MetObs.h"
#include "FrameCenter.h"
#include "NetCDF.h"
#include "NetCDF_RTZ.h"
#include "NetCDF_XYZ.h"
#include <iostream>
#include <vector>
#include <QHash>
#include <QDir>
#include <QList>
#include <QString>

using namespace std;

class VarDriverThermo : public VarDriver3D
{

public:
	VarDriverThermo();
	~VarDriverThermo();
	
	// ESMF type calls
	bool initialize(const QDomElement& configuration);
	bool run();
	bool finalize();
  bool testing(QList<Observation>* obVector);
  bool testing_rtz(QList<Observation>* obVector);

	
private:
   int test;
   CostFunctionThermo* obCost3D;
   bool loadObservations(QString& metFile, QList<Observation>* obVector);
   NetCDF* ncFile;
   QList<Observation> obVector;
   
   real* obs;
};

#endif
