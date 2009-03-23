/*
 *  Observation.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "Observation.h"

Observation::Observation()
{
	numVars = 6;
}

Observation::~Observation()
{
}

void Observation::setRadius(const double &r)
{
	radius = r;
}

void Observation::setAltitude(const double &alt)
{
	altitude = alt;
}

void Observation::setType(const int &t)
{
	type = t;
}

void Observation::setWeight(const double &wgt, const unsigned int& var)
{
	if (var < numVars) {
		weight[var] = wgt;
	}
}

void Observation::setError(const double &err)
{
	if (err != 0) {
		error = 1./err;
	} else {
		error = 1.;
	}
}


void Observation::setOb(const double &ob)
{
	obNet = ob;
}
