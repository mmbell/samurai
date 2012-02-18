/*
 *  Observation.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "Observation.h"

Observation::Observation()
{
	numVars = 7;
	radius = -999;
	theta = -999;
	altitude = -999;
	cartesianX = -999;
	cartesianY = -999;
	for (unsigned int i = 0; i < numVars; i++) {
        weight[i] = -999;
        dweight[i] = -999;
        d2weight[i] = -999;
    }
	error = -999;
	obNet = -999;
	type = -999;
	time = -999;
}

Observation::~Observation()
{
}

void Observation::setRadius(const double &r)
{
	radius = r;
}

void Observation::setTheta(const double &t)
{
	theta = t;
}

void Observation::setAltitude(const double &alt)
{
	altitude = alt;
}

void Observation::setCartesianX(const double &x)
{
	cartesianX = x;
}

void Observation::setCartesianY(const double &y)
{
	cartesianY = y;
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

void Observation::setDWeight(const double &dwgt, const unsigned int& var)
{
	if (var < numVars) {
		dweight[var] = dwgt;
	}
}

void Observation::setD2Weight(const double &d2wgt, const unsigned int& var)
{
	if (var < numVars) {
		d2weight[var] = d2wgt;
	}
}

void Observation::setError(const double &err)
{
	error = err;
}


void Observation::setOb(const double &ob)
{
	obNet = ob;
}

void Observation::setTime(const int &t)
{
	time = t;
}
