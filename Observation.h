/*
 *  Observation.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef OBSERVATION_H
#define OBSERVATION_H
#include "precision.h"


class Observation
{

public:
	Observation();
	~Observation();
	double inline getRadius() const { return radius; }
	void setRadius(const double &r);
	
	double inline getAltitude() const { return altitude; }
	void setAltitude(const double &alt);
	
	void setType(const int &t);
	int inline getType() const { return type; }
	
	double inline getWeight(const unsigned int& var) { return weight[var]; }
	void setWeight(const double &wgt, const unsigned int& var);
	
	// Note that you set this as error, but it returns 1/error or error
	double inline getInverseError() const { return (1./error); }
	double inline getError() const { return error; }
	void setError(const double &err);
	
	double inline getOb() const { return obNet; }
	void setOb(const double &ob);

private:
	double radius;
	double altitude;
	double weight[6];
	double error;
	double obNet;
	unsigned int numVars;
	int type;
};

#endif
