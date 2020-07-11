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
	
	double inline getTheta() const { return theta; }
	void setTheta(const double &t);
	
	double inline getAltitude() const { return altitude; }
	void setAltitude(const double &alt);
	
	double inline getCartesianX() const { return cartesianX; }
	void setCartesianX(const double &x);
	
	double inline getCartesianY() const { return cartesianY; }
	void setCartesianY(const double &y);
	
	void setType(const int &t);
	int inline getType() const { return type; }
	
	double getWeight(const unsigned int& var, const unsigned int& d);
	void setWeight(const double &wgt, const unsigned int& var, const unsigned int& derivative = 0);

	// Note that you set this as error, but it returns 1/error or error
    double getInverseError() const;
	double inline getError() const { return error; }
	void setError(const double &err);
	
	double inline getOb() const { return obNet; }
	void setOb(const double &ob);

	int64_t inline getTime() const { return time; }
	void setTime(const int64_t& t);
	
private:
	double radius;
	double theta;
	double altitude;
	double cartesianX;
	double cartesianY;
	double weight[7][4];
	double error;
	double obNet;
	unsigned int numVars;
	int type;
	int64_t time;
};

#endif
