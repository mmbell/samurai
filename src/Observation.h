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
	
	double inline getWeight(const unsigned int& var) { return weight[var]; }
	void setWeight(const double &wgt, const unsigned int& var);

    double inline getDWeight(const unsigned int& var) { return dweight[var]; }
	void setDWeight(const double &dwgt, const unsigned int& var);

    double inline getD2Weight(const unsigned int& var) { return d2weight[var]; }
	void setD2Weight(const double &d2wgt, const unsigned int& var);

	// Note that you set this as error, but it returns 1/error or error
	double inline getInverseError() const { return (1./error); }
	double inline getError() const { return error; }
	void setError(const double &err);
	
	double inline getOb() const { return obNet; }
	void setOb(const double &ob);

	int inline getTime() const { return time; }
	void setTime(const int& t);
	
private:
	double radius;
	double theta;
	double altitude;
	double cartesianX;
	double cartesianY;
	double weight[7];
    double dweight[7];
    double d2weight[7];
	double error;
	double obNet;
	unsigned int numVars;
	int type;
	int time;
};

#endif
