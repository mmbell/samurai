/*
 *  CostFunction.h
 *  barwind
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNC_H
#define COSTFUNC_H

using namespace std;

class CostFunction
{
	
public:
	
	CostFunction(const int& numObs = 0, const int& stateSize = 0);
	virtual ~CostFunction();
	void setNumObservations(const int& numObs);
	int getNumObservations();
	void setLengthStateVector(const int& stateSize);
	int getLengthStateVector();
	bool minimize();
	
protected:
	int mObs;
	int nState;
	double* currState;
	double* currGradient;
	double* tempState;
	double* tempGradient;
	double* xt;
	double* df;
	virtual double funcValue(double* state) = 0;
	virtual void funcGradient(double* state, double* gradient) = 0;
	
	void conjugateGradient(double* q, double* xi, const double ftol, double funcMin);
	void dlinmin(double* &p, double* &xi, double &fret);
	double f1dim(const double x);
	double df1dim(const double x);
	inline void mov3(double &a, double &b, double &c,
					 const double d, const double e, const double f);
	double dbrent(const double ax, const double bx, const double cx,
				  const double tol, double &xmin);
	inline void shft3(double &a, double &b, double &c, const double d);
	void mnbrack(double &ax, double &bx, double &cx, 
				 double &fa, double &fb, double &fc);
	
private:
	inline double MAX(const double &a, const double &b) 
	{return b > a ? (b) : (a); }
	inline double SIGN(const double &a, const double &b) 
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }
	inline void SWAP(double &a, double &b)
	{double dum=a; a=b; b=dum;}
	
};

#endif

