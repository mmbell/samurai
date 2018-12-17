/*
 *  CostFunction.h
 *  SAMURAI
 *
 *  Copyright 2010 Michael Bell. All rights reserved.
 *
 */

#ifndef COSTFUNC_H
#define COSTFUNC_H
#include "precision.h"
#include "Projection.h"

using namespace std;

class CostFunction
{
	
public:
	
  CostFunction(const Projection& proj, const int& numObs = 0, const int& stateSize = 0);
	virtual ~CostFunction();
	void setNumObservations(const int& numObs);
	int getNumObservations();
	void setLengthStateVector(const int& stateSize);
	int getLengthStateVector();
	bool minimize();
	
protected:
	int mObs;
	int nState;
	real* currState;
	real* currGradient;
	real* tempState;
	real* tempGradient;
	real* xt;
	real* df;
	const Projection& projection;
	       
	       
	virtual real funcValue(real* state) = 0;
	virtual void funcGradient(real* state, real* gradient) = 0;
	
	void conjugateGradient(real* q, real* xi, const real ftol, real funcMin);
	void dlinmin(real* &p, real* &xi, real &fret);
	real f1dim(const real x);
	real df1dim(const real x);
	inline void mov3(real &a, real &b, real &c,
					 const real d, const real e, const real f);
	real dbrent(const real ax, const real bx, const real cx,
				  const real tol, real &xmin);
	inline void shft3(real &a, real &b, real &c, const real d);
	void mnbrack(real &ax, real &bx, real &cx, 
				 real &fa, real &fb, real &fc);
	
private:
	inline real MAX(const real &a, const real &b) 
	{return b > a ? (b) : (a); }
	inline real SIGN(const real &a, const real &b) 
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }
	inline void SWAP(real &a, real &b)
	{real dum=a; a=b; b=dum;}
	
};

#endif

