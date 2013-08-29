/*
 *  CostFunction.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunction.h"
#include <cmath>
#include <iostream>
#include <limits>

CostFunction::CostFunction(const int& numObs, const int& stateSize)
{
	
	// Create a cost function
	mObs = numObs;
	nState = stateSize;
	
}

CostFunction::~CostFunction()
{
	
}

void CostFunction::setNumObservations(const int& numObs)
{
	mObs = numObs;
}

int CostFunction::getNumObservations()
{
	return mObs;
}

void CostFunction::setLengthStateVector(const int& stateSize)
{
	nState = stateSize;
}

int CostFunction::getLengthStateVector()
{
	return nState;
}

bool CostFunction::minimize()
{
	
	real ftol = 1.0e-5;
	real minimum = 1e34;
	cout << "\tInner Loop Conjugate Gradient" << endl;
	conjugateGradient(currState, currGradient, ftol, minimum);
	return true;
	
}


void CostFunction::conjugateGradient(real* q, real* xi, const real ftol, real fret)
{
	const int ITMAX = 2000;
	const real EPS = 1.0e-18;
	int j, its;
	real gg, gam, fq, dgg;
	
	real* g = new real[nState];
	real* h = new real[nState];

	fq = funcValue(q);
	funcGradient(q, xi);
	for (j=0; j<nState; j++) {
		g[j] = -xi[j];
		xi[j] = h[j] = g[j];
	}
	for (its=0; its<ITMAX; its++) {
		cout << "\t\tIteration: " << its << "\tJ: " << fq << endl;
		dlinmin(q, xi, fret);			
		if (2.0*fabs(fret-fq) <= ftol*(fabs(fret)+fabs(fq)+EPS)) {
			cout << "\tMinimum J: " << fret << endl;
			cout << "\tFound minimum in " << its << " iterations." << endl;
			delete[] g;
			delete[] h;
			return;
		}
		fq = fret;
		funcGradient(q, xi);
		dgg = gg = 0.0;
		for (j=0; j<nState; j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			return;
		}
		gam = dgg/gg;
		for (j=0; j<nState; j++) {
			g[j] = -xi[j];
			xi[j] = h[j] = g[j] + gam*h[j];
		}
	}
	// Got here there were too many iterations
	delete[] g;
	delete[] h;
	cout << "Iterations exceeded in inner minimization loop" << endl;
	return;
}

void CostFunction::dlinmin(real* &p, real* &xi, real &fret)
{
	const real TOL=2.0e-8;
	int j;
	real xx,xmin,fx,fb,fa,bx,ax;
	
	// Fill the temporary state vector
	for (j=0; j<nState; j++) {
		tempState[j] = p[j];
		tempGradient[j] = xi[j];
	}
	ax = 0.0;
	xx = 1.0;
	xmin = 0.0;
	mnbrack(ax,xx,bx,fa,fx,fb);
	fret = dbrent(ax,xx,bx,TOL,xmin);
	for (j=0; j<nState; j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
}

real CostFunction::f1dim(const real x)
{
	int j;
	real f1 = 0.0;
	for (j=0; j<nState; j++) xt[j] = tempState[j] + x*tempGradient[j];
	f1 = funcValue(xt);
	return f1;
}

real CostFunction::df1dim(const real x)
{
	int j;
	real df1 = 0.0;
	for (j=0; j<nState; j++) {
		xt[j] = tempState[j] + x*tempGradient[j];
		df[j] = 0;
	}
	funcGradient(xt, df);
	for (j=0; j<nState; j++) df1 += df[j]*tempGradient[j];
	return df1;
}

inline void CostFunction::mov3(real &a, real &b, real &c,
							   const real d, const real e, const real f)
{
	a=d; b=e; c=f;
}

real CostFunction::dbrent(const real ax, const real bx, const real cx,
							const real tol, real &xmin)
{
	const int ITMAX=100;
	const real ZEPS=numeric_limits<real>::epsilon()*1.0e-3; // Replace with actual machine epsilon
	bool ok1, ok2;
	int iter;
	real a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
	real fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
	
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=f1dim(x);
	dw=dv=dx=df1dim(x);
	for (iter=0; iter<ITMAX; iter++) {
		xm = 0.5*(a+b);
		tol1 = tol*fabs(x)+ZEPS;
		tol2 = 2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			xmin = x;
			return fx;
		}
		if (fabs(e) > tol1) {
			d1 = 2.0*(b-a);
			d2 = d1;
			if (dw != dx) d1 = (w-x)*dx/(dx-dw);
			if (dv != dx) d2 = (v-x)*dx/(dx-dv);
			u1 = x+d1;
			u2 = x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde = e;
			e = d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu=f1dim(u);
		} else {
			u=x+SIGN(tol1,d);
			fu = f1dim(u);
			if (fu > fx) {
				xmin = x;
				return fx;
			}
		}
		du = df1dim(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			mov3(v,fv,dv,w,fw,dw);
			mov3(w,fw,dw,x,fx,dx);
			mov3(x,fx,dx,u,fu,du);
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				mov3(v,fv,dv,w,fw,dw);
				mov3(w,fw,dw,u,fu,du);
			} else if (fu < fv || v == x || v == w) {
				mov3(v,fv,dv,u,fu,du);
			}
		}
	}
	// Too many iterations
	return 0.0;
}

inline void CostFunction::shft3(real &a, real &b, real &c, const real d)
{
	a=b;
	b=c;
	c=d;
}

void CostFunction::mnbrack(real &ax, real &bx, real &cx, 
						   real &fa, real &fb, real &fc)
{
	const real GOLD=1.618034, GLIMIT=100.0, TINY=1.0e-20;
	real ulim,u,r,q,fu;
	
	fa = f1dim(ax);
	fb = f1dim(bx);
	if (fb > fa) {
		SWAP(ax, bx);
		SWAP(fb, fa);
	}
	cx = bx+GOLD*(bx-ax);
	fc = f1dim(cx);
	while (fb > fc) {
		r = (bx - ax)*(fb - fc);
		q = (bx - cx)*(fb - fa);
		u = bx - ((bx-cx)*q-(bx-ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r), TINY), q-r));
		ulim = bx + GLIMIT*(cx-bx);
		if ((bx-u)*(u-cx) > 0.0) {
			fu = f1dim(u);
			if (fu < fc) {
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
				return;
			} else if (fu > fb) {
				cx = u;
				fc = fu;
				return;
			}
			u = cx + GOLD*(cx-bx);
			fu = f1dim(u);
		} else if ((cx-u)*(u-ulim) > 0.0) {
			fu = f1dim(u);
			if (fu < fc) {
				shft3(bx,cx,u,u+GOLD*(u-cx));
				shft3(fb,fc,fu,f1dim(u));
			}
		} else if ((u-ulim)*(ulim-cx) >= 0.0) {
			u = ulim;
			fu = f1dim(u);
		} else {
			u = cx + GOLD*(cx-bx);
			fu = f1dim(u);
		}
		shft3(ax,bx,cx,u);
		shft3(fa,fb,fc,fu);
	}
}


