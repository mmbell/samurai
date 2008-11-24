/*
 *  ParametricVortex.h
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef PARMVORTEX_H
#define PARMVORTEX_H

class ParametricVortex
{

  public:
	ParametricVortex(const double* Vin, const double* Radin,  const int& numradii);
	~ParametricVortex();
	float getWindAtRadiusKM(const float& radius);
	float getX1();
	float getn();
	float getA();
	
  private:
	float** vertex;
	float* VT;
	float* vertexSum;
	const double* winds;
	const double* radii;
	float Vg;
	float vmax, rmw, maxr;
	float X1fit, Afit, nfit, R1fit;
	float X2, deltar;
	float Lc, Zetac;
	
	void fitParameters();
	inline void getVertexSum();
	float simplexTest(int& low, float factor);
	float getS2(const float& X1, const float& A, const float& n);

};

#endif
