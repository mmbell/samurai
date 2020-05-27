/*
 *  RecursiveFilter.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef RECURSIVEF_H
#define RECURSIVEF_H
#include "precision.h"


class RecursiveFilter
{
	
public:
	RecursiveFilter(const int& fOrder);
	RecursiveFilter(const int& fOrder, const double& fLengthScale);
	RecursiveFilter(const double& fLengthScale, const double* tau, const int& arrLength );
	~RecursiveFilter();
    void setFilterLengthScale(const double& fLengthScale);
	#pragma acc routine seq 
	bool filterArray(double* array, double *q, double *s, const int& arrLength);
	bool filterArray(double* array, const int& arrLength);
	bool aniFilterArray(double* array, const int& arrLength);

private:
	int order;
	double lengthScale;
	double beta;
	double alpha[5];
	double* abeta;
	double* aalpha[5];
	double Sn[5][5];
	void getIsotropicFilterCoefficients();
	void getAnisotropicFilterCoefficients(const double* tau, const int& arrLength);
	double factorial(const double& max); 		
  #pragma acc routine seq
	void solveBC(double* A, double* B, double S[5][5]);	
};

#endif
