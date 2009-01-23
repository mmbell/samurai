/*
 *  RecursiveFilter.h
 *  tcvar
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
	RecursiveFilter(const int& fOrder, const int& fLengthScale);
	~RecursiveFilter();
	bool filterArray(double* array, const int& arrLength);
	
private:
	int order;
	int lengthScale;
	double beta;
	double alpha[5];
	double Sn[5][5];
	void getFilterCoefficients();
	double factorial(const double& max); 		
	void solveBC(double* A, double* B);	
};

#endif
