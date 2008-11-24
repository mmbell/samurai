/*
 *  RecursiveFilter.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "RecursiveFilter.h"
#include <iostream>
#include <cmath>

RecursiveFilter::RecursiveFilter(const int& fOrder, const int& fLengthScale)
{
	
	// Create a filter based on the specified order and lengthscale
	order = fOrder;
	lengthScale = fLengthScale;
	getFilterCoefficients();
	
}

RecursiveFilter::~RecursiveFilter()
{
}

void RecursiveFilter::getFilterCoefficients()
{
	double sigma = (lengthScale*lengthScale)/2;
	double b[5][5];
	double a[201][201];
	double p[201];
	double K[9];
	double D[9];
	
	// b coefficients
	b[1][1] = 1;
	b[1][2] = 1/12;
	b[1][3] = 1/90;
	b[1][4] = 1/560;
	b[2][2] = 1;
	b[2][3] = 1/6;
	b[2][4] = 7/240;
	b[3][3] = 1;
	b[3][4] = 1/4;
	b[4][4] = 1;
	
    // Compute polynomial coefficients
	K[2] = b[1][1]*sigma;
	K[4] = b[1][2]*sigma + b[2][2]*(pow(sigma,2))/factorial(2);
	K[6] = b[1][3]*sigma + b[2][3]*(pow(sigma,2))/factorial(2)
		+ b[3][3]*(pow(sigma,3))/factorial(3);
	K[8] = b[1][4]*sigma + b[2][4]*(pow(sigma,2))/factorial(2)
		+ b[3][4]*(pow(sigma,3))/factorial(3) + b[4][4]*(pow(sigma,4))/factorial(4);
	
	for (int i=0;i<=order*2;i++) {
		D[i] = 0;
	}
	
	for (int n=2;n<=order*2;n+=2) {
		for (int i=0;i<=n;i++) {
			int index = n/2 - i;
			double coeff = (pow((-1),(n/2)))*(pow((-1),i))*factorial(n)/(factorial(n-i)*factorial(i));
			D[index+order] += K[n]*coeff;
		}
	}
	D[order] += 1;
	
	for (int i=0;i<=200;i++) {
		for (int j=0;j<=200;j++) {
			if (abs(j-i) <= order) {
				a[i][j] = D[j-i+order];
			} else {
				a[i][j] = 0;
			}
		}
	}
	
	// Cholesky decomp
	for (int i=0;i<=200;i++) {
		for (int j=0;j<=200;j++) {
			double sum=a[i][j];
			for (int k=i-1;k>=0;k--) {
				sum -= a[i][k]*a[j][k];
			}
			if (i == j) {
				if (sum <= 0.0) { 
					std::cout << "cholesky failed at i,j sum\n";
					break;
				} else {
					p[i] = sqrt(sum);
				}
			} else {
				a[j][i]=sum/p[i];
			}
		}
	}
	beta = 1/p[200];
	alpha[1]=-beta*a[200][199];
	if (order > 1) {
		alpha[2]=-beta*a[200][198];
	} else {
		alpha[2]=0;
	}
	if (order > 2) {
		alpha[3]=-beta*a[200][197];
	} else {
		alpha[3]=0;
	}
	if (order > 3) {
		alpha[4]=-beta*a[200][196];
	} else {
		alpha[4]=0;
	}
	
	std::cout << "4th Order Recursive Filter coefficients computed for lengthscale " 
		<< lengthScale << " Delta R" << std::endl;
	std::cout << "\tBeta: " << beta << std::endl;
	std::cout << "\tAlpha[1]: " << alpha[1] << std::endl;
	std::cout << "\tAlpha[2]: " << alpha[2] << std::endl;
	std::cout << "\tAlpha[3]: " << alpha[3] << std::endl;
	std::cout << "\tAlpha[4]: " << alpha[4] << std::endl;
		
	
	// Boundary conditions
	double L[5][5];
	double U[5][5];
	double LT[5][5];
	double UT[5][5];
	double LI[5][5];
	double col[5];
	double colinv[5];
	
	for (int i=1;i<=order;i++) {
		for (int j=1;j<=order;j++) {
			L[i][j]=0.;
			U[i][j]=0.;
		}
	}
	
	for (int i=1;i<=order;i++) {
		for (int j=0;j<=order;j++) {
			int index = i+j;
			if (index <= order) {
				L[i+j][i] = -alpha[j];
				U[i][i+j] = alpha[order-j];
			}
		}
		L[i][i] = 1.;
	}
	
	// Invert L
	for (int j=1;j<=order;j++) {
		for (int i=1;i<=order;i++) { col[i] = 0.0; }
		col[j] = 1;
		colinv[0] = col[0];
		for (int i=1;i<=order;i++) {
			colinv[i] = col[i];
			for (int k=i-1;k>=1;k--) {
				colinv[i] -= L[i][k]*colinv[k];
			}
		}
		for (int i=1;i<=order;i++) {
			LI[i][j] = colinv[i];
		}
	}
	
	// Transpose L & T
	for (int i=1;i<=order;i++) {
		for (int j=1;j<=order;j++) {
			LT[i][j] = L[j][i];
			UT[i][j] = U[j][i];
		}
	}
	
	// Get Sn matrix
	double tmp[5][5];
	double tmp2[5][5];
	for (int i=1;i<=order;i++) {
		for (int j=1;j<=order;j++) {
			tmp[i][j] = 0.;
			for (int k=1;k<=order;k++) {
				tmp[i][j] += UT[i][k]*LI[k][j];
			}
		}
	}
	for (int i=1;i<=order;i++) {
		for (int j=1;j<=order;j++) {
			tmp2[i][j] = 0.;
			for (int k=1;k<=order;k++) {
				tmp2[i][j] += tmp[i][k]*U[k][j];
			}
		}
	}
	
	for (int i=1;i<=order;i++) {
		for (int j=1;j<=order;j++) {
			Sn[i-1][j-1] = (LT[i][j] - tmp2[i][j])/beta;
		}
	}
	
}

bool RecursiveFilter::filterArray(double* array, const int& arrLength)
{
	int maxi = arrLength-1;
	double* p = array;
	double* q = new double[arrLength];
	double* s = new double[arrLength];

	for (int i=0; i<= maxi; i++) {
		q[i] = 0;
		s[i] = 0;
	}
	
	q[0]=beta*p[0];
	q[1]=beta*p[1] + alpha[1]*q[0];
	q[2]=beta*p[2] + alpha[1]*q[1] + alpha[2]*q[0];
	q[3]=beta*p[3] + alpha[1]*q[2] + alpha[2]*q[1] + alpha[3]*q[0];	
	for (int i=order; i<= maxi; i++) {
		q[i] = beta*p[i] + alpha[1]*q[i-1]
		+ alpha[2]*q[i-2] + alpha[3]*q[i-3] + alpha[4]*q[i-4];
	}
	
    // Invert Sn
	double* A = new double[4];
	double* B = new double[4];
	for (int i=maxi-order+1; i<= maxi; i++) {
		B[i-(maxi-order+1)] = q[i];
	}
	solveBC(A, B);
	for (int i=maxi; i>= (maxi-order+1); i--) {
		s[i] = A[i-(maxi-order+1)];
	}
	for (int i=maxi-order;i>=0;i--) {
		s[i] = beta*q[i] + alpha[1]*s[i+1]
		+ alpha[2]*s[i+2] + alpha[3]*s[i+3] + alpha[4]*s[i+4];
		// std::cout << s[i] << std::endl;
	}
	delete[] A;
	delete[] B;
	
	
	for (int i=0; i<= maxi; i++) {
		// To get a 'true' Gaussian, need to scale by this factor
		// To preserve total quantity (needed for Variational analysis)
		// Do not scale resulting vector
		//double Pi = 3.141592653589793238462643;
		//array[i] = s[i]*sqrt(2*Pi)*lengthScale;
		array[i] = s[i];
	}
	delete[] q;
	delete[] s;
	
	return true;
	
}

double RecursiveFilter::factorial(const double& max) 
{
	double n = 1;
	for (double i=2;i<=max;i++) {
		n *= i;
	}
	return n;
	
}

void RecursiveFilter::solveBC(double* A, double*B) 
{
	
	int n = order-1;
	for(int j=0;j<=n-1;j++) {
		for (int i=j+1;i<=n;i++) {
			Sn[i][j]=Sn[i][j]/Sn[j][j];
		}
		for	(int i=j+1;i<=n;i++) {
			for (int k=j+1;k<=n;k++) {
				Sn[i][k]=Sn[i][k]-Sn[i][j]*Sn[j][k];
			}
			B[i]=B[i]-Sn[i][j]*B[j];
		}
	}
	A[n]=B[n]/Sn[n][n];
	for(int j=n-1;j>=0;j--) {
		A[j]=B[j];
		for(int k=n;k>=j+1;k--) {
			A[j]=A[j]-A[k]*Sn[j][k];
		}
		A[j]=A[j]/Sn[j][j];
	}
}
