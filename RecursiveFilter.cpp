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

RecursiveFilter::RecursiveFilter(const int& fOrder, const double& fLengthScale)
{
	
	// Create a filter based on the specified order and lengthscale
	order = fOrder;
	lengthScale = fLengthScale;
	getIsotropicFilterCoefficients();
	
}

RecursiveFilter::RecursiveFilter(const double& fLengthScale, const double* tau, const int& arrLength)
{
	
	// Create a 4th order filter based on the specified lengthscale and tau mapfactors
	order = 4;
	lengthScale = fLengthScale;
	getIsotropicFilterCoefficients();
	getAnisotropicFilterCoefficients(tau, arrLength);
	
}

RecursiveFilter::~RecursiveFilter()
{
}

void RecursiveFilter::getIsotropicFilterCoefficients()
{
	double sigma = (lengthScale*lengthScale)/2;
	double b[5][5];
	double a[201][201];
	double p[201];
	double K[9];
	double D[9];
	
	// b coefficients
	b[1][1] = 1;
	b[1][2] = 1./12.;
	b[1][3] = 1./90.;
	b[1][4] = 1./560.;
	b[2][2] = 1.;
	b[2][3] = 1./6.;
	b[2][4] = 7./240.;
	b[3][3] = 1.;
	b[3][4] = 1./4.;
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
			//std::cout << "Coeff: " << n << "\t" << K[n] << "\t" << coeff << "\n";
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
			//std::cout << a[i][j] << " ";
		}
		//std::cout << std::endl;
	}
	
	// Cholesky decomp
	for (int i=0;i<=200;i++) {
		for (int j=i;j<=200;j++) {
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
		<< lengthScale << " Delta X" << std::endl;
	std::cout << "\tCoefficients: " << beta << "\t";
	std::cout << alpha[1] << "\t";
	std::cout << alpha[2] << "\t";
	std::cout << alpha[3] << "\t";
	std::cout << alpha[4] << "\t";
	std::cout << "( " << beta+alpha[1]+alpha[2]+alpha[3]+alpha[4] << " )\n";
	
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
	// Create a temporary copy of Sn just in case it gets thrashed in the Gaussian Elimination
	double Stmp[5][5];
	for (int i=0;i<=4;i++) {
		for (int j=0;j<=4;j++) {
			Stmp[i][j] = Sn[i][j];
		}
	}

	double* A = new double[4];
	double* B = new double[4];
	for (int i=maxi-order+1; i<= maxi; i++) {
		B[i-(maxi-order+1)] = q[i];
	}
	solveBC(A, B, Stmp);
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

void RecursiveFilter::solveBC(double* A, double* B, double S[5][5]) 
{
	
	int n = order-1;
	for(int j=0;j<=n-1;j++) {
		for (int i=j+1;i<=n;i++) {
			S[i][j]=S[i][j]/S[j][j];
		}
		for	(int i=j+1;i<=n;i++) {
			for (int k=j+1;k<=n;k++) {
				S[i][k]=S[i][k]-S[i][j]*S[j][k];
			}
			B[i]=B[i]-S[i][j]*B[j];
		}
	}
	A[n]=B[n]/S[n][n];
	for(int j=n-1;j>=0;j--) {
		A[j]=B[j];
		for(int k=n;k>=j+1;k--) {
			A[j]=A[j]-A[k]*S[j][k];
		}
		A[j]=A[j]/S[j][j];
	}
}

void RecursiveFilter::getAnisotropicFilterCoefficients(const double* tau, const int& arr)
{
	double sigma = (lengthScale*lengthScale)/2;
	double b[5][5];
	
	// Extend the arrLength for the Cholesky decomp?
	int arrLength = arr;
	double* htau = new double[arrLength];
	for (int i=0; i< arrLength-1; i++) {
		htau[i] = (tau[i+1] + tau[i])/2;
		std::cout << htau[i] << std::endl;
	}
	double** K1 = new double*[arrLength];
	double** K2 = new double*[arrLength];
	double** K3 = new double*[arrLength];
	double** K4 = new double*[arrLength];
	double** D = new double*[arrLength];
	double** sqv = new double*[arrLength];
	for (int i=0; i<arrLength;i++) {
		K1[i] = new double[arrLength];
		K2[i] = new double[arrLength];
		K3[i] = new double[arrLength];
		K4[i] = new double[arrLength];
		D[i] = new double[arrLength];
		sqv[i] = new double[arrLength];
	}
	double* p= new double[arrLength];
	abeta = new double[arrLength];
	
	for (int a = 1; a <5; a++) {
		aalpha[a] = new double[arrLength];
	}
	
	double coeff[5];
	
	// b coefficients
	b[1][1] = 1;
	b[1][2] = 1./12.;
	b[1][3] = 1./90.;
	b[1][4] = 1./560.;
	b[2][2] = 1.;
	b[2][3] = 1./6.;
	b[2][4] = 7./240.;
	b[3][3] = 1.;
	b[3][4] = 1./4.;
	b[4][4] = 1;
	
    // Compute polynomial coefficients
	coeff[1] = b[1][1]*sigma;
	coeff[2] = b[1][2]*sigma + b[2][2]*(pow(sigma,2))/factorial(2);
	coeff[3] = b[1][3]*sigma + b[2][3]*(pow(sigma,2))/factorial(2)
	+ b[3][3]*(pow(sigma,3))/factorial(3);
	coeff[4] = b[1][4]*sigma + b[2][4]*(pow(sigma,2))/factorial(2)
	+ b[3][4]*(pow(sigma,3))/factorial(3) + b[4][4]*(pow(sigma,4))/factorial(4);
	
	// Set up K
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			K1[i][j] = 0.;
			sqv[i][j] = 0.;
		}
		
		/*if (i > 0) K1[i][i-1] = -(sigma*(tau[i-1]+tau[i])/2)/sqrt(sigma*tau[i]*sigma*tau[i-1]);
		if ((i > 0) and (i<(arrLength-1)))
			K1[i][i] = ((sigma*(tau[i-1]+tau[i])/2) + (sigma*(tau[i+1]+tau[i])/2))/(sigma*tau[i]);
		if (i<(arrLength-1)) K1[i][i+1] = -(sigma*(tau[i+1]+tau[i])/2)/sqrt(sigma*tau[i]*sigma*tau[i+1]);
		sqv[i][i] = sqrt(tau[i]); */
		
		if (i > 1) {
			K1[i][i-1] = -(sigma*(tau[i]-tau[i-1]))/sqrt(sigma*(htau[i]-htau[i-1])*sigma*(htau[i-1]-htau[i-2]));
		} else if (i == 1) {
			K1[i][i-1] = -(sigma*(tau[i]-tau[i-1]))/sqrt(sigma*(htau[i]-htau[i-1])*sigma*(htau[i-1]+htau[i-1]));
		}
		if ((i > 0) and (i<(arrLength-1)))
			K1[i][i] = (sigma*(tau[i]-tau[i-1]) + sigma*(tau[i+1]-tau[i]))/(sigma*(htau[i]-htau[i-1]));
		if (i<(arrLength-2)) {
			K1[i][i+1] = -(sigma*(tau[i+1]-tau[i]))/sqrt(sigma*(htau[i]-htau[i-1])*sigma*(htau[i+1]-htau[i]));
		} else {
			K1[i][i+1] = -(sigma*(tau[i]-tau[i-1]))/sqrt(sigma*(htau[i]-htau[i-1])*sigma*(htau[i]-htau[i-1]));
		}
		if (i == 0) {
			sqv[i][i] = sqrt(htau[i]*2);
		} else if (i == (arrLength-1)) {
			sqv[i][i] = sqrt(htau[i-1]*2);
		} else {
			sqv[i][i] = sqrt(htau[i]-htau[i-1]);
		}
	}
	K1[0][0] = 2.;
	K1[0][1] = K1[1][0];
	K1[arrLength-1][arrLength-1] = 2;
	K1[arrLength-1][arrLength-2] = K1[arrLength-2][arrLength-1];
	sqv[arrLength-1][arrLength-1] = sqv[arrLength-2][arrLength-2];
	// Square it
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			K2[i][j] = 0;
		}
	}
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			for (int r=0; r< arrLength; r++) {
				K2[i][j] += K1[i][r]*K1[r][j];
			}
		}
	}
	
	// Cube it
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			K3[i][j] = 0;

		}
	}
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			for (int r=0; r< arrLength; r++) {
				K3[i][j] += K1[i][r]*K2[r][j];
			}
		}
	}
	
	// 4th it
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			K4[i][j] = 0;
		}
	}
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			for (int r=0; r< arrLength; r++) {
				K4[i][j] += K1[i][r]*K3[r][j];
			}
		}
	}
	
	// 4th it
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			K4[i][j] = 0;
		}
	}
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			for (int r=0; r< arrLength; r++) {
				K4[i][j] += K1[i][r]*K3[r][j];
			}
		}
	}
	
	// Now calculate D
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			D[i][j] = coeff[1]*K1[i][j] + coeff[2]*K2[i][j] + coeff[3]*K3[i][j] + coeff[4]*K4[i][j];
			//std::cout << D[i][j] << " ";
		}
		//std::cout << std::endl;
	}
	
	// Get the filter coefficients
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			if (K1[i][j] != 0)
				std::cout << K1[i][j] << "\t";
		}
		std::cout << std::endl;
	}	
	
	// Cholesky decomp
	for (int i=0;i<arrLength;i++) {
		for (int j=i;j<arrLength;j++) {
			double sum=D[i][j];
			for (int k=i-1;k>=0;k--) {
				sum -= D[i][k]*D[j][k];
			}
			if (i == j) {
				if (sum <= 0.0) { 
					std::cout << "cholesky failed at i,j sum\n";
					break;
				} else {
					p[i] = sqrt(sum);
				}
			} else {
				D[j][i]=sum/p[i];
				if (p[i] == 0.) { 
					std::cout << "Problem! " << i << "\t" << j << "\n";
				}
			}
		}
	}
	
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			if (D[i][j] != 0)
				std::cout << D[i][j] << "\t";
		}
		std::cout << std::endl;
	}	
	
	
	/* Multiply by sqv (reuse K1)
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			K1[i][j] = 0.;
		}
	}
	
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			for (int r=0; r< arrLength; r++) {
				K1[i][j] += D[i][r]*sqv[r][j];
			}
		}
	}
	
	// Multiply by 1/sqv (reuse K1)
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			D[i][j] = 0.;
			if (i == 0) {
				sqv[i][i] = 1/sqrt(htau[i]*2);
			} else if (i == (arrLength-1)) {
				sqv[i][i] = 1/sqrt(htau[i-1]*2);
			} else {
				sqv[i][i] = 1/sqrt(htau[i]-htau[i-1]);
			}
		}
	}
	sqv[arrLength-1][arrLength-1] = sqv[arrLength-2][arrLength-2];
	
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			for (int r=0; r< arrLength; r++) {
				D[i][j] += sqv[i][r]*K1[r][j];
			}
		}
	} */
	
	
	for (int i=0; i<arrLength; i++) {
		for (int j=0; j<arrLength; j++) {
			if (D[i][j] != 0)
			std::cout << D[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	abeta[0]= 1/p[0];
	aalpha[1][0]=0;
	aalpha[2][0]=0;
	aalpha[3][0]=0;
	aalpha[4][0]=0;
	std::cout << "0\t" << abeta[0]<< "\t" << aalpha[1][0] << "\t" << aalpha[2][0] << "\t" 
	<< aalpha[3][0] << "\t" << aalpha[4][0]  << "\t(" << abeta[0] + aalpha[1][0] + aalpha[2][0] + aalpha[3][0] + aalpha[4][0] << ")\n";
		
	abeta[1]= 1/p[1];
	aalpha[1][1]=-abeta[1]*D[1][0];
	aalpha[2][1]=0;
	aalpha[3][1]=0;
	aalpha[4][1]=0;
	std::cout << "1\t" << abeta[1]<< "\t" << aalpha[1][1] << "\t" << aalpha[2][1] << "\t" 
	<< aalpha[3][1] << "\t" << aalpha[4][1]  << "\t(" << abeta[1] + aalpha[1][1] + aalpha[2][1] + aalpha[3][1] + aalpha[4][1] << ")\n";
	
	abeta[2]= 1/p[2];	
	aalpha[1][2]=-abeta[2]*D[2][1];
	aalpha[2][2]=-abeta[2]*D[2][0];
	aalpha[3][2]=0;
	aalpha[4][2]=0;
	std::cout << "2\t" << abeta[2]<< "\t" << aalpha[1][2] << "\t" << aalpha[2][2] << "\t" 
	<< aalpha[3][2] << "\t" << aalpha[4][2]  << "\t(" << abeta[2] + aalpha[1][2] + aalpha[2][2] + aalpha[3][2] + aalpha[4][2] << ")\n";
	
	abeta[3]= 1/p[3];
	aalpha[1][3]=-abeta[2]*D[3][2];
	aalpha[2][3]=-abeta[2]*D[3][1];
	aalpha[3][3]=-abeta[2]*D[3][0];
	aalpha[4][3]=0;
	std::cout << "3\t" << abeta[3]<< "\t" << aalpha[1][3] << "\t" << aalpha[2][3] << "\t" 
	<< aalpha[3][3] << "\t" << aalpha[4][3]  << "\t(" << abeta[3] + aalpha[1][3] + aalpha[2][3] + aalpha[3][3] + aalpha[4][3] << ")\n";

	for (int i=4;i<arrLength;i++) {
		abeta[i] = 1/p[i];
		aalpha[1][i]=-abeta[i]*D[i][i-1];
		aalpha[2][i]=-abeta[i]*D[i][i-2];
		aalpha[3][i]=-abeta[i]*D[i][i-3];
		aalpha[4][i]=-abeta[i]*D[i][i-4];
		double sum = abeta[i] + aalpha[1][i] + aalpha[2][i] + aalpha[3][i] + aalpha[4][i];
		std::cout << i<< "\t" << abeta[i]<< "\t" << aalpha[1][i] << "\t" << aalpha[2][i] << "\t" 
		<< aalpha[3][i] << "\t" << aalpha[4][i]  << "\t(" << sum << ")\n";
	}
	
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
				L[i+j][i] = -aalpha[i+j][arr-1];
				U[i][i+j] = aalpha[order-j][arr-1];
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
	
	for (int i=0; i<arrLength;i++) {
		delete[] K1[i];
		delete[] K2[i];
		delete[] K3[i];
		delete[] K4[i]; 
		delete[] D[i];
		delete[] sqv[i];
	}
	delete[] tau;
	delete[] K1;
	delete[] K2;
	delete[] K3;
	delete[] K4; 
	delete[] D;
	delete[] sqv;
	delete[] p;
	
}

bool RecursiveFilter::aniFilterArray(double* array, const int& arrLength)
{
	int maxi = arrLength-1;
	double* p = array;
	double* q = new double[arrLength];
	double* s = new double[arrLength];
	
	for (int i=0; i<= maxi; i++) {
		q[i] = 0;
		s[i] = 0;
	}
	
	q[0]=abeta[0]*p[0];
	q[1]=abeta[1]*p[1] + aalpha[1][1]*q[0];
	q[2]=abeta[2]*p[2] + aalpha[1][2]*q[1] + aalpha[2][2]*q[0];
	q[3]=abeta[3]*p[3] + aalpha[1][3]*q[2] + aalpha[2][3]*q[1] + aalpha[3][3]*q[0];	
	for (int i=order; i<= maxi; i++) {
		q[i] = abeta[i]*p[i] + aalpha[1][i]*q[i-1]
		+ aalpha[2][i]*q[i-2] + aalpha[3][i]*q[i-3] + aalpha[4][i]*q[i-4];
	}
	
	for (int i=0; i<= maxi; i++) {
		std::cout << i << "\t" << q[i] << "\n";
	}
	
    // Invert Sn
	// Create a temporary copy of Sn just in case it gets thrashed in the Gaussian Elimination
	double Stmp[5][5];
	for (int i=0;i<=4;i++) {
		for (int j=0;j<=4;j++) {
			Stmp[i][j] = Sn[i][j];
		}
	}
	
	double* A = new double[4];
	double* B = new double[4];
	for (int i=maxi-order+1; i<= maxi; i++) {
		B[i-(maxi-order+1)] = q[i];
	}
	solveBC(A, B, Stmp);
	for (int i=maxi; i>= (maxi-order+1); i--) {
		s[i] = A[i-(maxi-order+1)];
	}
	for (int i=maxi-order;i>=0;i--) {
		s[i] = abeta[i]*q[i] + aalpha[1][i]*s[i+1]
		+ aalpha[2][i]*s[i+2] + aalpha[3][i]*s[i+3] + aalpha[4][i]*s[i+4];
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

