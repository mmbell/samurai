/*
 *  ParametricVortex.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "ParametricVortex.h"
#include <iostream>
#include <fstream>
#include <cmath>

ParametricVortex::ParametricVortex(const double* Vin, const double* Radin, const int& numradii)
{
	
	// Allocate memory for the simplex vertices	
    vertex = new float*[3];
    vertex[0] = new float[3];
    vertex[1] = new float[3];
    vertex[2] = new float[3];
	vertex[3] = new float[3];
    VT = new float[3];
    vertexSum = new float[3]; 
	
	// Assign given winds to internal variables
	winds = Vin;
	radii = Radin;
	maxr = numradii;
	
	// Find Vmax and Rmax from the data
	vmax = 0;
	for (int r=0;r < maxr;r++) {
		if (winds[r] > vmax) {
			vmax = winds[r];
			rmw = radii[r];
		}
	}
	
	// Assign some constants
	X2 = 25;
	deltar = 20;
	
	// Perform a simplex search to fit the parameters	
	fitParameters();

	
}

ParametricVortex::~ParametricVortex()
{
	delete[] vertex[0];
	delete[] vertex[1];
	delete[] vertex[2];
	delete[] vertex[3];
	delete[] vertex;
	delete[] VT;
	delete[] vertexSum;
}

void ParametricVortex::fitParameters()
{

    std::cout << "Fitting winds to Willoughby Parametric Vortex\n";
	
		
	// X1, A, n
	// First guess from regressions without latitude
	float X1 = 317.1 - 2.026*vmax;
	float n = 0.4067 + 0.0144*vmax;
	float A = 0.0696 + 0.0049*vmax;
	vertex[0][0] = X1;
	vertex[0][1] = A;
	vertex[0][2] = n;
	vertex[1][0] = X1+X1*.20;
	vertex[1][1] = A;
	vertex[1][2] = n;
	vertex[2][0] = X1;
	vertex[2][1] = A+A*.20;
	vertex[2][2] = n;
	vertex[3][0] = X1;
	vertex[3][1] = A;
	vertex[3][2] = n+n*.20;
	
	vertexSum[0] = 0;
	vertexSum[1] = 0;
	vertexSum[2] = 0;
	
	for (int v=0; v <= 3; v++) {	
		VT[v] = getS2(vertex[v][0],vertex[v][1],vertex[v][2]);
	}
	
	// Run the simplex search loop
	getVertexSum();
	int numIterations = 0;
	int low = 0;
	int mid = 0;
	int high = 0;
	
	for(;;) {
		
		low = 0;
		// Sort the initial guesses
		high = VT[0] > VT[1] ? (mid = 1,0) : (mid = 0,1);
		for (int v=0; v<=3; v++) {
			if (VT[v] <= VT[low]) low = v;
			if (VT[v] > VT[high]) {
				mid = high;
				high = v;
			} else if (VT[v] > VT[mid] && v != high) mid = v;
		}
		
		// Check convergence
		float epsilon = 2.0 * fabs(VT[high]-VT[low])/(fabs(VT[high]) + fabs(VT[low]) + 1.0e-10);
		if (epsilon < 0.05) {
			
			// Converged
			std::cout << "Parameters Converged in " << numIterations << " iterations" << std::endl;
			std::cout << "\tS2 = " << VT[high]  << std::endl;
			std::cout << "\tX1 = " << vertex[high][0]  << std::endl;
			std::cout << "\tX2 = " << X2  << std::endl;
			std::cout << "\tA = " << vertex[high][1]  << std::endl;
			std::cout << "\tn = " << vertex[high][2]  << std::endl;
			
			X1fit = vertex[high][0];
			Afit = vertex[high][1];
			nfit = vertex[high][2];
			
			
			// Determine R1 from n and X1
			float R1, eps, w, dw;
			// float residual = (n*X1)/(n*X1 + rmw);
			float residual = (n*((1-A)*X1 + A*X2))/(n*((1-A)*X1 + A*X2) + rmw);
			R1=rmw/2;
			eps = (rmw-R1)/(deltar);
			for (int j=0;j<20;j++) {
				w = 126*pow(eps,5) - 420*pow(eps,6) + 540*pow(eps,7) - 315*pow(eps,8) + 70*pow(eps,9)-residual;
				dw = 5*126*pow(eps,4) - 6*420*pow(eps,5) + 7*540*pow(eps,6) - 8*315*pow(eps,7) + 9*70*pow(eps,8);
				eps = eps - w/dw;
				if (fabs(w/dw) < .001) { 
					break;
				}
			}
			std::cout << "\tR1 = " << R1  << std::endl;
			std::cout << "\tR2 - R1 = " << deltar << std::endl;
			R1fit = rmw-(eps*deltar);
			
			if (n < 1) {
				// Remove the singularity
				Lc = sqrt((rmw*rmw - (pow(0.25,1-nfit)*rmw*rmw/16))/(1-pow(0.25,1-nfit)));
				Zetac = 1/((Lc/2)*(rmw/Lc-pow((rmw/Lc),3)/2));
			} else { Lc = 0; Zetac = 0; }
			break;
		}
		
		// Check iterations
		if (numIterations > 100) {
			std::cout << "Too many iterations (>100)";
			break;
		}
		
		numIterations += 2;
		// Reflection
		float VTtest = simplexTest(low, -1.0);
		if (VTtest >= VT[high])
			// Better point than highest, so try expansion
			VTtest = simplexTest(low, 2.0);
		else if (VTtest <= VT[mid]) { 
			// Worse point than second highest, so try contraction
			float VTsave = VT[low];
			VTtest = simplexTest(low, 0.5);
			if (VTtest <= VTsave) {
				for (int v=0; v<=3; v++) {
					if (v != high) {
						for (int i=0; i<=2; i++)
							vertex[v][i] = vertexSum[i] = 0.5*(vertex[v][i] + vertex[high][i]);
						VT[v] = getS2(vertex[v][0],vertex[v][1],vertex[v][2]);
					}
				}
				numIterations += 2;
				getVertexSum();
			}
		} else --numIterations;
		
	}


	
}

inline void ParametricVortex::getVertexSum()
{
	
	float sum;
	int v;
	for (int i=0; i<=2; i++) {
		for (sum = 0.0, v=0; v<=3; v++)
			sum += vertex[v][i];
		vertexSum[i] = sum;
	}
}

float ParametricVortex::simplexTest(int& low, float factor)
{
	
	// Test a simplex vertex
	float VTtest = -999;
	float* vertexTest = new float[2];
	float factor1 = (1.0 - factor)/2;
	float factor2 = factor1 - factor;
	for (int i=0; i<=2; i++)
		vertexTest[i] = vertexSum[i]*factor1 - vertex[low][i]*factor2;
	VTtest = getS2(vertexTest[0],vertexTest[1],vertexTest[2]);
	
	// If its a better point than the worst, replace it
	if (VTtest > VT[low]) {
		VT[low] = VTtest;
		for (int i=0; i<=2; i++) {
			vertexSum[i] += vertexTest[i]-vertex[low][i];
			vertex[low][i] = vertexTest[i];
		}
	}
	delete[] vertexTest;
	return VTtest;
	
}

float ParametricVortex::getS2(const float& X1, const float& A, const float& n)
{
	
	// Determine R1 from n and X1
	float R1, R2, eps, w, dw;
	// float residual = (n*X1)/(n*X1 + rmw);
	float residual = (n*((1-A)*X1 + A*X2))/(n*((1-A)*X1 + A*X2) + rmw);
	R1=rmw/2;
	eps = (rmw-R1)/(deltar);
	for (int j=0;j<20;j++) {
		w = 126*pow(eps,5) - 420*pow(eps,6) + 540*pow(eps,7) - 315*pow(eps,8) + 70*pow(eps,9)-residual;
		dw = 5*126*pow(eps,4) - 6*420*pow(eps,5) + 7*540*pow(eps,6) - 8*315*pow(eps,7) + 9*70*pow(eps,8);
		eps = eps - w/dw;
		if (fabs(w/dw) < .001) { 
			// std::cout << "Convergence reached on R1 at " << R1;
			break;
		}
	}
	R1 = rmw-(eps*deltar);
	R2 = R1 + deltar;
	
	// Loop through the radii
	float S2 = 0;
	for (int rad=0; rad < maxr; rad++) {
		if (winds[rad] != -999) {
			float r = radii[rad];
			if (r <= R1) {
				Vg = vmax*pow((r/rmw),n);
			} else if ((r >= R1) and (r <= R2)) {
				eps = (r -R1)/(deltar);
				w = 126*pow(eps,5) - 420*pow(eps,6) + 540*pow(eps,7) - 315*pow(eps,8) + 70*pow(eps,9);
				Vg = vmax*pow((r/rmw),n)*(1-w) + vmax*((1-A)*exp(-(r-rmw)/X1) + A*exp(-(r-rmw)/X2))*w;
			} else {
				Vg = vmax*((1-A)*exp(-(r-rmw)/X1) + A*exp(-(r-rmw)/X2));
			}
		}
		S2 = S2 + pow((winds[rad] - Vg),2);
	}
	
	return (-S2);	
}


float ParametricVortex::getWindAtRadiusKM(const float& radius)
{
	
	float R2 = R1fit + deltar;
	float Vg;
	
	if (radius <= R1fit) {
		if (nfit >= 1) {
			Vg = vmax*pow((radius/rmw),nfit);
		} else {
			Vg = (Lc*Zetac*0.5)*(radius/Lc-pow((radius/Lc),3)/2);
		}
	} else if ((radius >= R1fit) and (radius <= R2)) {
		float eps = (radius -R1fit)/(deltar);
		float w = 126*pow(eps,5) - 420*pow(eps,6) + 540*pow(eps,7) - 315*pow(eps,8) + 70*pow(eps,9);
		if (nfit >= 1) {
			Vg = vmax*pow((radius/rmw),nfit)*(1-w) 
				+ vmax*((1-Afit)*exp(-(radius-rmw)/X1fit)
				+ Afit*exp(-(radius-rmw)/X2))*w;
		} else {
			Vg = (Lc*Zetac*0.5)*(radius/Lc-pow((radius/Lc),3)/2)*(1-w) 
				+ vmax*((1-Afit)*exp(-(radius-rmw)/X1fit)
				+ Afit*exp(-(radius-rmw)/X2))*w;
		}
	} else {
		Vg = vmax*((1-Afit)*exp(-(radius-rmw)/X1fit)
			+ Afit*exp(-(radius-rmw)/X2));
	}

	return Vg;

}

float ParametricVortex::getX1()
{
	return X1fit;
}

float ParametricVortex::getn()
{
	return nfit;
}

float ParametricVortex::getA()
{
	return Afit;
}

