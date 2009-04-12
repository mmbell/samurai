/*
 *  CostFunctionRZ_CPU.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionRZ_CPU.h"
#include <cmath>

/* First, compute HCq and Hxb via recursive filter on q, followed by loop over M observations
 yielding 2 y-estimate M x 1 vectors. Solve for innovation vector d by subtracting from y.
 For function, subtract from HCq and take inner product (with R-1), add to i.p. of q vector.
 For gradient, loop over M again, this time summing product of operator*basis*oberror-1*y-hat
 and operator*basis*oberror1*d for each coefficient N, yielding an N x 1 vector. Filter this
 vector and the innovation vector. Add 1+vector1+vector2 to get the gradient. Treat the penalty
 constraints as observations. */

CostFunctionRZ_CPU::CostFunctionRZ_CPU(const int& numObs, const int& stateSize)
	: CostFunction(numObs, stateSize)
{
}

CostFunctionRZ_CPU::~CostFunctionRZ_CPU()
{
}

void CostFunctionRZ_CPU::finalize()
{	

	delete iFilter;
	delete jFilter;
	delete[] currState;
	delete[] currGradient;
	delete[] innovation;
	delete[] tempState;
	delete[] tempGradient;
	delete[] HCq;
	delete[] xt;
	delete[] df;
	delete[] stateA;
	delete[] stateB;
	delete[] CTHTd;
	delete[] iTemp;
	delete[] jTemp;
	delete[] iL;
	delete[] jL;
	
}

void CostFunctionRZ_CPU::initialize(const real& imin, const real& imax, const int& idim,
									const real& jmin, const real& jmax, const int& jdim,
									const real* ia, const real* ja, real* bgU, real* obs,
									const unsigned int* IXdim, const vector<real>* IX,
									const unsigned int* JXdim, const vector<real>* JX)
{

	// Initialize number of control variables -- one less than physical variables
	varDim = 5;
	
	// Initialize background errors and filter scales
	bgError[0] = 25.;
	bgError[1] = 500.;
	bgError[2] = 150.;
	bgError[3] = 5.;
	bgError[4] = 0.025;	

	rhoBase = 1.1646;
	rhoInvScaleHeight = 1.068e-4;
	
	// Assign local object pointers
	bgFields = bgU;
	obsVector = obs;
	iMin = imin;
	iMax = imax;
	iDim = idim;
	jMin = jmin;
	jMax = jmax;
	jDim = jdim;
	iA = ia;
	jA = ja;
	IXDim = IXdim;
	IXform = IX;
	if (JXdim) JXDim = JXdim;
	if (JX) JXform = JX;
	DI = (iMax - iMin) / (iDim - 1);
	DIrecip = 1./DI;
	DJ = (jMax - jMin) / (jDim - 1);
    DJrecip = 1./DJ;
	
	// Set up the recursive filter
	real iFilterScale = 2;
	real jFilterScale = 2;
	iFilter = new RecursiveFilter(4,iFilterScale);
	jFilter = new RecursiveFilter(4,jFilterScale);

	/* Test the filter
	real* p = new real[iDim];
	for (int i = 0; i < iDim; i++) {
		p[i] = 0;
	}
	
	p[iDim/2] = 1.;
	real Pi = 3.141592653589793;
	iFilter->filterArray(p, iDim);
	for (int i = 0; i < iDim; i++) {
		real ibasis =  Basis(i, iDim/2, iDim-1, iMin, DI, DIrecip, 2);
		cout << i << "\t" << p[i] << "\t" << ibasis << endl;
	} */
	
	/*iFilterCoeff = iFilter->getFilterCoefficients();
	iFilterBC = iFilter->getSn();
	jFilterCoeff = jFilter->getFilterCoefficients();
	jFilterBC = jFilter->getSn(); */

	
	// Allocate memory for the needed arrays
	// These are common to all CostFunctions
	currState = new real[nState];
	currGradient = new real[nState];
	tempState = new real[nState];
	tempGradient = new real[nState];
	xt = new real[nState];
	df = new real[nState];
	innovation = new real[mObs];	

	// These are local to this one
	HCq = new real[mObs];
	bgState = new real[nState];
	stateB = new real[nState];
	stateA = new real[nState];
	CTHTd = new real[nState];
	iTemp = new real[iDim];
	jTemp = new real[jDim];
	
	// Set up the spline matrices
	setupSplines();
	initState();
}	

void CostFunctionRZ_CPU::initState()
{

	// Clear the state vector
	cout << "Initializing state vector..." << endl;
	for (int n = 0; n < nState; n++) {
		currState[n] = 0.0;
		currGradient[n] = 0.0;
		tempState[n] = 0.0;
		tempGradient[n] = 0.0;
		xt[n] = 0.0;
		df[n] = 0.0;
		bgState[n] = 0.0;
		stateA[n] = 0.0;
		stateB[n] = 0.0;
	}
	
	
	// SB Transform
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
			for (int imu = -1; imu <= 1; imu += 2) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
				int ii = (int)((i - iMin)*DIrecip);
				for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
					for (int jmu = -1; jmu <= 1; jmu += 2) {
						real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5));
						int jj = (int)((j - jMin)*DJrecip);
						for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
							if ((iNode < 0) or (iNode >= iDim)) continue;
							for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
								if ((jNode < 0) or (jNode >= jDim)) continue;
								real iBC, jBC;
								if (var <= 1) {
									iBC = 4;
								} else {
									iBC = 2;
								}
								if (var == 1) {
									jBC = 4;
								} else {
									jBC = 2;
								}
								real im = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, iBC);
								real jm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, jBC);
								int bgJ = jIndex*2 + (jmu+1)/2;
								int bgI = iIndex*2 + (imu+1)/2;
								bgState[varDim*iDim*jNode +varDim*iNode + var] += 
									0.25 * bgFields[varDim*(iDim-1)*2*bgJ +varDim*bgI + var] * im * jm;
							}
						}	
					}
				}
			}
		}	
	}
	
	// Compute and display the variable RMS and BG errors for reference
	for (int var = 0; var < varDim; var++) {
		double varScale = 0;
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				varScale += bgState[varDim*iDim*jIndex +varDim*iIndex + var] * bgState[varDim*iDim*jIndex +varDim*iIndex + var];
			}
		}
		varScale = sqrt(varScale/(iDim*jDim));
		if (varScale) {
			double errPct = 100*bgError[var]/varScale;
			cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var] 
			<< " ( " << errPct << " %)" << endl;
		} else {
			cout << "Variable " << var << " RMS = " << varScale << "\t BG Error = " << bgError[var] 
			<< " ( Infinite! %)" << endl;
		}
	}
	
	// SA transform = bg B's -> bg A's
	SAtransform(bgState, stateA);
	
	// d = y - HXb
	calcInnovation();
	
	// HTd
	calcHTranspose(innovation, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtransform(stateA, stateB);
	
	for (int var = 0; var < varDim; var++) {
		//FI & D
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iTemp[iIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var] * bgError[var];
			}
			iFilter->filterArray(iTemp, iDim);
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex]; 
			}
		}
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				CTHTd[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
	} 
	
	/* D
	for (int var = 0; var < varDim; var++) {
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				CTHTd[varDim*iDim*jIndex +varDim*iIndex + var] = stateB[varDim*iDim*jIndex +varDim*iIndex + var] * bgError[var];
			}
		}
	} */
	
}	

double CostFunctionRZ_CPU::funcValue(double* state)
{
	// Update the Y hat vector
	updateHCq(state);

	double qIP, obIP;
	qIP = 0;
	obIP = 0;
	// Compute inner product of state vector
	for (int n = 0; n < nState; n++) {
		qIP += state[n]*state[n];
	}
		
	// Subtract d from HCq to yield mObs length vector and compute inner product
	#pragma omp parallel for reduction(+:obIP)
	for (int m = 0; m < mObs; m++) {
		obIP += (HCq[m]-innovation[m])*(obsVector[m*10+1])*(HCq[m]-innovation[m]);
	}
	double J = 0.5*qIP + 0.5*obIP;
	return J;
	
}

void CostFunctionRZ_CPU::funcGradient(double* state, double* gradient)
{
	
	updateHCq(state);
	
	// HTHCq
	calcHTranspose(HCq, stateA);
	
	// S^T (Inverse SA transform) yield B's, put it in the tempState
	SAtransform(stateA, stateB);
		
	for (int var = 0; var < varDim; var++) {
		//FI & D
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iTemp[iIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var] * bgError[var];
			}
			iFilter->filterArray(iTemp, iDim);
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex]; 
			}
		}
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
	}
	
	/* D
	for (int var = 0; var < varDim; var++) {
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] *= bgError[var];
			}
		}
	} */
	
	for (int n = 0; n < nState; n++) {
		gradient[n] = state[n] + stateB[n] - CTHTd[n];
	}
	
	
}

void CostFunctionRZ_CPU::updateHCq(double* state)
{
		
	for (int var = 0; var < varDim; var++) {
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = state[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
		//FI
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iTemp[iIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			iFilter->filterArray(iTemp, iDim);
			// D
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex] * bgError[var]; 
			}
		}
	}
	
	/* D
	for (int var = 0; var < varDim; var++) {
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = state[varDim*iDim*jIndex +varDim*iIndex + var] * bgError[var];
			}
		}
	} */
	

	// S (SA transform) yield A's, put it in the tempState
	SAtransform(stateB, stateA);

	// H
	#pragma omp parallel for
	for (int m = 0; m < mObs; m++) {
		int mi = m*10;
		real w1 = obsVector[mi+2];
		real w2 = obsVector[mi+3];
		real w3 = obsVector[mi+4];
		real w4 = obsVector[mi+5];
		real w5 = obsVector[mi+6];
		real w6 = obsVector[mi+7];		
		real i = obsVector[mi+8];
		real invI;
		if (i != 0) {
			invI = 1./i;
		} else {
			invI = 0.;
		}
		real j = obsVector[mi+9];
		real rhoBar = rhoBase*exp(-rhoInvScaleHeight*j);
		real tempsum = 0;
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		real ibasis = 0;
		real jbasis = 0;
		real idbasis = 0;
		real jdbasis = 0;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {				
				if ((iNode < 0) or (iNode >= iDim) or (jNode < 0) or (jNode >= jDim)) continue;
				// If conditions may cause parallel warps to diverge -- maybe actually faster to just calculate them all but need to test
				if(w1) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * w1 * invI * rhoBar;
				}
				if(w2 or w3) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					idbasis = DBasis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jdbasis = DBasis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					float coeff = stateA[varDim*iDim*jNode +varDim*iNode + 1];
					tempsum += coeff * ibasis * (-jdbasis) * w2 * 1e3 * invI * rhoBar;
					tempsum += coeff * idbasis * jbasis * w3 * invI * rhoBar;
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis * w4 * rhoBar;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * w5 * rhoBar;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * w6;
				}
			}
		}
		HCq[m] = tempsum;
	}	

	
}

void CostFunctionRZ_CPU::updateBG()
{

	for (int var = 0; var < varDim; var++) {
		//FJ
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jTemp[jIndex] = currState[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			jFilter->filterArray(jTemp, jDim);
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = jTemp[jIndex];
			}
		}
		//FI
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iTemp[iIndex] = stateB[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			iFilter->filterArray(iTemp, iDim);
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				// D
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = iTemp[iIndex] * bgError[var]; 
			}
		}
	} 
	
	/* D
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				stateB[varDim*iDim*jIndex +varDim*iIndex + var] = currState[varDim*iDim*jIndex +varDim*iIndex + var] * bgError[var];
			}
		}
	} */
	
	// In BG update we are directly summing B's
	ofstream bstream("BAnalysis.out");
	bstream << "Variable\tI\tJ\tBackground\tAnalysis\tIncrement\n";
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				bstream << var << "\t" << iIndex << "\t" << jIndex << "\t"; 
				bstream << bgState[varDim*iDim*jIndex +varDim*iIndex + var] << "\t";
				bgState[varDim*iDim*jIndex +varDim*iIndex + var] += stateB[varDim*iDim*jIndex +varDim*iIndex + var];
				bstream << bgState[varDim*iDim*jIndex +varDim*iIndex + var] << "\t"; 
				bstream << stateB[varDim*iDim*jIndex +varDim*iIndex + var] << endl;
			}
		}
	}
	
	// S (SA transform) yield A's, put it in the tempState
	SAtransform(bgState, stateA);
	
	// H --> to Mish for output
	ofstream fluxstream("FluxAnalysis.out");
	fluxstream << "Radius\tHeight\trhoM\trhoE\tu\tv\tw\tqv\trho\tPsi\tBG1\tBG2\n";
	fluxstream.precision(10);
	real rhoBase = 1.1646;
	real rhoInvScaleHeight = 1.068e-4;
	real CoriolisF = 6e-5;
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		for (int ihalf = 0; ihalf <=1; ihalf++) {
			for (int imu = -ihalf; imu <= ihalf; imu++) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5*ihalf));
				real invI;
				if (i != 0) {
					invI = 1./i;
				} else {
					invI = 0.;
				}
				if (i > (iDim-1)*DI) continue;
				for (int jIndex = 0; jIndex < jDim; jIndex++) {
					for (int jhalf =0; jhalf <=1; jhalf++) {
						for (int jmu = -jhalf; jmu <= jhalf; jmu++) {
							real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5*jhalf));
							real rhoBar = rhoBase*exp(-rhoInvScaleHeight*j);
							if (j > (jDim-1)*DJ) continue;	
							int ii = (int)((i - iMin)*DIrecip);
							int jj = (int)((j - jMin)*DJrecip);
							real ibasis = 0;
							real jbasis = 0;
							real idbasis = 0;
							real jdbasis = 0;
							real rhov = 0;
							real rhou = 0;
							real rhow = 0;
							real rhoe = 0;
							real rhoqv = 0;
							real rhoprime = 0;
							real psi = 0;
							
							for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
								for (int jNode = jj-1; jNode <= jj+2; ++jNode) {				
									if ((iNode < 0) or (iNode >= iDim) or (jNode < 0) or (jNode >= jDim)) continue;
									ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
									jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
									rhov += stateA[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * invI  * rhoBar;
									
									ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
									jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
									idbasis = DBasis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
									jdbasis = DBasis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
									float coeff = stateA[varDim*iDim*jNode +varDim*iNode + 1];
									rhou += coeff * ibasis * (-jdbasis) * 1e3 * invI * rhoBar;
									rhow += coeff * idbasis * jbasis * invI * rhoBar;
									psi += coeff * ibasis * jbasis;
									
									ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 2);
									jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
									rhoe += stateA[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis * rhoBar;
									rhoqv += stateA[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * rhoBar;
									rhoprime += stateA[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis;
								}
							}
							
							// Get it to relevant variables for the flux calculation
							real rho = rhoBar + rhoprime;
							real v = rhov / rho;
							real u = rhou / rho;
							real w = rhow / rho;
							real KE = 0.5*rho*(v*v + u*u + w*w);
							real rhoM = (i*rhov*1000 + 0.5*CoriolisF*i*i*1e6);
							real rhoE = rhoe*1000 + KE;
							real qv = rhoqv / rho;
							
							// Output it
							fluxstream << scientific << i << "\t" << j << "\t" << rhoM << "\t" << rhoE
							<< "\t" << u << "\t" << v << "\t" << w << "\t" << qv << "\t" << rho << "\t" << psi << "\n";
						}
					}
				}
			}
		}
	}
	
	
	// H -> to QC
	// Write the Obs to a summary text file
	ofstream qcstream("AnalysisQC.out");
	ifstream obstream("./Observations.out");
	ostream_iterator<string> os(qcstream, "\t ");
	*os++ = "Type";
	*os++ = "r";
	*os++ = "z";
	*os++ = "NULL";
	*os++ = "Observation";
	*os++ = "Inverse Error";
	*os++ = "Weight 1";
	*os++ = "Weight 2";
	*os++ = "Weight 3";
	*os++ = "Weight 4";
	*os++ = "Weight 5";
	*os++ = "Weight 6";
	*os++ = "Analysis";
	*os++ = "Background";
	qcstream << endl;
	double temp;
	
	ostream_iterator<double> od(qcstream, "\t ");
	for (int m = 0; m < mObs; m++) {
		int mi = m*10;
		real w1 = obsVector[mi+2];
		real w2 = obsVector[mi+3];
		real w3 = obsVector[mi+4];
		real w4 = obsVector[mi+5];
		real w5 = obsVector[mi+6];
		real w6 = obsVector[mi+7];		
		real i = obsVector[mi+8];
		real invI;
		if (i != 0) {
			invI = 1./i;
		} else {
			invI = 0.;
		}
		real j = obsVector[mi+9];
		real rhoBar = rhoBase*exp(-rhoInvScaleHeight*j);
		real tempsum = 0;
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		real ibasis = 0;
		real jbasis = 0;
		real idbasis = 0;
		real jdbasis = 0;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {				
				if ((iNode < 0) or (iNode >= iDim) or (jNode < 0) or (jNode >= jDim)) continue;
				// If conditions may cause parallel warps to diverge -- maybe actually faster to just calculate them all but need to test
				if(w1) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * w1 * invI * rhoBar;
				}
				if(w2 or w3) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					idbasis = DBasis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jdbasis = DBasis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					float coeff = stateA[varDim*iDim*jNode +varDim*iNode + 1];
					tempsum += coeff * ibasis * (-jdbasis) * w2 * 1e3 * invI * rhoBar;
					tempsum += coeff * idbasis * jbasis * w3 * invI * rhoBar;
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis * w4 * rhoBar;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * w5 * rhoBar;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * w6;
				}
			}
		}
		
		for (int t=0; t<12; t++) {
			obstream >> temp;
			*od++ = temp;
		}
		*od++ = tempsum;
		*od++ = obsVector[mi]-innovation[m];
		qcstream << endl;
		
	}	
	
	
}

void CostFunctionRZ_CPU::calcInnovation()
{
	// Initialize and fill the innovation vector
	cout << "Initializing innovation vector..." << endl;
	for (int m = 0; m < mObs; m++) {
		HCq[m] = 0.0;
		innovation[m] = obsVector[m*10];
	}
	
	real innovationRMS = 0;
	#pragma omp parallel for reduction(+:innovationRMS)
	for (int m = 0; m < mObs; m++) {
		int mi = m*10;
		real w1 = obsVector[mi+2];
		real w2 = obsVector[mi+3];
		real w3 = obsVector[mi+4];
		real w4 = obsVector[mi+5];
		real w5 = obsVector[mi+6];
		real w6 = obsVector[mi+7];		
		real i = obsVector[mi+8];
		real invI;
		if (i != 0) {
			invI = 1./i;
		} else {
			invI = 0.;
		}
		real j = obsVector[mi+9];
		real rhoBar = rhoBase*exp(-rhoInvScaleHeight*j);
		real tempsum = 0;
		int ii = (int)((i - iMin)*DIrecip);
		int jj = (int)((j - jMin)*DJrecip);
		real ibasis = 0;
		real jbasis = 0;
		real idbasis = 0;
		real jdbasis = 0;
		
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			for (int jNode = jj-1; jNode <= jj+2; ++jNode) {				
				if ((iNode < 0) or (iNode >= iDim) or (jNode < 0) or (jNode >= jDim)) continue;
				// If conditions may cause parallel warps to diverge -- maybe actually faster to just calculate them all but need to test
				if(w1) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode] * ibasis * jbasis * w1 * invI * rhoBar;
				}
				if(w2 or w3) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					idbasis = DBasis(iNode, i, iDim-1, iMin, DI, DIrecip, 4);
					jdbasis = DBasis(jNode, j, jDim-1, jMin, DJ, DJrecip, 4);
					float coeff = stateA[varDim*iDim*jNode +varDim*iNode + 1];
					tempsum += coeff * ibasis * (-jdbasis) * w2 * 1e3 * invI * rhoBar;
					tempsum += coeff * idbasis * jbasis * w3 * invI * rhoBar;
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 2);
					jbasis = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, 2);
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 2] * ibasis * jbasis * w4 * rhoBar;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 3] * ibasis * jbasis * w5 * rhoBar;
					tempsum += stateA[varDim*iDim*jNode +varDim*iNode + 4] * ibasis * jbasis * w6;
				}
			}
		}
		innovation[m] -= tempsum;
		innovationRMS += (innovation[m]*innovation[m]);
	}
	if (mObs) innovationRMS /= mObs;
	innovationRMS = sqrt(innovationRMS);
	cout << "Innovation RMS : " << innovationRMS << endl;
	
}	

void CostFunctionRZ_CPU::calcHTranspose(const real* yhat, real* Astate) 
{
	
	// Clear the Astate
	for (int var = 0; var < varDim; var++) {
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = 0.;
			}
		}
	}
	
	// Calculate H Transpose	
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			#pragma omp parallel for
			for (int m = 0; m < mObs; m++) {
				// Sum over obs this time
				// Multiply state by H weights
				int mi = m*10;
				real invError = obsVector[mi+1];
				real w1 = obsVector[mi+2];
				real w2 = obsVector[mi+3];
				real w3 = obsVector[mi+4];
				real w4 = obsVector[mi+5];
				real w5 = obsVector[mi+6];
				real w6 = obsVector[mi+7];		
				real i = obsVector[mi+8];
				real invI;
				if (i != 0) {
					invI = 1./i;
				} else {
					invI = 0.;
				}
				real j = obsVector[mi+9];
				real rhoBar = rhoBase*exp(-rhoInvScaleHeight*j);
				int ii = (int)((i - iMin)*DIrecip);
				if ((iIndex < ii-1) or (iIndex > ii+2)) continue;
				int jj = (int)((j - jMin)*DJrecip);
				if ((jIndex < jj-1) or (jIndex > jj+2)) continue;
				real ibasis = 0;
				real jbasis = 0;
				real idbasis = 0;
				real jdbasis = 0;
				
				if(w1) {
					ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 2);
					Astate[varDim*iDim*jIndex +varDim*iIndex] 
						+= yhat[m] * ibasis * jbasis * w1 * invError * invI * rhoBar;
				}
				if(w2 or w3) {
					ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 4);
					jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 4);
					idbasis = DBasis(iIndex, i, iDim-1, iMin, DI, DIrecip, 4);
					jdbasis = DBasis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 4);
					Astate[varDim*iDim*jIndex +varDim*iIndex + 1] 
						+= yhat[m] * ibasis * (-jdbasis) * w2 * 1e3 * invError * invI * rhoBar;
					Astate[varDim*iDim*jIndex +varDim*iIndex + 1]
						+= yhat[m] * idbasis * jbasis * w3 * invError * invI * rhoBar;
				}
				if (w4 or w5 or w6) {
					ibasis = Basis(iIndex, i, iDim-1, iMin, DI, DIrecip, 2);
					jbasis = Basis(jIndex, j, jDim-1, jMin, DJ, DJrecip, 2);
					Astate[varDim*iDim*jIndex +varDim*iIndex + 2] 
						+= yhat[m] * ibasis * jbasis * w4 * invError * rhoBar;
					Astate[varDim*iDim*jIndex +varDim*iIndex + 3] 
						+= yhat[m] * ibasis * jbasis * w5 * invError * rhoBar;
					Astate[varDim*iDim*jIndex +varDim*iIndex + 4]
						+= yhat[m] * ibasis * jbasis * w6 * invError;
				}
				
			}
		}
	}
	
}

bool CostFunctionRZ_CPU::SAtransform(real* Bstate, real* Astate)
{
	int k;
	for (int var = 0; var < varDim; var++) {
		real* jB = new real[jDim];
		real* x = new real[jDim];
		for (int iIndex = 0; iIndex < iDim; iIndex++) {
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				jB[jIndex] = Bstate[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			// Solve for A's using compact storage
			real sum = 0;
			for (int j = 0; j < jDim; j++) {
				for (sum=jB[j], k=-1;k>=-3;k--) {
					if ((j+k >= 0) and ((j*4-k) >= 0))
						sum -= jL[jDim*4*var + j*4-k]*x[j+k];
				}
				x[j] = sum/jL[jDim*4*var + j*4];
			}	
			for (int j=jDim-1;j>=0;j--) {
				for (sum=x[j], k=1;k<=3;k++) {
					if ((j+k < jDim) and (((j+k)*4+k) < jDim*4))
						sum -= jL[jDim*4*var + (j+k)*4+k]*x[j+k];
				}
				x[j] = sum/jL[jDim*4*var + j*4];
			}
			
			for (int jIndex = 0; jIndex < jDim; jIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = x[jIndex]; 
			}
		}
		delete[] jB;
		delete[] x;
		
		real* iB = new real[iDim];
		x = new real[iDim];
		for (int jIndex = 0; jIndex < jDim; jIndex++) {
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				iB[iIndex] = Astate[varDim*iDim*jIndex +varDim*iIndex + var];
			}
			// Solve for A's using compact storage
			real sum = 0;
			for (int i = 0; i < iDim; i++) {
				for (sum=iB[i], k=-1;k>=-3;k--) {
					if ((i+k >= 0) and ((i*4-k) >= 0))
						sum -= iL[iDim*4*var + i*4-k]*x[i+k];
				}
				x[i] = sum/iL[iDim*4*var + i*4];
			}	
			for (int i=iDim-1;i>=0;i--) {
				for (sum=x[i], k=1;k<=3;k++) {
					if ((i+k < iDim) and (((i+k)*4+k) < iDim*4))
						sum -= iL[iDim*4*var + (i+k)*4+k]*x[i+k];
				}
				x[i] = sum/iL[iDim*4*var + i*4];
			}
			
			for (int iIndex = 0; iIndex < iDim; iIndex++) {
				Astate[varDim*iDim*jIndex +varDim*iIndex + var] = x[iIndex]; 
			}
			
		}
		delete[] iB;
		delete[] x;
	}

	return true;
}		

bool CostFunctionRZ_CPU::setupSplines()
{
	
	// Do the spline via a Cholesky decomposition
	// and manipulate the DC Filter
	real Pi = 3.141592653589793;
	real** P = new real*[iDim];
	real* p = new real[iDim];
	iL = new real[varDim*iDim*4];
	real* iBL = new real[varDim*iDim*4];
	for (int i = 0; i < iDim; i++) {
		P[i] = new real[iDim];
		p[i] = 0.;
	}
		
	for (int i = 0; i < varDim*iDim*4; i++) {
		iL[i] = 0;
		iBL[i] = 0;
	}
	int BC;
	real cutoff_wl = 4;
	real eq = pow( (cutoff_wl/(2*Pi)) , 6);
	for (int var = 0; var < varDim; var++) {
		if (var <= 1) {
			BC = 4;
		} else {
			BC = 2;
		}
		
		for (int i = 0; i < iDim; i++) {
			for (int j = 0; j < iDim; j++) {
				P[i][j] = 0;

			}
		}
				
		for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
			for (int imu = -1; imu <= 1; imu += 2) {
				real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
				int ii = (int)((i - iMin)*DIrecip);
				for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
					if ((iNode < 0) or (iNode >= iDim)) continue;
					real pm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, BC);
					real qm = DDDBasis(iNode, i, iDim-1, iMin, DI, DIrecip, BC);
					real pn, qn;
					P[iNode][iNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((iNode+1) < iDim) {
						pn = Basis(iNode+1, i, iDim-1, iMin, DI, DIrecip, BC);
						qn = DDDBasis(iNode+1, i, iDim-1, iMin, DI, DIrecip, BC);
						P[iNode][iNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+1][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((iNode+2) < iDim) {
						pn = Basis(iNode+2, i, iDim-1, iMin, DI, DIrecip, BC);
						qn = DDDBasis(iNode+2, i, iDim-1, iMin, DI, DIrecip, BC);
						P[iNode][iNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+2][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((iNode+3) < iDim) {
						pn = Basis(iNode+3, i, iDim-1, iMin, DI, DIrecip, BC);
						qn = DDDBasis(iNode+3, i, iDim-1, iMin, DI, DIrecip, BC);
						P[iNode][iNode+3] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[iNode+3][iNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
				}
			}
		}
				
				
		// Cholesky decomp of P+Q
		for (int i=0;i<iDim;i++) {
			for (int j=i;j<iDim;j++) {
				double sum=P[i][j];
				for (int k=i-1;k>=0;k--) {
					sum -= P[i][k]*P[j][k];
				}
				if (i == j) {
					if (sum <= 0.0) { 
						std::cout << "cholesky failed at i,j sum\n";
						break;
					} else {
						p[i] = sqrt(sum);
					}
				} else {
					P[j][i]=sum/p[i];
					if (p[i] == 0.) { 
						std::cout << "Problem! " << i << "\t" << j << "\n";
					}
				}
			}
		}
		
		for (int i = 0; i < iDim; i++) {
			iL[iDim*4*var + i*4] = p[i];
			//cout << iL[iDim*4*var + i*4] << " ";
			for (int n=1;n<4;n++) {
				if ((i-n) >= 0) {
					iL[iDim*4*var + i*4+n] = P[i][i-n];
				}
				//cout << iL[iDim*4*var + i*4 + n] << " ";
			} //cout << endl;
		} //cout << endl;
	}
	
	// Clear the temporary memory and reallocate for jDim
	for (int i = 0; i < iDim; i++) {
		delete[] P[i];
	}
	delete[] P;
	delete[] p;
	
	P = new real*[jDim];
	p = new real[jDim];
	jL = new real[varDim*jDim*4];
	for (int j = 0; j < jDim; j++) {
		P[j] = new real[jDim];
		p[j] = 0.;
	}
	
	for (int j = 0; j < varDim*jDim*4; j++) {
		jL[j] = 0;
	}
	
	cutoff_wl = 4;
	eq = pow( (cutoff_wl/(2*Pi)) , 6);
	for (int var = 0; var < varDim; var++) {
		if (var == 1) {
			BC = 4;
		} else {
			BC = 2;
		}
		for (int i = 0; i < jDim; i++) {
			for (int j = 0; j < jDim; j++) {
				P[i][j] = 0;
			}
		}		
		for (int jIndex = 0; jIndex < (jDim-1); jIndex++) {
			for (int jmu = -1; jmu <= 1; jmu += 2) {
				real j = jMin + DJ * (jIndex + (0.5*sqrt(1./3.) * jmu + 0.5));
				int jj = (int)((j - jMin)*DJrecip);
				for (int jNode = jj-1; jNode <= jj+2; ++jNode) {
					if ((jNode < 0) or (jNode >= jDim)) continue;
					real pm = Basis(jNode, j, jDim-1, jMin, DJ, DJrecip, BC);
					real qm = DDDBasis(jNode, j, jDim-1, jMin, DJ, DJrecip, BC);
					real pn, qn;
					P[jNode][jNode] += 0.5 * ((pm * pm) + eq * (qm * qm));
					if ((jNode+1) < jDim) {
						pn = Basis(jNode+1, j, jDim-1, jMin, DJ, DJrecip, BC);
						qn = DDDBasis(jNode+1, j, jDim-1, jMin, DJ, DJrecip, BC);
						P[jNode][jNode+1] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+1][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));
					}
					if ((jNode+2) < jDim) {
						pn = Basis(jNode+2, j, jDim-1, jMin, DJ, DJrecip, BC);
						qn = DDDBasis(jNode+2, j, jDim-1, jMin, DJ, DJrecip, BC);
						P[jNode][jNode+2] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+2][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));						
					}
					if ((jNode+3) < jDim) {
						pn = Basis(jNode+3, j, jDim-1, jMin, DJ, DJrecip, BC);
						qn = DDDBasis(jNode+3, j, jDim-1, jMin, DJ, DJrecip, BC);
						P[jNode][jNode+3] += 0.5 * ((pm * pn) + eq * (qm * qn));
						P[jNode+3][jNode] += 0.5 * ((pm * pn) + eq * (qm * qn));						
					}
				}
			}
		}
		
		// Cholesky decomp
		for (int i=0;i<jDim;i++) {
			for (int j=i;j<jDim;j++) {
				double sum=P[i][j];
				for (int k=i-1;k>=0;k--) {
					sum -= P[i][k]*P[j][k];
				}
				if (i == j) {
					if (sum <= 0.0) { 
						std::cout << "cholesky failed at i,j sum\n";
						break;
					} else {
						p[i] = sqrt(sum);
					}
				} else {
					P[j][i]=sum/p[i];
					if (p[i] == 0.) { 
						std::cout << "Problem! " << i << "\t" << j << "\n";
					}
				}
			}
		}
		
		for (int j = 0; j < jDim; j++) {
			jL[jDim*4*var + j*4] = p[j];
			//cout << jL[jDim*4*var + j*4] << " ";
			for (int n=1;n<4;n++) {
				if ((j-n) >= 0) {
					jL[jDim*4*var + j*4+n] = P[j][j-n];
				}
				//cout << jL[jDim*4*var + j*4 + n] << " ";
			} //cout << endl;
		} //cout << endl;
	}
	
	// Clear the temporary memory
	for (int j = 0; j < jDim; j++) {
		delete[] P[j];
	}
	delete[] P;
	delete[] p;
	
	/*p = new real[iDim];
	//Analytic Test set
	for (int i = 0; i < iDim; i++) {
		p[i] = 0;
	}
	int var = 2;
	for (int iIndex = 0; iIndex < (iDim-1); iIndex++) {
		for (int imu = -1; imu <= 1; imu += 2) {
			real i = iMin + DI * (iIndex + (0.5*sqrt(1./3.) * imu + 0.5));
			int ii = (int)((i - iMin)*DIrecip);
			for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
				if ((iNode < 0) or (iNode >= iDim)) continue;
				real pm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 2);
				p[iNode] += 0.5*pm*10.;
			}
		}
	}
	real* x = new real[iDim];
	real* y = new real[iDim];
	// Solve for A's
	real sum = 0;
	int k;
	/*for (int i = 0; i < iDim; i++) {
		for (sum=p[i], k=i-1;k>=0;k--)
			sum -= P[i][k]*x[k];
		x[i] = sum/iL[iDim*4*var + i*4];
	}	
	for (int i=iDim-1;i>=0;i--) {
		for (sum=x[i], k=i+1;k<iDim;k++)
			sum -= P[k][i]*x[k];
		x[i] = sum/iL[iDim*4*var + i*4];
	} */
	 
	/* Solve for A's using compact storage
	sum = 0;
	for (int i = 0; i < iDim; i++) {
		for (sum=p[i], k=-1;k>=-3;k--) {
			if ((i+k >= 0) and ((i*4-k) >= 0))
				sum -= iL[iDim*4*var + i*4-k]*y[i+k];
		}
		y[i] = sum/iL[iDim*4*var + i*4];
	}	
	for (int i=iDim-1;i>=0;i--) {
		for (sum=y[i], k=1;k<=3;k++) {
			if ((i+k < iDim) and (((i+k)*4+k) < iDim*4))
				sum -= iL[iDim*4*var + (i+k)*4+k]*y[i+k];
		}
		y[i] = sum/iL[iDim*4*var + i*4];
	}
	//for (int i = 0; i < iDim; i++) {
	//	cout << i << "\t" << x[i] << "\t" << y[i] << endl;
	//} cout << endl;
	
	real si;
	for (int iIndex = 0; iIndex < iDim; iIndex++) {
		real i = iMin + DI * iIndex;
		int ii = (int)((i - iMin)*DIrecip);
		si = 0;
		for (int iNode = ii-1; iNode <= ii+2; ++iNode) {
			if ((iNode < 0) or (iNode >= iDim)) continue;
			real pm = Basis(iNode, i, iDim-1, iMin, DI, DIrecip, 2);
			si += y[iNode]*pm; // + sin(2*i*2*Pi/iDim) + sin(4*i*2*Pi/iDim) + sin(8*i*2*Pi/iDim) + sin(16*i*2*Pi/iDim);
		}
		//real analytic = sin(i*2*Pi/iDim) + sin(2*i*2*Pi/iDim) + sin(4*i*2*Pi/iDim) + sin(8*i*2*Pi/iDim) + sin(16*i*2*Pi/iDim); 
		real analytic = 10.;
		cout << i << "\t" << p[iIndex] << "\t" << y[iIndex] << "\t" << si << "\t" << analytic << "\t" << si - analytic << endl;
	} */
	
	return true;
	
}

// Basis Functions
float CostFunctionRZ_CPU::Basis(int m, float x, int M, float xmin, 
								float DX, float DXrecip, int C)
{
	
	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
	float ONESIXTH = 0.16666666666666666666666666667;
	float BC[9][4] = 
	//	0		1		M-1		M
	{{	-4,		-1,		-1,		-4 },
		{	0,		1,		1,		0 },
		{	2,		-1,		-1,		2 },
		{   -4,     -1,     1,      0 },
		{   -4,     -1,     -1,     2 },
		{   0,      1,      -1,     -4 },
		{   0,      1,      -1,     2 },
		{   2,      -1,     -1,     -4 },
		{   2,      -1,     1,      0 }};
	
	if (z < 2.0)
	{
		z = 2 - z;
		b = (z*z*z) * ONESIXTH;
		z -= 1.0;
		if (z > 0)
			b -= (z*z*z) * 4 * ONESIXTH;
	}
	
	// Boundary conditions, if any, are an additional addend.
	if (C >= 0) {
		if (m == 0 || m == 1) {
			float l = 0;
			xm = xmin + (-1 * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			
			if (z < 2.0)
			{
				z = 2 - z;
				l = (z*z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					l -= (z*z*z) * 4 * ONESIXTH;
			}
			b += BC[C][m] * l;
		} else if (m == M-1 || m == M) {
			float r = 0;
			xm = xmin + ((M+1) * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			
			if (z < 2.0)
			{
				z = 2 - z;
				r = (z*z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					r -= (z*z*z) * 4 * ONESIXTH;
			}
			b += BC[C][m+3-M] * r;
		}
	}
	return b;
}

float CostFunctionRZ_CPU::DBasis(int m, float x, int M, float xmin, 
								 float DX, float DXrecip, int C)
{
	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
	float ONESIXTH = 0.16666666666666666666666666667;
	float BC[9][4] = 
	//	0		1		M-1		M
	{{	-4,		-1,		-1,		-4 },
		{	0,		1,		1,		0 },
		{	2,		-1,		-1,		2 },
		{   -4,     -1,     1,      0 },
		{   -4,     -1,     -1,     2 },
		{   0,      1,      -1,     -4 },
		{   0,      1,      -1,     2 },
		{   2,      -1,     -1,     -4 },
		{   2,      -1,     1,      0 }};
	
	if (z < 2.0)
	{
		z = 2.0 - z;
		b = (z*z) * ONESIXTH;
		z -= 1.0;
		if (z > 0)
			b -= (z*z) * 4 * ONESIXTH;
		b *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
	}
	
	// Boundary conditions, if any, are an additional addend.
	if (C >= 0) {
		if (m == 0 || m == 1) {
			float l = 0;
			xm = xmin + (-1 * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			
			if (z < 2.0)
			{
				z = 2 - z;
				l = (z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					l -= (z*z) * 4 * ONESIXTH;
				l *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;
			}
			
			b += BC[C][m] * l;
		} else if (m == M-1 || m == M) {
			float r = 0;
			xm = xmin + ((M+1) * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			
			if (z < 2.0)
			{
				z = 2 - z;
				r = (z*z) * ONESIXTH;
				z -= 1.0;
				if (z > 0)
					r -= (z*z) * 4 * ONESIXTH;
				r *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * DXrecip;	
			}
			b += BC[C][m+3-M] * r;
		}
	}
	return b;
}


float CostFunctionRZ_CPU::DDBasis(int m, float x, int M, float xmin, 
								  float DX, float DXrecip, int C)
{
	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
	float BC[9][4] = 
	//	0		1		M-1		M
	{{	-4,		-1,		-1,		-4 },
		{	0,		1,		1,		0 },
		{	2,		-1,		-1,		2 },
		{   -4,     -1,     1,      0 },
		{   -4,     -1,     -1,     2 },
		{   0,      1,      -1,     -4 },
		{   0,      1,      -1,     2 },
		{   2,      -1,     -1,     -4 },
		{   2,      -1,     1,      0 }};
	
	if (z < 2.0)
	{
		z = 2.0 - z;
		b = z;
		z -= 1.0;
		if (z > 0)
			b -= z * 4;
		b *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
	}
	
	// Boundary conditions, if any, are an additional addend.
	if (C >= 0) {
		if (m == 0 || m == 1) {
			float l = 0;
			xm = xmin + (-1 * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			if (z < 2.0)
			{
				z = 2.0 - z;
				l = z;
				z -= 1.0;
				if (z > 0)
					l -= z * 4;
				l *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			}
			b += BC[C][m] * l;
		} else if (m == M-1 || m == M) {
			float r = 0;
			xm = xmin + ((M+1) * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			if (z < 2.0)
			{
				z = 2.0 - z;
				r = z;
				z -= 1.0;
				if (z > 0)
					r -= z * 4;
				r *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			}
			b += BC[C][m+3-M] * r;
		}
	}
	return b;
}

float CostFunctionRZ_CPU::DDDBasis(int m, float x, int M, float xmin, 
								   float DX, float DXrecip, int C)
{
	float b = 0;
	float xm = xmin + (m * DX);
	float delta = (x - xm) * DXrecip;
	float z = fabsf(delta);
	float BC[9][4] = 
	//	0		1		M-1		M
	{{	-4,		-1,		-1,		-4 },
		{	0,		1,		1,		0 },
		{	2,		-1,		-1,		2 },
		{   -4,     -1,     1,      0 },
		{   -4,     -1,     -1,     2 },
		{   0,      1,      -1,     -4 },
		{   0,      1,      -1,     2 },
		{   2,      -1,     -1,     -4 },
		{   2,      -1,     1,      0 }};
	
	if ((z > 1.0) and (z < 2.0)) {
		b = 1;
	} else if (z < 1.0) {
		b = -3.;
	}
	b *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
	
	// Boundary conditions, if any, are an additional addend.
	if (C >= 0) {
		if (m == 0 || m == 1) {
			float l = 0;
			xm = xmin + (-1 * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			if ((z > 1.0) and (z < 2.0)) {
				l = 1;
			} else if (z < 1.0) {
				l = -3.;
			}
			l *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			
			b += BC[C][m] * l;
		} else if (m == M-1 || m == M) {
			float r = 0;
			xm = xmin + ((M+1) * DX);
			delta = (x - xm) * DXrecip;
			z = fabsf(delta);
			if ((z > 1.0) and (z < 2.0)) {
				r = 1;
			} else if (z < 1.0) {
				r = -3.;
			}
			r *= ((delta > 0) ? -1.0 : 1.0) * DXrecip;
			b += BC[C][m+3-M] * r;
		}
	}
	return b;
}
