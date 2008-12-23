/*
 *  CostFunctionRZ.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionRZ.h"
#include <cmath>

// Global CUDA declarations
extern "C" {
	int cstGPU_init (float* obs_h, int mObs);
	void loadSplineCoeffs_GPU(float* coeffHost, int numCoeffs);
	void HCq_GPU(int mObs, float rmax, float rmin, float zmax, float zmin, float* HCq_h, int pState, int zState);
}

/* First, compute HCq and Hxb via recursive filter on q, followed by loop over M observations
 yielding 2 y-estimate M x 1 vectors. Solve for innovation vector d by subtracting from y.
 For function, subtract from HCq and take inner product (with R-1), add to i.p. of q vector.
 For gradient, loop over M again, this time summing product of operator*basis*oberror-1*y-hat
 and operator*basis*oberror1*d for each coefficient N, yielding an N x 1 vector. Filter this
 vector and the innovation vector. Add 1+vector1+vector2 to get the gradient. Treat the penalty
 constraints as observations. */

CostFunctionRZ::CostFunctionRZ(const int& numObs, const int& stateSize)
	: CostFunction(numObs, stateSize)
{
	/* Initialize background errors
	bgError[0] = 2.;
	bgError[1] = 2.e7;
	bgError[2] = 20.;
	bgError[3] = 2.;
	bgError[4] = 0.1; */

	bgError[0] = 0.1;
	bgError[1] = 0.2;
	bgError[2] = 0.1;
	bgError[3] = 0.1;
	bgError[4] = 0.1;
}

CostFunctionRZ::~CostFunctionRZ()
{
}

void CostFunctionRZ::finalize()
{	
	delete filterR;
	delete filterZ;
	delete[] currState;
	delete[] currGradient;
	delete[] innovation;
	delete[] tempState;
	delete[] tempGradient;
	delete[] HCq;
	delete[] xt;
	delete[] df;
	delete[] stateMod;
	delete[] HT;
	delete[] zHT;
	delete[] rCTHT;
	delete[] zCTHT;
	unsigned int heightSize = bgheights->size();
	for (unsigned int z = 0; z < heightSize; z++) {
		delete[] rzCTHT[z];
		delete[] rzSTHT[z];
	}
	delete[] rzCTHT;
	delete[] rzSTHT;
	for (unsigned int var = 0; var < numVars; var++) {
		for (unsigned int z = 0; z < zState; z++) {
			delete[] CTHTd[var][z];
		}
		delete[] CTHTd[var];
	}
	delete[] CTHTd;
		
	for (unsigned int z = 0; z < zState; z++) {
		delete[] field[z];
	}
	delete[] field;
	delete[] fieldR;
	delete[] fieldZ;
	delete[] zBuffer;
	
#ifdef CSTGPU
	// GPU variables
	delete[] coeffHost;
#endif
	
}

void CostFunctionRZ::initialize(SplineD* vecs, SplineD* scalars, SplineD* ctrls, SplineD* zs, SplineD* zpsi, 
								vector<double>* bgr, vector<double>* bgz, vector<double>** bgf, 
								unsigned int* ctrlR, vector<double>* RX, Observation* obs)
{

	// Initialize number of control variables -- one less than physical variables
	numVars = 5;
	
	// Assign local object pointers
	vecSpline = vecs;
	scalarSpline = scalars;
	bgradii = bgr;
	bgheights = bgz;
	bgFields = bgf;
	ctrlSpline = ctrls;
	RXform = RX;
    zSpline = zs;
	zSplinePsi = zpsi;
	obsVector = obs;
	pState = vecSpline[0].nNodes();
	zState = zSpline->nNodes();
	RState = ctrlR;
		
	// Set up the recursive filter
	filterR = new RecursiveFilter(4,8);
	filterZ = new RecursiveFilter(4,4);
	
	// Allocate memory for the needed arrays
	currState = new double[nState];
	currGradient = new double[nState];
	tempState = new double[nState];
	tempGradient = new double[nState];
	xt = new double[nState];
	df = new double[nState];
	stateMod = new double[nState];

	innovation = new double[mObs];

	HT = new double[pState];
	zHT = new double[zState];
	
	unsigned int heightSize = bgheights->size();
	maxRadius = 0;
	for (unsigned int z = 0; z < heightSize; z++)
		if (RState[z] > maxRadius) maxRadius = RState[z];
	
	rCTHT = new double[maxRadius];
	zCTHT= new double[zState];
	rzCTHT= new double*[heightSize];
	CTHTd = new double**[numVars];
	rzSTHT = new double*[heightSize];
	for (unsigned int var = 0; var < numVars; var++) {
		CTHTd[var] = new double*[zState];
		for (unsigned int z = 0; z < zState; z++) {
			CTHTd[var][z] = new double[maxRadius];
		}
	}
	
	for (unsigned int z = 0; z < heightSize; z++) {
		rzCTHT[z] = new double[maxRadius];
		rzSTHT[z] = new double[maxRadius];
	}
	
	fieldR = new double[maxRadius];
	fieldZ = new double[zState];
	field = new double*[zState];
	for (unsigned int z = 0; z < zState; z++) {
		field[z] = new double[maxRadius];
	}
	zBuffer = new double[heightSize];
	
#ifdef CSTGPU
	// GPU variables
	HCq = new float[mObs];
	coeffHost = new float[numVars*pState*zState];
	// Load the Observations
	float* obs_h = new float[mObs*9];
	for (int m = 0; m < mObs; m++) {
		int n = m*9;
		for (int var = 0; var < 6; var++) {
			obs_h[n+var] = obs[m].getWeight(var);
		}
		obs_h[n+6] = obs[m].getRadius();
		obs_h[n+7] = obs[m].getAltitude();
		obs_h[n+8] = obs[m].getInverseError();
	}
	cstGPU_init(obs_h, mObs);
	
#else
	HCq = new float[mObs];
	//HCq = new double[mObs];	
#endif
	
	initState();
}	

void CostFunctionRZ::initState() {

	// Clear the state vector
	for (int n = 0; n < nState; n++) {
		currState[n] = 0.0;
		currGradient[n] = 0.0;
		tempState[n] = 0.0;
		tempGradient[n] = 0.0;
		xt[n] = 0.0;
		df[n] = 0.0;
		stateMod[n] = 0.0;
	}
	
	// Fill the innovation vector
	cout << "Initializing innovation vector..." << endl;
	for (int m = 0; m < mObs; m++) {
		HCq[m] = 0.0;
		innovation[m] = obsVector[m].getOb();
	}
	
	// Compute the scaling factors
	for (unsigned int var = 0; var < numVars; var++) {
		varScale[var] = 0;
		for (unsigned int z = 0; z < bgheights->size(); z++) {
			for (unsigned int r = 0; r < bgradii->size(); r++) {
				varScale[var] += bgFields[z][var].at(r)*bgFields[z][var].at(r);
			}
		}
		varScale[var] = sqrt(varScale[var]/(bgradii->size()*bgheights->size()));
		if ((var == 1) and (varScale[var] == 0)) varScale[var] = 240;
		cout << "Variable RMS scale factor = " << var << "\t" << varScale[var] << endl;
	}
	
	for (unsigned int var = 0; var < numVars; var++) {
		SplineD* bgSpline;
		if (var > 1) {
			bgSpline = scalarSpline;
		} else {
			bgSpline = vecSpline;
		}		
		for (unsigned int z = 0; z < zState; z++) {
			bgSpline[z].solveGQ(&bgFields[z][var].front());
		}
		for (int m = 0; m < mObs; m++) {
			if (var == 0) {
				double weight =  obsVector[m].getWeight(var);
				if (!weight) continue;
				for (unsigned int z = 0; z < zState; z++) {
					zSpline->setCoefficient(z, bgSpline[z].evaluate(obsVector[m].getRadius()));
				}
				double yhat =  zSpline->evaluate(obsVector[m].getAltitude());
				innovation[m] -= yhat*weight;
			} else if (var == 1) {
				double weight =  obsVector[m].getWeight(var);
				if (weight) {
					for (unsigned int z = 0; z < zState; z++) {								
						zSplinePsi->setCoefficient(z, bgSpline[z].evaluate(obsVector[m].getRadius()));
					}
					double yhat =  1e6 * -(zSplinePsi->slope(obsVector[m].getAltitude()) / (obsVector[m].getRadius() * 1000));
					innovation[m] -= yhat*weight;
				}
				weight =  obsVector[m].getWeight(var+1);
				if (!weight) continue;
				for (unsigned int z = 0; z < zState; z++) {								
					zSplinePsi->setCoefficient(z, (bgSpline[z].slope(obsVector[m].getRadius()) / 1000));
				}
				double yhat =  1e6 * zSplinePsi->evaluate(obsVector[m].getAltitude()) / (obsVector[m].getRadius() * 1000);
				innovation[m] -= yhat*weight;
				
			} else {
				double weight =  obsVector[m].getWeight(var+1);
				if (!weight) continue;
				for (unsigned int z = 0; z < zState; z++) {								
					zSpline->setCoefficient(z, bgSpline[z].evaluate(obsVector[m].getRadius()));
				}
				double yhat =  zSpline->evaluate(obsVector[m].getAltitude());
				innovation[m] -= yhat*weight;
			}				
		}
	}
	
	// Calculate HTd	
	for (unsigned int var = 0; var < numVars; var++) {
		SplineD* bgSpline;
		if (var > 1) {
			bgSpline = scalarSpline;
		} else {
			bgSpline = vecSpline;
		}		
		for (int z = 0; z < zSpline->nNodes(); z++) {
			for (int p = 0; p < vecSpline[z].nNodes(); p++) {
				HT[p] = 0;
				double HTsum = 0;
				#pragma omp parallel for reduction(+:HTsum)
				for (int m = 0; m < mObs; m++) {
					// Sum over obs this time
					// Multiply state by H weights
					if (var == 0) {
						double weight = obsVector[m].getWeight(var);
						if (!weight) continue;
						double rbasis = bgSpline[z].getBasis(p, obsVector[m].getRadius());
						if (!rbasis) continue;
						double zbasis = zSpline->getBasis(z, obsVector[m].getAltitude());
						if (!zbasis) continue;
						double invError = obsVector[m].getInverseError();
						double nodeWeight = rbasis*zbasis*weight*invError;
						HTsum += innovation[m]*(nodeWeight);
					} else if (var == 1) {
						double weight = obsVector[m].getWeight(var);
						if (weight) {
							double rbasis = bgSpline[z].getBasis(p, obsVector[m].getRadius());
							if (rbasis) {
								double zbasis = zSplinePsi->getDBasis(z, obsVector[m].getAltitude());
								if (zbasis) {
									double invError = obsVector[m].getInverseError();
									double nodeWeight = 1e6 * -(rbasis*zbasis*weight*invError/(obsVector[m].getRadius() * 1000));
									HTsum += innovation[m]*(nodeWeight);
								}
							}
						}
						
						weight = obsVector[m].getWeight(var+1);
						if (!weight) continue;
						double rbasis = bgSpline[z].getDBasis(p, obsVector[m].getRadius())  / 1000;
						if (!rbasis) continue;
						double zbasis = zSplinePsi->getBasis(z, obsVector[m].getAltitude());
						if (!zbasis) continue;
						double invError = obsVector[m].getInverseError();
						double nodeWeight = 1e6 * rbasis*zbasis*weight*invError/(obsVector[m].getRadius() * 1000);
						HTsum += innovation[m]*(nodeWeight);
					} else {
						double weight = obsVector[m].getWeight(var+1);
						if (!weight) continue;
						double rbasis = bgSpline[z].getBasis(p, obsVector[m].getRadius());
						if (!rbasis) continue;
						double zbasis = zSpline->getBasis(z, obsVector[m].getAltitude());
						if (!zbasis) continue;
						double invError = obsVector[m].getInverseError();
						double nodeWeight = rbasis*zbasis*weight*invError;
						HTsum += innovation[m]*(nodeWeight);
					}
				}
				HT[p] = HTsum;
			}
			
			// S^T
			const double* STHT = bgSpline[z].solveInverseGQ(&HT[0]);

			// P^T
			for (unsigned int R = 0; R < RState[z]; R++) {
				rCTHT[R] = 0;
				for (unsigned int r = 0; r < bgradii->size(); r++) {
					double potRad = RXform[z].at(r);
					rCTHT[R] += STHT[r]*ctrlSpline[z].getBasis(R, potRad);
				}
				rCTHT[R] *= (bgError[var] * varScale[var]);
			}
			
			// FR
			filterR->filterArray(rCTHT, RState[z]);
			for (unsigned int R = 0; R < RState[z]; R++) {
				if ((var == 1) and (z == 0)) {
					// Force Psi delta to zero
					rCTHT[R]= 0;
				}			
				rzCTHT[z][R] = rCTHT[R];
			}
		}
		
		for (unsigned int R = 0; R < maxRadius; R++) {
			for (unsigned int z = 0; z < zState; z++) {
				// Pad the field with zeroes if it is outside the domain
				if (R < RState[z]) {
					zCTHT[z] = rzCTHT[z][R];
				} else {
					zCTHT[z] = 0;
				}
			}
			
			// FZ
			filterZ->filterArray(zCTHT, zState);
			for (unsigned int z = 0; z < zState; z++) {
				CTHTd[var][z][R] = zCTHT[z];
			}
		}
	}
	
}	

double CostFunctionRZ::funcValue(double* state)
{
	// Update the Y hat vector
#ifdef CSTGPU
	updateHCq_GPU(state);
#else
	updateHCq(state);
#endif
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
		obIP += (HCq[m]-innovation[m])*(obsVector[m].getInverseError())*(HCq[m]-innovation[m]);
		//cout << m << "\t" << HCq[m] << "\t" << innovation[m] << endl;
	}
	double J = 0.5*qIP + 0.5*obIP;
	return J;
	
}

void CostFunctionRZ::funcGradient(double* state, double* gradient)
{
	
	// Update the Y hat vector
#ifdef CSTGPU
	updateHCq_GPU(state);
#else
	updateHCq(state);
#endif

	// Calculate HTd	
	for (unsigned int var = 0; var < numVars; var++) {
		SplineD* bgSpline;
		if (var > 1) {
			bgSpline = scalarSpline;
		} else {
			bgSpline = vecSpline;
		}		
		for (int z = 0; z < zSpline->nNodes(); z++) {
			for (int p = 0; p < vecSpline[z].nNodes(); p++) {
				HT[p] = 0;
				double HTsum = 0;
				#pragma omp parallel for reduction(+:HTsum)
				for (int m = 0; m < mObs; m++) {
					// Sum over obs this time
					// Multiply state by H weights
					if (var == 0) {
						double weight = obsVector[m].getWeight(var);
						if (!weight) continue;
						double rbasis = bgSpline[z].getBasis(p, obsVector[m].getRadius());
						if (!rbasis) continue;
						double zbasis = zSpline->getBasis(z, obsVector[m].getAltitude());
						if (!zbasis) continue;
						double invError = obsVector[m].getInverseError();
						double nodeWeight = rbasis*zbasis*weight*invError;
						HTsum += HCq[m]*(nodeWeight);
					} else if (var == 1) {
						double weight = obsVector[m].getWeight(var);
						if (weight) {
							double rbasis = bgSpline[z].getBasis(p, obsVector[m].getRadius());
							if (rbasis) {
								double zbasis = zSplinePsi->getDBasis(z, obsVector[m].getAltitude());
								if (zbasis) {
									double invError = obsVector[m].getInverseError();
									double nodeWeight = 1e6 * -(rbasis*zbasis*weight*invError)/(obsVector[m].getRadius() * 1000);
									HTsum += HCq[m]*(nodeWeight);
								}
							}
						}
						
						weight = obsVector[m].getWeight(var+1);
						if (!weight) continue;
						double rbasis = bgSpline[z].getDBasis(p, obsVector[m].getRadius())  / 1000;
						if (!rbasis) continue;
						double zbasis = zSplinePsi->getBasis(z, obsVector[m].getAltitude());
						if (!zbasis) continue;
						double invError = obsVector[m].getInverseError();
						double nodeWeight = 1e6 * (rbasis*zbasis*weight*invError)/(obsVector[m].getRadius() * 1000);
						HTsum += HCq[m]*(nodeWeight);
					} else {
						double weight = obsVector[m].getWeight(var+1);
						if (!weight) continue;
						double rbasis = bgSpline[z].getBasis(p, obsVector[m].getRadius());
						if (!rbasis) continue;
						double zbasis = zSpline->getBasis(z, obsVector[m].getAltitude());
						if (!zbasis) continue;
						double invError = obsVector[m].getInverseError();
						double nodeWeight = rbasis*zbasis*weight*invError;
						HTsum += HCq[m]*(nodeWeight);
					}
				}
				HT[p] = HTsum;
			}
			
			// S^T
			const double* STHT = bgSpline[z].solveInverseGQ(&HT[0]);
			
			// P^T
			for (unsigned int R = 0; R < RState[z]; R++) {
				rCTHT[R] = 0;
				for (unsigned int r = 0; r < bgradii->size(); r++) {
					double potRad = RXform[z].at(r);
					rCTHT[R] += STHT[r]*ctrlSpline[z].getBasis(R, potRad);
				}
				rCTHT[R] *= (bgError[var] * varScale[var]);
			}
			
			// FR
			filterR->filterArray(rCTHT, RState[z]);
			for (unsigned int R = 0; R < RState[z]; R++) {
				if ((var == 1) and (z == 0)) {
					// Force Psi delta to zero
					rCTHT[R]= 0;
				}							
				rzCTHT[z][R] = rCTHT[R];
			}
		}
		
		for (unsigned int R = 0; R < maxRadius; R++) {
			for (unsigned int z = 0; z < zState; z++) {
				// Pad the field with zeroes if it is outside the domain
				if (R < RState[z]) {
					zCTHT[z] = rzCTHT[z][R];
				} else {
					zCTHT[z] = 0;
				}
			}
			
			// FZ
			filterZ->filterArray(zCTHT, zState);
						
			// Assign the state modifier
			unsigned int zi = 0;
			for (unsigned int z = 0; z < zState; z++) {
				unsigned int Ri = R + var*RState[z] + zi;
				stateMod[Ri] = zCTHT[z] - CTHTd[var][z][R];
				// Increment the state array index
				zi += numVars*RState[z];
			}
		}
	}
	
	for (int n = 0; n < nState; n++) {
		gradient[n] = state[n] + stateMod[n];
	}
	
}

void CostFunctionRZ::updateHCq(double* state)
{
		
	// Clear the HCq variable
	for (int m = 0; m < mObs; m++) {
		HCq[m] = 0;
	}
	
	unsigned int radSize = bgradii->size();
	double* Cq = new double[radSize];
	for (unsigned int var = 0; var < numVars; var++) {
		unsigned int zi = 0;
		SplineD* bgSpline;
		if (var > 1) {
			bgSpline = scalarSpline;
		} else {
			bgSpline = vecSpline;
		}		
		for (unsigned int z = 0; z < zState; z++) {
			for (unsigned int R = 0; R < RState[z]; R++) {
				unsigned int Ri = R + var*RState[z] + zi;
				field[z][R] = state[Ri];
			}			
			// Increment the state array index
			zi += numVars*RState[z];
		}
		for (unsigned int R = 0; R < maxRadius; R++) {
			for (unsigned int z = 0; z < zState; z++) {
				// Pad the field with zeroes if it is outside the domain
				if (R < RState[z]) {
					fieldZ[z] = field[z][R];
				} else {
					fieldZ[z] = 0;
				}
			}
			
			// FZ
			filterZ->filterArray(fieldZ, zState);
			for (unsigned int z = 0; z < zState; z++) {
				field[z][R] = fieldZ[z];
			}
			
		}
				
		for (unsigned int z = 0; z < zState; z++) {
			for (unsigned int R = 0; R < RState[z]; R++) {
				if ((var == 1) and (z == 0)) {
					// Force Psi delta to zero
					field[z][R]= 0;
				}			
				fieldR[R] = field[z][R];
			}
			// FR
			filterR->filterArray(fieldR, RState[z]);
			
			// D
			for (unsigned int R = 0; R < RState[z]; R++) {
				double coeff = fieldR[R] * bgError[var] * varScale[var];
				ctrlSpline[z].setCoefficient(R, coeff);
			}

			// P
			for (unsigned int r = 0; r < radSize; r++) {
				double potRad = RXform[z].at(r);
				Cq[r] = ctrlSpline[z].evaluate(potRad);
				//cout << z << "\t" << r << "\t" << var << "\t" << Cq[r] << endl;
			}
			
			// S
			bgSpline[z].solveGQ(Cq);
		}
		
		// H
		for (int m = 0; m < mObs; m++) {
			if (var == 0) {
				double weight =  obsVector[m].getWeight(var);
				if (!weight) continue;
				for (unsigned int z = 0; z < zState; z++) {
					zSpline->setCoefficient(z, bgSpline[z].evaluate(obsVector[m].getRadius()));
				}
				double yhat =  zSpline->evaluate(obsVector[m].getAltitude());
				HCq[m] += yhat*weight;
			} else if (var == 1) {
				double weight =  obsVector[m].getWeight(var);
				if (weight) {
					for (unsigned int z = 0; z < zState; z++) {								
						zSplinePsi->setCoefficient(z, bgSpline[z].evaluate(obsVector[m].getRadius()));
					}
					double yhat =  1e3 * -(zSplinePsi->slope(obsVector[m].getAltitude()) / ( obsVector[m].getRadius() ));
					HCq[m] += yhat*weight;
				}
				weight =  obsVector[m].getWeight(var+1);
				if (!weight) continue;
				for (unsigned int z = 0; z < zState; z++) {								
					zSplinePsi->setCoefficient(z, (bgSpline[z].slope(obsVector[m].getRadius()) ));
				}
				double yhat =  zSplinePsi->evaluate(obsVector[m].getAltitude()) / ( obsVector[m].getRadius() );
				HCq[m] += yhat*weight;
				
			} else {
				double weight =  obsVector[m].getWeight(var+1);
				if (!weight) continue;
				for (unsigned int z = 0; z < zState; z++) {								
					zSpline->setCoefficient(z, bgSpline[z].evaluate(obsVector[m].getRadius()));
				}
				double yhat =  zSpline->evaluate(obsVector[m].getAltitude());
				HCq[m] += yhat*weight;
			}				
			//cout << var << ", HCq: " << m << "\t" << HCq[m] << endl;
		}
	}

	delete[] Cq;
	
}

void CostFunctionRZ::updateHCq_GPU(double* state)
{
		
	unsigned int radSize = bgradii->size();
	double* Cq = new double[radSize];
	for (unsigned int var = 0; var < numVars; var++) {
		unsigned int zi = 0;
		SplineD* bgSpline;
		if (var > 1) {
			bgSpline = scalarSpline;
		} else {
			bgSpline = vecSpline;
		}		
		for (unsigned int z = 0; z < zState; z++) {
			for (unsigned int R = 0; R < RState[z]; R++) {
				unsigned int Ri = R + var*RState[z] + zi;
				field[z][R] = state[Ri];
			}			
			// Increment the state array index
			zi += numVars*RState[z];
		}
		for (unsigned int R = 0; R < maxRadius; R++) {
			for (unsigned int z = 0; z < zState; z++) {
				// Pad the field with zeroes if it is outside the domain
				if (R < RState[z]) {
					fieldZ[z] = field[z][R];
				} else {
					fieldZ[z] = 0;
				}
			}
			
			// FZ
			filterZ->filterArray(fieldZ, zState);
			for (unsigned int z = 0; z < zState; z++) {
				field[z][R] = fieldZ[z];
			}
			
		}
		
		for (unsigned int z = 0; z < zState; z++) {
			for (unsigned int R = 0; R < RState[z]; R++) {
				if ((var == 1) and (z == 0)) {
					// Force Psi delta to zero
					field[z][R]= 0;
				}			
				fieldR[R] = field[z][R];
			}
			// FR
			filterR->filterArray(fieldR, RState[z]);
			
			// D
			for (unsigned int R = 0; R < RState[z]; R++) {
				double coeff = fieldR[R] * bgError[var] * varScale[var];
				ctrlSpline[z].setCoefficient(R, coeff);
			}
			
			// P
			for (unsigned int r = 0; r < radSize; r++) {
				double potRad = RXform[z].at(r);
				Cq[r] = ctrlSpline[z].evaluate(potRad);
				//cout << z << "\t" << r << "\t" << var << "\t" << Cq[r] << endl;
			}
			
			// S
			bgSpline[z].solveGQ(Cq);
		}
		
		// Load the spline coefficients onto the GPU
		for (unsigned int z = 0; z < zState; z++) {
			for (unsigned int p = 0; p < pState; p++) {
				coeffHost[var*zState*pState +z*pState + p] = bgSpline[z].getCoefficient(p);
			}
		}
	}
	
	loadSplineCoeffs_GPU(coeffHost, numVars*pState*zState);
	delete[] Cq;
	
	// H
	HCq_GPU(mObs, bgradii->back(), bgradii->front(), bgheights->back(), bgheights->front(), HCq, pState, zState);
}


void CostFunctionRZ::getCq(double* Cq)
{
	// Apply series of transforms to copy of state vector
	for (unsigned int var = 0; var < numVars; var++) {
		unsigned int zi = 0;
		for (unsigned int z = 0; z < zState; z++) {
			for (unsigned int R = 0; R < RState[z]; R++) {
				unsigned int Ri = R + var*RState[z] + zi;
				field[z][R] = currState[Ri];
			}			
			// Increment the state array index
			zi += numVars*RState[z];
		}
		for (unsigned int R = 0; R < maxRadius; R++) {
			for (unsigned int z = 0; z < zState; z++) {
				// Pad the field with zeroes if it is outside the domain
				if (R < RState[z]) {
					fieldZ[z] = field[z][R];
				} else {
					fieldZ[z] = 0;
				}
			}
			
			// FZ
			filterZ->filterArray(fieldZ, zState);
			for (unsigned int z = 0; z < zState; z++) {
				field[z][R] = fieldZ[z];
			}
			
		}
				
		for (unsigned int z = 0; z < zState; z++) {
			for (unsigned int R = 0; R < RState[z]; R++) {
				if ((var == 1) and (z == 0)) {
					// Force Psi delta to zero
					field[z][R]= 0;
				}			
				fieldR[R] = field[z][R];
			}
			
			// FR
			filterR->filterArray(fieldR, RState[z]);
			
			//D
			for (unsigned int R = 0; R < RState[z]; R++) {
				double coeff = fieldR[R] * bgError[var] * varScale[var];
				ctrlSpline[z].setCoefficient(R, coeff);
			}
					
			// P
			for (unsigned int r = 0; r < bgradii->size(); r++) {
				double potRad = RXform[z].at(r);
				unsigned int i = r + var*bgradii->size() + z*numVars*bgradii->size();
				Cq[i] = ctrlSpline[z].evaluate(potRad);
				//cout << i << "\t" << r << "\t" << var << "\t" << Cq[i] << endl;
			}
			// S handled by VarDriver2d
		}
	}
		
}



