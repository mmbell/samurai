/*
 *  CostFunctionR.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunctionR.h"
#include <cmath>

/* First, compute HCq and Hxb via recursive filter on q, followed by loop over M observations
 yielding 2 y-estimate M x 1 vectors. Solve for innovation vector d by subtracting from y.
 For function, subtract from HCq and take inner product (with R-1), add to i.p. of q vector.
 For gradient, loop over M again, this time summing product of operator*basis*oberror-1*y-hat
 and operator*basis*oberror1*d for each coefficient N, yielding an N x 1 vector. Filter this
 vector and the innovation vector. Add 1+vector1+vector2 to get the gradient. Treat the penalty
 constraints as observations. */

CostFunctionR::CostFunctionR(const int& numObs, const int& stateSize)
	: CostFunction(numObs, stateSize)
{
	// Initialize background errors
	bgError[0] = 2.;
	bgError[1] = 2.;
	bgError[2] = 2.;
	bgError[3] = 20.;
	bgError[4] = 2.;
	bgError[5] = 0.1;
}

CostFunctionR::~CostFunctionR()
{
}

void CostFunctionR::finalize()
{	
	delete filter1d;
	delete[] currState;
	delete[] currGradient;
	delete[] innovation;
	delete[] tempState;
	delete[] tempGradient;
	delete[] HCq;
	delete[] xt;
	delete[] df;
	delete[] stateMod;
	delete[] HTHCq;
	delete[] HTd;
}

void CostFunctionR::initialize(SplineD* bgs, vector<real>* bgr, vector<real>* bgf,
							   SplineD* ctrl, vector<real>* ctrlR,
							   vector<real>* RX, vector<real>* rX, 
							   Observation* obs)
{

	// Initialize number of variables
	numVars = 6;
	
	// Assign local object pointers
	bgSpline = bgs;
	bgradii = bgr;
	bgFields = bgf;
	ctrlSpline = ctrl;
	ctrlRadii = ctrlR;
	RXform = RX;
	rXform = rX;
	obsVector = obs;
	
	// Set up the recursive filter
	filter1d = new RecursiveFilter(4,4);

	// Allocate memory for the needed arrays
	currState = new double[nState];
	currGradient = new double[nState];
	tempState = new double[nState];
	tempGradient = new double[nState];
	xt = new double[nState];
	df = new double[nState];
	stateMod = new double[nState];

	HCq = new real[mObs];
	innovation = new double[mObs];
	int pState = bgSpline->nNodes();
	HTHCq = new real[pState];
	HTd = new real[pState];
	
	initState();
}	

void CostFunctionR::initState() {

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
	
	for (unsigned int var = 0; var < numVars; var++) {
		bgSpline->solve(&bgFields[var].front());
		for (int m = 0; m < mObs; m++) {
			innovation[m] -= bgSpline->evaluate(obsVector[m].getRadius()) * obsVector[m].getWeight(var);
		}
		/* if (fabs(innovation[m]) > 20) 
			cout << m << "\t" << innovation[m] << endl; */
	}
	
}	

double CostFunctionR::funcValue(double* state)
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
	for (int m = 0; m < mObs; m++) {
		obIP += (HCq[m]-innovation[m])*(obsVector[m].getInverseError())*(HCq[m]-innovation[m]);
	}	
	
	double J = 0.5*qIP + 0.5*obIP;
	return J;
	
}

void CostFunctionR::funcGradient(double* state, double* gradient)
{
	
	// Update the Y hat vector
	updateHCq(state);

	// Apply H^T
	int pState = bgSpline->nNodes();	
	unsigned int fieldSize = ctrlRadii->size();
	double* CHTHCq = new double[fieldSize];
	double* CHTd = new double[fieldSize];
	for (unsigned int var = 0; var < numVars; var++) {
		for (int p = 0; p < pState; p++) {
			HTHCq[p] = 0;
			HTd[p] = 0;
			
			for (int m = 0; m < mObs; m++) {
				// Sum over obs this time
				// Multiply state by H weights
				double basis = bgSpline->getBasis(p, obsVector[m].getRadius());
				double weight = obsVector[m].getWeight(var);
				double R = obsVector[m].getInverseError();
				double nodeWeight = basis*weight*R;
				HTHCq[p] += HCq[m]*(nodeWeight);
				HTd[p] += innovation[m]*(nodeWeight);
			}
		}
		
		// P & S & D
		bgSpline->solve(&HTHCq[0]);
		for (unsigned int R = 0; R < fieldSize; R++) {
			double r = rXform->at(R);
			CHTHCq[R] = bgSpline->evaluate(r) * bgError[var];
		}
		
		bgSpline->solve(&HTd[0]);
		for (unsigned int R = 0; R < fieldSize; R++) {
			double r = rXform->at(R);
			CHTd[R] = bgSpline->evaluate(r) * bgError[var];
		}
		
		//F
		filter1d->filterArray(CHTHCq, fieldSize);
		filter1d->filterArray(CHTd, fieldSize);
		
		for (unsigned int R = 0; R < fieldSize; R++) {
			unsigned int Ri = R + var*fieldSize;
			stateMod[Ri] = CHTHCq[R] - CHTd[R];
		}
		
	}
	
	for (int n = 0; n < nState; n++) {
		gradient[n] = state[n] + stateMod[n];
	}

	delete[] CHTHCq;
	delete[] CHTd;
	
}

void CostFunctionR::updateHCq(double* state)
{
		
	// Clear the HCq variable
	for (int m = 0; m < mObs; m++) {
		HCq[m] = 0;
	}
	
	unsigned int fieldSize = ctrlRadii->size();
	double* field = new double[fieldSize];
	
	unsigned int bgSize = bgradii->size();
	real* Cq = new real[bgSize];
	
	for (unsigned int var = 0; var < numVars; var++) {
		for (unsigned int R = 0; R < fieldSize; R++) {
			unsigned int Ri = R + var*fieldSize;
			field[R] = state[Ri];
		}
		
		// F
		filter1d->filterArray(field, fieldSize);
		
		// D
		//ctrlSpline->solve(field);
		for (unsigned int R = 0; R < fieldSize; R++) {
			double coeff = field[R] * bgError[var];
			ctrlSpline->setCoefficient(R, coeff);
		}
		
		// P
		for (unsigned int r = 0; r < bgradii->size(); r++) {
			double R = RXform->at(r);
			Cq[r] = ctrlSpline->evaluate(R);
		}
		
		// H
		bgSpline->solve(Cq);
		for (int m = 0; m < mObs; m++) {
			double yhat =  bgSpline->evaluate(obsVector[m].getRadius());
			double weight =  obsVector[m].getWeight(var);
			HCq[m] += yhat*weight;
		}
	}
	
	delete[] field;
	delete[] Cq;
	
}

void CostFunctionR::getCq(double* Cq)
{
	// Apply series of transforms to copy of state vector	
	unsigned int fieldSize = ctrlRadii->size();
	double* field = new double[fieldSize];

	for (unsigned int var = 0; var < numVars; var++) {
		for (unsigned int R = 0; R < fieldSize; R++) {
			unsigned int Ri = R + var*fieldSize;
			field[R] = currState[Ri];
		}
		
		// F
		filter1d->filterArray(field, fieldSize);
		
		// D		
		//ctrlSpline->solve(field);
		for (unsigned int R = 0; R < fieldSize; R++) {
			double coeff = field[R] * bgError[var];
			ctrlSpline->setCoefficient(R, coeff);
		}
		
		// P
		for (unsigned int r = 0; r < bgradii->size(); r++) {
			double R = RXform->at(r);
			unsigned int i = r + var*bgradii->size();
			Cq[i] = ctrlSpline->evaluate(R);
		}
		
	}
	
	delete[] field;
}



