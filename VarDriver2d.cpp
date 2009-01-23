/*
 *  VarDriver2d.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriver2d.h"
#include "Dorade.h"
#include <iterator>
#include <fstream>
#include <cmath>
#include <QTextStream>
#include <QFile>

VarDriver2d::VarDriver2d()
	: VarDriver()
{
	numVars = 6;
}

VarDriver2d::~VarDriver2d()
{
	unsigned int maxHeights = 35; // Can I make this dynamic?
	for (unsigned int zi = 0; zi < maxHeights; zi++) {
		delete[] BG[zi];
		delete[] BGsave[zi];
	}
	delete[] BG;
	delete[] BGsave;
	delete[] scalarSpline;
	delete[] vecSpline;
	delete[] ctrlSpline;
	delete zSpline;
	delete zSplinePsi;
	delete[] RnumGridpts;
	delete[] RXform;
	delete cost2d;

}

bool VarDriver2d::run()
{
	enum VarCases { 
		analyticParaboloid,
		vortexBG,
	};
	
	int runCase = vortexBG;
	switch (runCase) {
		case analyticParaboloid:
			// Test a analytic Paraboloid
			if (!initParaboloid()) { return false; }
			runParaboloid();
			if (!finalizeParaboloid()) { return false; }
			break;
		case vortexBG:
			// Run a 2d Vortex BG field
			if (!initVortexBG()) { return false; }
			runVortexBG();
			if (!finalizeVortexBG()) { return false; }
			break;
	};
	
	cout << "Analysis complete!" << endl;
	return true;
	
}

bool VarDriver2d::initVortexBG()
{
	// Run a 2d vortex background field
	cout << "Initializing Vortex Background" << endl;

	unsigned int maxHeights = 35; // Can I make this dynamic?
	BG = new vector<real>*[maxHeights];
	BGsave = new vector<real>*[maxHeights];
	for (unsigned int zi = 0; zi < maxHeights; zi++) {
		BG[zi] = new vector<real>[numVars];
		BGsave[zi] = new vector<real>[numVars];
	}
	
	// Read in the background state
	// Read the r and v pairs from a file
	double height, radius, v, psi, h, q, rho;
	int zi = 0;
	vector<real> vIn, vBG, psiBG, hBG, qBG, rpBG;
	vector<real>* vInit = new vector<real>[maxHeights];
	ifstream vdata("./v2data.txt");
	vdata.width(14);
	while (vdata >> height >> radius >> v >> psi >> h >> q >> rho)
	{
		if (z.empty()) z.push_back (height);
		if (height != z.back()) {
			// Assign the initial background fields
			BG[zi][0] =  vBG;
			BG[zi][1] =  psiBG;
			BG[zi][2] =  hBG;
			BG[zi][3] =  qBG;
			BG[zi][4] =  rpBG;
			BG[zi][5] =  vIn;
			r.clear(); vIn.clear();
			vBG.clear(); psiBG.clear();
			hBG.clear(); qBG.clear(); rpBG.clear();
			z.push_back (height);
			zi++;
		}		
		r.push_back (radius);
		// No negative v for potential radius transform!!!
		if (v < 0)  v = 0;
		vIn.push_back (v);
		vBG.push_back (v*rho);
		psiBG.push_back (psi/1e6);
		hBG.push_back (h/1e3);
		// cout << "H: " << h/1e3 << endl;
		qBG.push_back (q);
		rpBG.push_back (rho-1.1646*exp(-1.068e-4*height));
		// cout << x.back() << " " << y.back() << endl;
	}
	
	// Assign the final strip
	BG[zi][0] =  vBG;
	BG[zi][1] =  psiBG;
	BG[zi][2] =  hBG;
	BG[zi][3] =  qBG;
	BG[zi][4] =  rpBG;
	BG[zi][5] =  vIn;
	
	// Check that z.size is not bigger than allocated array
	if (z.size() > maxHeights) {
		cerr << "Memory overflow in z direction" << z.size() << "\t" << numHeights << endl;
		return false;
	}
	
	// Set up a vertical spline on the original grid
	double wl = 2; 
	bc = SplineBase::BC_ZERO_SECOND;
	int num_nodes = z.size();
	vector<real> zBuffer;
	zBuffer.assign(z.size(), 0.0);	
	zSpline = new SplineD(&z.front(), z.size(), &zBuffer.front(), wl, SplineBase::BC_ZERO_SECOND, num_nodes);
	zSplinePsi = new SplineD(&z.front(), z.size(), &zBuffer.front(), wl, SplineBase::BC_LZERO_RSECOND, num_nodes);
	cerr << "Vertical Spline Cutoff frequency " << wl
	<< ", number of nodes " << num_nodes
	<< ", and boundary condition type " << bc << "\n";
	if (!zSpline->ok())
		cerr << "Vertical spline has a problem :(" << endl;
	numHeights = z.size();

	/* doubleign the vertical grid on the Gaussian points
	double zincr = 250;
	double zMin = 0.;
	double zMax = z.back();
	z.clear();
	num_nodes = (int)(zMax - zMin)/zincr + 1;
	int zinit = 1;
	for (unsigned int i = 0; i< numVars; i++) {
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			for (unsigned int zi = 0; zi < numHeights; zi++) {
				zBuffer[zi] = BG[zi][i].at(ri);
			}
			zSpline->solve(&zBuffer.front());
			unsigned int zi = 0;
			for (double height = zMin; height < zMax; height += zincr) {
				for (int mu = -1; mu <= 1; mu += 2) {
					double gheight = height + (0.5*sqrt(1./3.)* mu + 0.5) * zincr;
					// cout << gheight << "\t" << zSpline->evaluate(gheight) << endl;
					if (zi >= maxHeights) {
						cerr << "Memory overflow in z direction" << z.size() << "\t" << numHeights << endl;
						return false;
					}
					if (zi < numHeights) {
						BG[zi][i].at(ri) = zSpline->evaluate(gheight);
					} else {
						BG[zi][i].push_back(zSpline->evaluate(gheight));
					}
					if (zinit) z.push_back(gheight);
					zi++;
					
				}
			}
			zinit = 0;
		}
	}
	
	zSpline->setDomainGQ(&z.front(), z.size(), wl, bc, num_nodes); */

	numHeights = z.size();
	// Set up the background vectors in r space
	vecSpline = new SplineD[numHeights];
	scalarSpline = new SplineD[numHeights];
	num_nodes = r.size();
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		vecSpline[zi].setDomain(&r.front(), r.size(), wl, SplineBase::BC_LZERO_RSECOND, num_nodes);
		scalarSpline[zi].setDomain(&r.front(), r.size(), wl, SplineBase::BC_ZERO_SECOND, num_nodes);
	}

	// Iterate through the new radius grid for each variable
	double rincr = 1;
	double rMin = 0.;
	double rMax = (double)r.back();
	num_nodes = (int)((rMax - rMin)/rincr) + 1;
	r.clear();
	// And Set up the transform array for physical to potential radius
	vector<real> parmWinds;	
	RXform = new vector<real>[numHeights];
	for (unsigned int i = 0; i< numVars; i++) {
		SplineD* bgSpline;
		if (i > 1) {
			bgSpline = scalarSpline;
		} else {
			bgSpline = vecSpline;
		}		
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			bgSpline[zi].solve(&BG[zi][i].front());
			// Get the parametric winds from the input winds if data does not extend through domain
			//ParametricVortex parmVortex (&BG[zi][5].front(), &r[0],  r.size());
			/*if (i == 0) {
				scalarSpline->solve(&BG[zi][4].front());
			}*/
			for (double rad = rMin; rad< rMax; rad+= rincr)
			{
				// Gaussian quadrature points
				for (int mu = -1; mu <= 1; mu += 2) {
					double grad = rad + 0.5*sqrt(1./3.)* mu + 0.5;
					BGsave[zi][i].push_back(bgSpline[zi].evaluate(grad));
					if (i == 5) { // Raw tangential wind
						double vParm = BGsave[zi][i].back();
						/* double rp = scalarSpline->evaluate (rad);
						double rhoBar = 1.1646*exp(-1.068e-4*z[zi]);
						double newv = vParm/(rhoBar + rp); */
						//double vParm = parmVortex.getWindAtRadiusKM(grad);
						double potRad = (grad*vParm*1e3 + CoriolisF*grad*grad*1e6/2)*2/CoriolisF;
						if (potRad < 0) {
							// Inertially unstable!
							cerr << "Inertially unstable background field! Removing instability" << endl;
							potRad = 0;
						} else {
							potRad = sqrt(potRad)/1e3;
						}
						if (zi == 0) r.push_back(grad);
						RXform[zi].push_back(potRad);
						//cout << zi << ": " << grad << "\t\t" << potRad << endl;
					}
				}
			}
			BG[zi][i].clear();
			BG[zi][i] = BGsave[zi][i];
		}
	}

	// Set up the background vector in Gaussian r space
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		vecSpline[zi].setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_LZERO_RSECOND, num_nodes);
		scalarSpline[zi].setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_ZERO_SECOND, num_nodes);
	}
	cerr << "Background Spline Cutoff frequency " << wl
	<< ", number of nodes " << num_nodes
	<< ", and boundary condition type " << bc << "\n";
	
	/* Test the inverse spline
	SplineD::Debug(1);
	bgSpline[0].solveGQ(&BG[0][0].front());
	double* bgCoeff = new double[bgSpline[0].nNodes()];
	for (int i = 0; i < bgSpline[0].nNodes(); i++) {
		bgCoeff[i] = bgSpline[0].getCoefficient(i);
		//cout << "\t" << bgCoeff[i];
	}
	//cout << endl;
	const double* bgInverse = bgSpline[0].solveInverseGQ(bgCoeff);
	for (unsigned int ri = 0; ri < r.size(); ri++) {
		cout << bgInverse[ri] << "\t" ;
	}
	cout << endl; */
	
	// Set up Control space splines
	ctrlSpline = new SplineD[numHeights];
	RnumGridpts = new unsigned int[numHeights];
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		double minR = 0;
		double maxR = (int)RXform[zi].back() + 1;
		double Rwidth = maxR - minR + 1;
		RnumGridpts[zi] = (int)Rwidth;
		cout << "Z: " << zi << "\t Min / Max R: " << minR << " / " << maxR << endl;
		cout << "# control gridpoints: " << RnumGridpts[zi] << "\t" << numHeights << endl;
		R.clear();		
		for (double potRad = minR; potRad <= maxR; potRad += Rwidth/RnumGridpts[zi]) {
			R.push_back(potRad);
		}
		ctrlSpline[zi].setDomain(&R.front(), R.size(), wl, bc, RnumGridpts[zi]);
		vector<real> initCtrl(R.size(), 0.);
		ctrlSpline[zi].solve(&initCtrl.front());
	}
	
	// Read in the TC centers, create a time-based spline
	readTCcenters();
	
	// Read in the observations, process them into weights and positions
	// Don't forget R^-1 weighting
	processMetObs();
	
	cout << "Number of Observations: " << obVector.size() << endl;
	int stateSize = 0;
	for (unsigned int zi = 0; zi < z.size(); zi++)
		stateSize += RnumGridpts[zi] * numVars;
	
	cost2d = new CostFunctionRZ(obVector.size(), stateSize);
	cost2d->initialize(vecSpline, scalarSpline, ctrlSpline, zSpline, zSplinePsi,
					   &r, &z, BG,  RnumGridpts, RXform, &obVector.front()); 
	delete[] vInit;
	
	return true;
}

void VarDriver2d::runVortexBG()
{
	double CQTOL = 0.5;
	double CQRMS = 999;
	int iter=0;
	while ((CQRMS > CQTOL) and (iter < 1)) {
		iter++;
		cout << "Outer Loop Iteration: " << iter << endl;
		cost2d->minimize();
		CQRMS = updateXforms();
		cost2d->initState();
	}	
	cout << "Increment RMS Tolerance of " << CQTOL << " reached in "
		<< iter << " iterations. Writing analysis results..." << endl;
}

bool VarDriver2d::finalizeVortexBG()
{
	double wl = 2; 
	int num_nodes = vecSpline[0].nNodes();
	vector<real> empty(r.size(), 0.); 
	SplineD* vecAnalysisSpline = new SplineD(&r.front(), r.size(), &empty.front(), wl, SplineBase::BC_LZERO_RSECOND, num_nodes);
	vecAnalysisSpline->setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_LZERO_RSECOND, num_nodes);
	SplineD* scalarAnalysisSpline = new SplineD(&r.front(), r.size(), &empty.front(), wl, SplineBase::BC_ZERO_SECOND, num_nodes);
	scalarAnalysisSpline->setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_ZERO_SECOND, num_nodes);
	if (!scalarAnalysisSpline->ok())
	{
		cerr << "Analysis has a problem :(" << endl;
		return false;
	}

	vector<real>* temp = new vector<real>[numHeights];
	vector<real>* press = new vector<real>[numHeights];
	vector<real>* hstarfwd = new vector<real>[numHeights];
	vector<real>* rhoz = new vector<real>[numHeights];
	vector<real>* uz = new vector<real>[numHeights];
	vector<real>* vz = new vector<real>[numHeights];
	vector<real>* wz = new vector<real>[numHeights];
	
	for (unsigned int T = 235; T < 315; T++) {
		temp[0].push_back(T);
		hstarfwd[0].push_back(0.0);
	}
	SplineD* tSpline = new SplineD(&temp[0].front(), temp[0].size(), &hstarfwd[0].front(), wl, bc, temp[0].size());
	temp[0].clear();
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		// Get the temperature, pressure, and E
		double alt = z[zi];
		scalarSpline[zi].solveGQ(&BG[zi][2].front());
		scalarAnalysisSpline->solveGQ(&BG[zi][4].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			double hstar = scalarSpline[zi].evaluate(rad);
			double rhoPrime = scalarAnalysisSpline->evaluate(rad);
			double rho = 1.1646*exp(-1.068e-4*alt) + rhoPrime;
			hstarfwd[0].clear();
			for (unsigned int T = 235; T < 315; T++) {
				double es = (6.1078 * exp (17.2694 * (1.0 - 237.3 / (T - 273.15 + 237.3))));
				double qs = (1000.0 * 0.622 * es / (2.87*T*rho - es));
				double hs = 1005.7*T +  2.5e3*qs + 9.81*alt;
				hstarfwd[0].push_back(hstar*1000 - hs);
			}
		
			tSpline->solve(&hstarfwd[0].front());
			
			// Newton's method initial guess is ballpark
			//double T =  ((9.81 * 1.1646*exp(-1.068e-4*alt) / 1.068e-4) - 5460.)/(rho*287.);
			double T = 280;
			double hmin = 1e34;
			int iter = 0;
			while ((fabs(hmin) > 0.01) and (iter < 1000)) {
				double h = tSpline->evaluate(T);
				double hprime = tSpline->slope(T);
				if (hprime != 0) {
					T = T - h/hprime;
					hmin = h;
				}
				iter++;
			}
			
			// Now push back temp, press, and rho
			temp[zi].push_back(T);
			press[zi].push_back(T*rho*287./100);
			rhoz[zi].push_back(rho);
		}
	}
	
	// Winds
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		vecSpline[zi].solveGQ(&BG[zi][0].front());
		scalarAnalysisSpline->solveGQ(&rhoz[zi].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			double vwind = vecSpline[zi].evaluate(rad) / scalarAnalysisSpline->evaluate(rad);
			vz[zi].push_back(vwind);
		}
		vecSpline[zi].solveGQ(&BG[zi][1].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			if (rad != 0) {
				double wwind = 1e6 * (vecSpline[zi].slope(rad) / 1000) / (scalarAnalysisSpline->evaluate(rad) * rad * 1000);
				wz[zi].push_back(wwind);
			} else {
				wz[zi].push_back(0.0);
			}
		}
	}
	for (unsigned int ri = 0; ri < r.size(); ri++) {
		double rad = r[ri];
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			vecSpline[zi].solveGQ(&BG[zi][1].front());
			zSplinePsi->setCoefficient(zi, vecSpline[zi].evaluate(rad));
		}
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			double alt = z[zi];
			scalarAnalysisSpline->solveGQ(&rhoz[zi].front());
			if (rad != 0) {
				double uwind = 1e6 * -zSplinePsi->slope(alt) / (scalarAnalysisSpline->evaluate(rad) * rad * 1000);
				uz[zi].push_back(uwind);
			} else {
				uz[zi].push_back(0.0);
			}
		}
	}
	
	vector<real>** final = new vector<real>*[numHeights];
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		final[zi] = new vector<real>[11];
		for (unsigned int var = 0; var < numVars-1; var++) {
			final[zi][var] = BG[zi][var];
		}
		final[zi][5] = temp[zi];
		final[zi][6] = press[zi];
		final[zi][7] = rhoz[zi];
		final[zi][8] = uz[zi];
		final[zi][9] = vz[zi];
		final[zi][10] = wz[zi];
	}
	
	ofstream fcurve("Coefficients.out");
	ostream_iterator<string> os(fcurve, "\t ");
	*os++ = "z";
	*os++ = "r";
	*os++ = "radius";
	*os++ = "rhov";
	*os++ = "psi";
	*os++ = "hstar";
	*os++ = "qv";
	*os++ = "rhoprime";
	*os++ = "temperature";
	*os++ = "pressure";
	*os++ = "rho";
	*os++ = "u";
	*os++ = "v";
	*os++ = "w";
	
	fcurve << endl;
	ostream_iterator<double> of(fcurve, "\t ");
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		for (int ri = 0; ri < scalarAnalysisSpline->nNodes(); ri++) {
			double rad = r[ri];
			*of++ = zi;
			*of++ = ri;
			*of++ = rad;
			for (unsigned int var = 0; var < 11; var++) {
				if ((var > 1) and (var < 8)) {
					scalarAnalysisSpline->solveGQ(&final[zi][var].front());
					*of++ = scalarAnalysisSpline->getCoefficient(ri);
				} else {
					vecAnalysisSpline->solveGQ(&final[zi][var].front());
					*of++ = vecAnalysisSpline->getCoefficient(ri);
				}
			}
			fcurve << endl;	
		}
	}
	
	// Write to an asi file for plotting
	QString asifile("tcvar_analysis");
	writeAsi(asifile, final, scalarAnalysisSpline, vecAnalysisSpline);
	
	// Get the increments
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		// Get the temperature, pressure, and E
		double alt = z[zi];
		scalarSpline[zi].solveGQ(&BGsave[zi][2].front());
		scalarAnalysisSpline->solveGQ(&BGsave[zi][4].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			double hstar = scalarSpline[zi].evaluate(rad);
			double rhoPrime = scalarAnalysisSpline->evaluate(rad);
			double rho = 1.1646*exp(-1.068e-4*alt) + rhoPrime;
			hstarfwd[0].clear();
			for (unsigned int T = 235; T < 315; T++) {
				double es = (6.1078 * exp (17.2694 * (1.0 - 237.3 / (T - 273.15 + 237.3))));
				double qs = (1000.0 * 0.622 * es / (2.87*T*rho - es));
				double hs = 1005.7*T +  2.5e3*qs + 9.81*alt;
				hstarfwd[0].push_back(hstar*1000 - hs);
			}
			
			tSpline->solve(&hstarfwd[0].front());
			
			// Newton's method initial guess is ballpark
			//double T =  ((9.81 * 1.1646*exp(-1.068e-4*alt) / 1.068e-4) - 5460.)/(rho*287.);
			double T = 280;
			double hmin = 1e34;
			int iter = 0;
			while ((fabs(hmin) > 0.01) and (iter < 1000)) {
				double h = tSpline->evaluate(T);
				double hprime = tSpline->slope(T);
				if (hprime != 0) {
					T = T - h/hprime;
					hmin = h;
				}
				iter++;
			}
			
			// Now push back temp, press, and rho
			temp[zi].at(ri) -= T;
			press[zi].at(ri) -= (T*rho*287./100);
			rhoz[zi].at(ri) -= rho;
			for (unsigned int var = 0; var < numVars-1; var++) {
				BG[zi][var].at(ri) -= BGsave[zi][var].at(ri);
			}
		}
	}
	
	// Winds
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		vecSpline[zi].solveGQ(&BGsave[zi][0].front());
		vecAnalysisSpline->solveGQ(&BGsave[zi][4].front());
		double alt = z[zi];
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			double rho = scalarAnalysisSpline->evaluate(rad) + 1.1646*exp(-1.068e-4*alt);
			double vwind = vecSpline[zi].evaluate(rad) / (rho);
			vz[zi].at(ri) -= vwind;
		}
		vecSpline[zi].solveGQ(&BGsave[zi][1].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			if (rad != 0) {
				double rho = scalarAnalysisSpline->evaluate(rad) + 1.1646*exp(-1.068e-4*alt);													
				double wwind = 1e6 * (vecSpline[zi].slope(rad) / 1000) / (rho * rad * 1000);
				wz[zi].at(ri) -= wwind;
			} 
		}
	}
	for (unsigned int ri = 0; ri < r.size(); ri++) {
		double rad = r[ri];
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			vecSpline[zi].solveGQ(&BGsave[zi][1].front());
			zSplinePsi->setCoefficient(zi, vecSpline[zi].evaluate(rad));
		}
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			double alt = z[zi];
			vecAnalysisSpline->solveGQ(&BGsave[zi][4].front());
			if (rad != 0) {
				double rho = scalarAnalysisSpline->evaluate(rad) + 1.1646*exp(-1.068e-4*alt);
				double uwind = 1e6 * -zSplinePsi->slope(alt) / (rho * rad * 1000);
				uz[zi].at(ri) -= uwind;
			} 
		}
	}
	vector<real>** increment = new vector<real>*[numHeights];
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		increment[zi] = new vector<real>[11];
		for (unsigned int var = 0; var < numVars-1; var++) {
			increment[zi][var] = BG[zi][var];
		}
		increment[zi][5] = temp[zi];
		increment[zi][6] = press[zi];
		increment[zi][7] = rhoz[zi];
		increment[zi][8] = uz[zi];
		increment[zi][9] = vz[zi];
		increment[zi][10] = wz[zi];
	}
	
	
	// Write to an asi file for plotting
	asifile = "tcvar_increment";
	writeAsi(asifile, increment, scalarAnalysisSpline, vecAnalysisSpline);
	
	// Write some info to a file for plotting
	ofstream acurve("Analysis.out");
	ostream_iterator<string> as(acurve, "\t ");
	*as++ = "Height";
	*as++ = "Radius";
	*as++ = "rhov Background";
	*as++ = "rhov Increment";
	*as++ = "rhov Analysis";
	*as++ = "psi Background";
	*as++ = "psi Increment";
	*as++ = "psi Analysis";
	*as++ = "hstar Background";
	*as++ = "hstar Increment";
	*as++ = "hstar Analysis";
	*as++ = "qv Background";
	*as++ = "qv Increment";
	*as++ = "qv Analysis";
	*as++ = "rho prime Background";
	*as++ = "rho prime Increment";
	*as++ = "rho prime Analysis";
	*as++ = "temperature Background";
	*as++ = "temperature Increment";
	*as++ = "temperature Analysis";
	*as++ = "pressure Background";
	*as++ = "pressure Increment";
	*as++ = "pressure Analysis";
	*as++ = "rho Background";
	*as++ = "rho Increment";
	*as++ = "rho Analysis";
	*as++ = "u Background";
	*as++ = "u Increment";
	*as++ = "u Analysis";	
	*as++ = "v Background";
	*as++ = "v Increment";
	*as++ = "v Analysis";
	*as++ = "w Background";
	*as++ = "w Increment";	
	*as++ = "w Analysis";
	
	acurve << endl;
	ostream_iterator<double> af(acurve, "\t ");
	double xs = (vecSpline[0].Xmax() - vecSpline[0].Xmin()) / 100.0;
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		for (double x = vecSpline[zi].Xmin(); x <= vecSpline[zi].Xmax(); x += xs)
		{
			*af++ = z[zi];
			*af++ = x;
			for (unsigned int var = 0; var < numVars-1; var++) {
				if (var > 1) {
					scalarSpline[zi].solveGQ(&BGsave[zi][var].front());
					scalarAnalysisSpline->solveGQ(&increment[zi][var].front());
					double bg = scalarSpline[zi].evaluate (x);
					double incr = scalarAnalysisSpline->evaluate (x);
					*af++ = bg;
					*af++ = incr;
					*af++ = bg+incr;
				} else {
					vecSpline[zi].solveGQ(&BGsave[zi][var].front());
					vecAnalysisSpline->solveGQ(&increment[zi][var].front());
					double bg = vecSpline[zi].evaluate (x);
					double incr = vecAnalysisSpline->evaluate (x);
					*af++ = bg;
					*af++ = incr;
					*af++ = bg+incr;
				}
			}
			for (unsigned int var = numVars-1; var < 11; var++) {
				if (var < 8) {
					scalarSpline[zi].solveGQ(&final[zi][var].front());
					scalarAnalysisSpline->solveGQ(&increment[zi][var].front());
					double fnl = scalarSpline[zi].evaluate (x);
					double incr = scalarAnalysisSpline->evaluate (x);
					*af++ = fnl-incr;
					*af++ = incr;
					*af++ = fnl;
				} else {
					vecSpline[zi].solveGQ(&final[zi][var].front());
					vecAnalysisSpline->solveGQ(&increment[zi][var].front());
					double fnl = vecSpline[zi].evaluate (x);
					double incr = vecAnalysisSpline->evaluate (x);
					*af++ = fnl-incr;
					*af++ = incr;
					*af++ = fnl;
				}					
			}			
			acurve << endl;
		}
	}
	
	delete scalarAnalysisSpline;
	delete vecAnalysisSpline;
	delete tSpline;
	delete[] temp;
	delete[] press;
	delete[] hstarfwd;
	delete[] rhoz;
	delete[] uz;
	delete[] vz;
	delete[] wz;
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		delete[] final[zi];
		delete[] increment[zi];
	}
	delete[] final;
	delete[] increment;
	cost2d->finalize();
	
	return true;
}

bool VarDriver2d::initParaboloid()
{
	
	costAnalytic = CostFunctionAnalytic();
	costAnalytic.initialize();
	return true;
}

void VarDriver2d::runParaboloid()
{
	costAnalytic.minimize();
}

bool VarDriver2d::finalizeParaboloid()
{
	
	costAnalytic.finalize();
	return true;
}

double VarDriver2d::updateXforms()
{
	int wl = 2; 
	int num_nodes = vecSpline[0].nNodes();
	double *Cq = new double[r.size()*z.size()*numVars-1];
	cost2d->getCq(Cq);
	double RMS = 0;
	vector<real> incr(r.size(), 0.); 
	// Need to update rho with perturbation
	SplineD* vecIncrSpline = new SplineD(&r.front(), r.size(), &incr.front(), wl, SplineBase::BC_LZERO_RSECOND, num_nodes);
	vecIncrSpline->setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_LZERO_RSECOND, num_nodes);
	SplineD* scalarIncrSpline = new SplineD(&r.front(), r.size(), &incr.front(), wl, SplineBase::BC_ZERO_SECOND, num_nodes);
	scalarIncrSpline->setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_ZERO_SECOND, num_nodes);

	if (!vecIncrSpline->ok())
	{
		cerr << "Analysis has a problem :(" << endl;
		return false;
	} 
	
	// Increment the variables
	double RxRMS = 0;
	for (unsigned int zi = 0; zi < z.size(); zi++) {	
		double rhoBar = 1.1646*exp(-1.068e-4*z[zi]);
		for (unsigned int var = 0; var < numVars-1; var++) {
			if (var > 1) {
				scalarSpline[zi].solveGQ(&BG[zi][var].front());
				BG[zi][var].clear();
				
				incr.clear();
				for (unsigned int p = 0; p < r.size(); p++) {
					unsigned int i = p + var*r.size() + zi*(numVars-1)*r.size();		
					incr.push_back(Cq[i]);
				}
				scalarIncrSpline->solveGQ(&incr.front());
				for (unsigned int ri = 0; ri < r.size(); ri++) {
					double rad = r[ri];
					double bg = scalarSpline[zi].evaluate (rad);
					double cq = scalarIncrSpline->evaluate (rad);
					//cout << "Incr:" << ri << "\t" << bg << "\t" << cq << endl;
					double newbg = bg + cq;
					BG[zi][var].push_back(newbg);
				}
			} else {
				vecSpline[zi].solveGQ(&BG[zi][var].front());
				BG[zi][var].clear();
				
				incr.clear();
				for (unsigned int p = 0; p < r.size(); p++) {
					unsigned int i = p + var*r.size() + zi*(numVars-1)*r.size();
					incr.push_back(Cq[i]);
				}
				vecIncrSpline->solveGQ(&incr.front());
				for (unsigned int ri = 0; ri < r.size(); ri++) {
					double rad = r[ri];
					double bg = vecSpline[zi].evaluate (rad);
					double cq = vecIncrSpline->evaluate (rad);
					//cout << "Incr:" << ri << "\t" << bg << "\t" << cq << endl;
					double newbg = bg + cq;
					BG[zi][var].push_back(newbg);
				}
			}
		}
		
		// Compute an updated transform
		vecSpline[zi].solveGQ(&BG[zi][0].front());
		scalarIncrSpline->solveGQ(&BG[zi][4].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			double rhov = vecSpline[zi].evaluate (rad);
			double rp = scalarIncrSpline->evaluate (rad);
			double newv = rhov/(rhoBar + rp);
			double potRad = (rad*newv*1e3 + CoriolisF*rad*rad*1e6/2)*2/CoriolisF;
			if (potRad < 0) {
				// Inertially unstable!
				cerr << "Inertially unstable winds! Removing instability" << endl;
				potRad = 0;
			} else {
				potRad = sqrt(potRad)/1e3;
			}
			RxRMS += ((RXform[zi][ri] - potRad)*(RXform[zi][ri] - potRad));
			//cout << RxRMS << endl;
			RXform[zi][ri] = potRad;
		}
	}
	RxRMS = sqrt(RxRMS/(r.size()*z.size()));
	delete vecIncrSpline;
	delete scalarIncrSpline;
	
	// Compute the RMS of the increment vector
	double vMax = 0;
	double bgRMS = 0;
	for (unsigned int var = 0; var < numVars-1; var++) {
		RMS = 0;
		bgRMS = 0;
		double maxIncr = 0;
		double maxPercent = 0;
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			for (unsigned int p = 0; p < r.size(); p++) {
				unsigned int i = p + var*r.size() + zi*(numVars-1)*r.size();
				if (Cq[i] > maxIncr) maxIncr = Cq[i];
				if (BG[zi][var].at(p) != 0) {
					if (100*abs(Cq[i]/BG[zi][var].at(p)) > maxPercent) maxPercent = 100*abs(Cq[i]/BG[zi][var].at(p));
				}
				RMS += Cq[i]*Cq[i];
				//cout << i << "\t" << p << "\t" << var << "\t" << Cq[i] << endl;
				bgRMS += BG[zi][var].at(p)*BG[zi][var].at(p);
			}
		}
		RMS = sqrt(RMS/(r.size()*z.size()));
		bgRMS = sqrt(bgRMS/(r.size()*z.size()));
		cout << "Variable " << var << " RMS = " << bgRMS << endl;
		cout << "Max Increment = " << maxIncr << endl;
		cout << "Max % Increment = " << maxPercent << " %" << endl;
		cout << "RMS Increment = " << RMS << "\t( " << 100*RMS/bgRMS << " %)" << endl << endl;
		if (var == 0) vMax = 100*RMS/bgRMS;
	}
	
	RMS = 0;
	for (unsigned int ri = 0; ri < (numVars-1)*r.size()*z.size(); ri++) {
		RMS += Cq[ri]*Cq[ri];
		//cout << "Cq:\t" << r[ri] << "\t" << Cq[ri] << endl;
	}
	RMS = sqrt(RMS/(numVars*r.size()*z.size()));
	cout << "Cq (delta) RMS = " << RMS << endl;
	delete[] Cq;
	
	// Use Rx criteria
	cout << "R transform RMS difference: " << RxRMS << endl << endl;
	return RxRMS;
	// Use the maximum Vt as the criterion
	//return vMax;
	// Use total RMS (not very good measure)
	//return RMS;
	
}	

void VarDriver2d::EvalSpline (SplineD* spline, ostream &out)
{
	ostream_iterator<double> of(out, "\t ");
	
	double x = spline->Xmin();
	double xs = (spline->Xmax() - x) / 2000.0;
	for (; x <= spline->Xmax(); x += xs)
	{
		*of++ = x;
		*of++ = spline->evaluate (x);
		*of++ = spline->slope (x);
		out << endl;
	}
}


void VarDriver2d::processMetObs()
{

	vector<real> rhoP;
	
	// Check the data directory for files
	QDir dataPath("./vardata");
	dataPath.setFilter(QDir::Files);
	dataPath.setSorting(QDir::Name);
	QStringList filenames = dataPath.entryList();

	int processedFiles = 0;
	QList<MetObs>* metData = new QList<MetObs>;
	cout << "Found " << filenames.size() << " data files to read..." << endl;
	for (int i = 0; i < filenames.size(); ++i) {
		metData->clear();
		QString file = filenames.at(i);
		QStringList fileparts = file.split(".");
		if (fileparts.isEmpty()) {
			cout << "Unknown file! " << file.toAscii().data() << endl;
			continue;
		}
		QString suffix = fileparts.last();
		QString prefix = fileparts.first();
		if (prefix == "swp") {
			// Switch it to suffix
			suffix = "swp";
		}
		cout << "Processing " << file.toAscii().data() << " of type " << suffix.toAscii().data() << endl;
		QFile metFile(dataPath.filePath(file));

		// Read different types of files
		switch (dataSuffix.value(suffix)) {
			case (frd):
				if (!read_frd(metFile, metData))
					cout << "Error reading frd file" << endl;
				break;
			case (cls):
				if (!read_cls(metFile, metData))
					cout << "Error reading frd file" << endl;
				break;
			case (sec):
				if (!read_sec(metFile, metData))
					cout << "Error reading min file" << endl;
				break;
			case (ten):
				if (!read_ten(metFile, metData))
					cout << "Error reading ten file" << endl;
				break;
			case (swp):
				if (!read_dorade(metFile, metData))
					cout << "Error reading swp file" << endl;
				break;
			default:
				cout << "Unknown data type, skipping..." << endl;
				continue;
		}
		
		processedFiles++;
		
		// Process the metObs into Observations
		QDateTime startTime = tcVector.front().getTime();
		QDateTime endTime = tcVector.back().getTime();
		for (int i = 0; i < metData->size(); ++i) {
			MetObs metOb = metData->at(i);
			QDateTime obTime = metOb.getTime();
			/* cout << startTime.toString(Qt::ISODate).toAscii().data()
			<< "\t" << obTime.toString(Qt::ISODate).toAscii().data()
			<< "\t" << endTime.toString(Qt::ISODate).toAscii().data() << endl; */
			if ((obTime < startTime) or (obTime > endTime)) continue;
			int tci = startTime.secsTo(obTime);
			if ((tci < 0) or (tci > (int)tcVector.size())) {
				cout << "Time problem with observation " << tci << endl;
				continue;
			}
			// Our generic observation
			Observation varOb;

			// Get the radius
			double latrad = tcVector[tci].getLat() * Pi/180.0;
			double fac_lat = 111.13209 - 0.56605 * cos(2.0 * latrad)
			+ 0.00012 * cos(4.0 * latrad) - 0.000002 * cos(6.0 * latrad);
			double fac_lon = 111.41513 * cos(latrad)
			- 0.09455 * cos(3.0 * latrad) + 0.00012 * cos(5.0 * latrad);
			double y = (metOb.getLat() - tcVector[tci].getLat())*fac_lat;
			double x = (metOb.getLon() - tcVector[tci].getLon())*fac_lon;
			double rad = sqrt(x*x + y*y);
			// Make sure the ob is in the domain
			if ((rad < r.front()) or (rad > r.back()) or
				(metOb.getAltitude() > z.back()))
				continue;

			varOb.setRadius(rad);
			float Um = tcVector[tci].getUmean();
			float Vm = tcVector[tci].getVmean();
			varOb.setAltitude(metOb.getAltitude());
			
			// Use the background rho prime to correct radar velocities -- this could be done more efficiently
			rhoP.clear();
			for (unsigned int zi = 0; zi < z.size(); zi++) {
				scalarSpline[zi].solveGQ(&BG[zi][4].front());
				rhoP.push_back(scalarSpline[zi].evaluate(rad));
			}
			zSpline->solve(&rhoP.front());
			double rhoBar = 1.1646*exp(-1.068e-4*metOb.getAltitude());
			double rhoPrimeBG = zSpline->evaluate(metOb.getAltitude());

			// Initialize the weights
			varOb.setWeight(0., 0);
			varOb.setWeight(0., 1);
			varOb.setWeight(0., 2);
			varOb.setWeight(0., 3);
			varOb.setWeight(0., 4);
			varOb.setWeight(0., 5);
			double u, v, w, rho, qv, hstar, rhov, rhou, rhow; 
			switch (metOb.getObType()) {
				case (MetObs::dropsonde):
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					w = metOb.getVerticalVelocity();
					rho = metOb.getDryDensity();
					qv = metOb.getQv();
					hstar = metOb.getMoistSaturationStaticEnergy()/1e3;
					
					// Separate obs for each measurement
					// rho v 1 m/s error
					if ((u != -999) and (rho != -999)) {
						varOb.setWeight(1., 0);
						rhov = rho*(-(u - Um)*y + (v-Vm)*x)/rad;
						varOb.setOb(rhov);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						// rho u 1 m/s error
						varOb.setWeight(1., 1);
						rhou = rho*((u - Um)*x + (v-Vm)*y)/rad;
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);
					}
					if ((w != -999) and (rho != -999)) {
						// rho w 1.5 m/s error
						varOb.setWeight(1., 2);
						rhow = rho*w;
						varOb.setOb(rhow);
						varOb.setError(1.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 2);
					}
					if (hstar != -999) {
						// hstar 20 kJ error
						varOb.setWeight(1., 3);
						varOb.setOb(hstar);
						varOb.setError(5.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 3);
					}
					if (qv != -999) {
						// Qv 2 g/kg error
						varOb.setWeight(1., 4);
						varOb.setOb(qv);
						varOb.setError(2.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 4);
					}
					if (rho != -999) {
						// Rho prime .1 kg/m^3 error
						varOb.setWeight(1., 5);
						varOb.setOb(rho-rhoBar);
						varOb.setError(0.01);
						obVector.push_back(varOb);
						varOb.setWeight(0., 5);
					}
					
					break;
					
				case (MetObs::flightlevel):
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					w = metOb.getVerticalVelocity();
					rho = metOb.getDryDensity();
					qv = metOb.getQv();
					hstar = metOb.getMoistSaturationStaticEnergy()/1e3;
					
					// Separate obs for each measurement
					// rho v 1 m/s error
					if ((u != -999) and (rho != -999)) {
						varOb.setWeight(1., 0);
						rhov = rho*(-(u - Um)*y + (v-Vm)*x)/rad;
						varOb.setOb(rhov);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						// rho u 1 m/s error
						varOb.setWeight(1., 1);
						rhou = rho*((u - Um)*x + (v-Vm)*y)/rad;
						varOb.setOb(rhou);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);
					}
					if ((w != -999) and (rho != -999)) {
						// rho w 1 dm/s error
						varOb.setWeight(1., 2);
						rhow = rho*w;
						varOb.setOb(rhow);
						varOb.setError(0.25);
						obVector.push_back(varOb);
						varOb.setWeight(0., 2);
					}
					if (hstar != -999) {
						// hstar 20 kJ error
						varOb.setWeight(1., 3);
						varOb.setOb(hstar);
						varOb.setError(5.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 3);
					}
					if (qv != -999) {
						// Qv 2 g/kg error
						varOb.setWeight(1., 4);
						varOb.setOb(qv);
						varOb.setError(2.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 4);
					}
					if (rho != -999) {
						// Rho prime .1 kg/m^3 error
						varOb.setWeight(1., 5);
						varOb.setOb(rho-rhoBar);
						varOb.setError(0.01);
						obVector.push_back(varOb);
						varOb.setWeight(0., 5);
					}
					
					break;
					
				case (MetObs::radar):
					// rho v 1 m/s error
					double az = metOb.getAzimuth()*Pi/180.;
					double el = metOb.getElevation()*Pi/180.;
					double uWgt = (x*sin(az)*cos(el) + y*cos(az)*cos(el))/rad;
					double vWgt = (x*cos(az)*cos(el) - y*sin(az)*cos(el))/rad;
					double wWgt = sin(el);
					
					// Fall speed
					double Z = metOb.getReflectivity();
					double H = metOb.getAltitude();
					double ZZ=pow(10.0,(Z*0.1));
					double hlow= 5600 - 1000 * .5; 
					double hhi= hlow + 1000;
					/* density correction term (rhoo/rho)*0.45 [rho(Z)=rhoo exp-(z/H), where 
						C  H is the scale height = 9.58125 from Gray's inner 2 deg composite] 
					C 0.45 density correction from Beard (1985, JOAT pp 468-471) 
					 Adjusted to use Jordan hydrostatic scale height -MB */
					double DCOR=exp(0.45*metOb.getAltitude()*0.0001068);
					//C The snow relationship (Atlas et al., 1973) --- VT=0.817*Z**0.063  (m/s) 
					double VTS=-DCOR * (0.817*pow(ZZ,(double)0.063));
					/* if(irsw.gt.0) then
					C The rain relationship --- from Willis analytical-gamma distribution 
					TERM1=7.331/ZZ**0.010022 
					TERM2=0.14034*ZZ**0.095238 
					VTR=-DCOR * (5.5011E+09/(TERM1+TERM2)**10.5) 
					else 
					C The rain relationship (Joss and Waldvogel,1971) --- VT=2.6*Z**.107 (m/s) */
					double VTR=-DCOR * (2.6*pow(ZZ,(double).107));
					/* endif 
					C test if height is in the transition region between SNOW and RAIN
					C  defined as hlow in km < H < hhi in km
					C  if in the transition region do a linear weight of VTR and VTS */
					if ((Z > 20) and (Z <= 30)) {
						double WEIGHTR=(Z-20)/(10);
						double WEIGHTS=1.-WEIGHTR;
						VTS=(VTR*WEIGHTR+VTS*WEIGHTS)/(WEIGHTR+WEIGHTS);
					} else if (Z > 30) {
						VTS=VTR;
					}
					double w_term=VTR*(hhi-H)/1000 + VTS*(H-hlow)/1000;  
					if (H < hlow) w_term=VTR; 
					if (H > hhi) w_term=VTS;
								
					double Vdopp = metOb.getRadialVelocity() - w_term*sin(el) - Um*sin(az)*cos(el) - Vm*cos(az)*cos(el);

					varOb.setWeight(vWgt, 0);
					varOb.setWeight(uWgt, 1);
					varOb.setWeight(wWgt, 2);
					// double rhopWgt = -Vdopp;
					//varOb.setWeight(rhopWgt, 5);
					varOb.setError(metOb.getSpectrumWidth());
					varOb.setOb((rhoBar+rhoPrimeBG)*Vdopp);
					obVector.push_back(varOb);
						
					break;
			}
		} 
		cout << obVector.size() << " total observations." << endl;
	}
	
	delete metData;
	if (!processedFiles) {
		cout << "No files processed, nothing to do :(" << endl;
		// return 0;
	} else {
		cout << "Finished processing " << processedFiles << " files." << endl;
	}
		
}

bool VarDriver2d::writeAsi(const QString& fileName, vector<real>** fields, SplineD* scalar, SplineD* vectorSpline)
{
	QString outFileName;
	if(QDir::isAbsolutePath(fileName)) {
		outFileName = fileName;
	}
	else {
		outFileName = QDir::current().filePath(fileName);
	}
	
	// Write out the CAPPI to an asi file
	
	// Initialize header
	int id[511];
	for (int n = 1; n <= 510; n++) {
		id[n]=-999;
	}
	
	// Calculate headers
	QStringList fieldNames;
	fieldNames << "RV" << "SF" << "HS" << "QV" << "RP" << "T" << "P" << "RO" << "U" << "V" << "W";
	id[175] = fieldNames.size();
    for(int n = 0; n < id[175]; n++) {
		QString name_1 = fieldNames.at(n).left(1);
        QString name_2 = fieldNames.at(n).mid(1,1);
		int int_1 = *name_1.toAscii().data();
		int int_2 = *name_2.toAscii().data();
		id[176 + (5 * n)] = (int_1 * 256) + int_2;
		id[177 + (5 * n)] = 8224;
		id[178 + (5 * n)] = 8224;
		id[179 + (5 * n)] = 8224;
		id[180 + (5 * n)] = 1;
	}
	
	// Cartesian file
	id[16] = 17217;
	id[17] = 21076;
	
	/* Lat and Lon
	id[33] = (int)latReference;
	id[34] = (int)((latReference - (float)id[33]) * 60.);
	id[35] = (int)((((latReference - (float)id[33]) * 60.) - (float)id[34]) * 60.) * 100;
	if (lonReference < 0) {
		lonReference += 360.;
	}
	id[36] = (int)lonReference;
	id[37] = (int)((lonReference - (float)id[36]) * 60.);
	id[38] = (int)((((lonReference - (float)id[36]) * 60.) - (float)id[37]) * 60.) * 100; */

	id[33] = 0;
	id[34] = 0;
	id[35] = 0;
	id[36] = 0;
	id[37] = 0;
	id[38] = 0;
	id[40] = 90;
	
	// Scale factors
	id[68] = 100;
	id[69] = 64;
	
	// X Header
	double rincr = 1;
	id[160] = (int)(0);
	id[161] = (int)(scalar->nNodes() * 100);
	id[162] = (int)scalar->nNodes();
	id[163] = (int)rincr * 1000;
	id[164] = 1;
	
	// Y Header
	double zincr = 0.25;
	id[165] = (int)(z.front());
	id[166] = (int)(z.back()/10);
	id[167] = (int)z.size();
	id[168] = (int)(zincr * 1000);
	id[169] = 2;
	
	// Z Header
	id[170] = 1000;
	id[171] = 1000;
	id[172] = 1;
	id[173] = 1000;
	id[174] = 3;
	
	// Number of radars
	id[303] = 1;
	
	// Index of center
	id[309] = (int)(1);
	id[310] = (int)((1 - z.front()) * 100);
	id[311] = 0;
	
	// Write ascii file for grid2ps
	//Message::toScreen("Trying to write cappi to "+outFileName);
	outFileName += ".asi";
	QFile asiFile(outFileName);
	if(!asiFile.open(QIODevice::WriteOnly)) {
		cout << "Can't open CAPPI file for writing" << endl;
		return false;
	}
	
	QTextStream out(&asiFile);
	
	// Write header
    int line = 0;
	for (int n = 1; n <= 510; n++) {
		line++;
		out << qSetFieldWidth(8) << id[n];
		if (line == 10) {
			out << endl;
            line = 0;
		}
	}
	
	// Write data
	for(int k = 0; k < 1; k++) {
		out << reset << "level" << qSetFieldWidth(2) << k+1 << endl;
		for(int j = 0; j < int(z.size()); j++) {
			out << reset << "azimuth" << qSetFieldWidth(3) << j+1 << endl;
			for(int n = 0; n < fieldNames.size(); n++) {
				if (n > 1) {
					scalar->solveGQ(&fields[j][n].front());
				} else {
					vectorSpline->solveGQ(&fields[j][n].front());
				}
				out << reset << left << fieldNames.at(n) << endl;
				int line = 0;
				const real* curve;
				if (n > 1) {
					curve = scalar->curve();
				} else {
					curve = vectorSpline->curve();
				}
				for (int i = 0; i < int(scalar->nNodes());  i++){
					if (n > 1) {
						out << reset << qSetRealNumberPrecision(3) << scientific << qSetFieldWidth(10) << curve[i];
					} else {
						out << reset << qSetRealNumberPrecision(3) << scientific << qSetFieldWidth(10) << curve[i];
					}						
					line++;
					if (line == 8) {
						out << endl;
						line = 0;
					}
				}
				if (line != 0) {
					out << endl;
				}
			}
		}
	}
	
	return true;
	
}     

