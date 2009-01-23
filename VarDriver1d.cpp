/*
 *  VarDriver1d.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriver1d.h"
#include "CostFunctionAnalytic.h"
#include "Dorade.h"
#include <iterator>
#include <fstream>
#include <cmath>
#include <QTextStream>

VarDriver1d::VarDriver1d()
	: VarDriver()
{
	numVars = 6;
	BG = new vector<real>[numVars];
	BGsave = new vector<real>[numVars];
		
}

VarDriver1d::~VarDriver1d()
{
	delete[] BG;
	delete[] BGsave;
}

bool VarDriver1d::run()
{
	enum VarCases { 
		analyticParaboloid,
		parametricBG,
	};
	
	int runCase = parametricBG;
	switch (runCase) {
		case analyticParaboloid:
			// Test a analytic Paraboloid
			if (!initParaboloid()) { return false; }
			runParaboloid();
			if (!finalizeParaboloid()) { return false; }
			break;
		case parametricBG:
			// Test a 1d Parametric Vortex BG field with synthetic obs
			if (!initParametricBG()) { return false; }
			runParametricBG();
			if (!finalizeParametricBG()) { return false; }
			break;
	};
	
	cout << "Analysis complete!" << endl;
	return true;
	
}

bool VarDriver1d::initParametricBG()
{
	// Run a 1d test with a parametric background field and synthetic observations
	cout << "Initializing Parametric Background Test" << endl;
	
	
	// Read in the background state
	// Read the r and v pairs from a file
	double f, v, u, w, h, q, rho;
	double base = 0;
	int i = 0;
	vector<real> vIn, vBG, uBG, wBG, hBG, qBG, rpBG;
	ifstream vdata("/Users/mbell/Development/tcvar/vdata.txt");
	vdata.width(14);
	while (vdata >> f >> v >> u >> w >> h >> q >> rho)
	{
		if (++i == 1)
		{
			base = f;
		}
		f -= base;
		
		r.push_back (f);
		vIn.push_back (v);
		vBG.push_back (v*rho);
		uBG.push_back (u*rho);
		wBG.push_back (w*rho);
		hBG.push_back (h/1e3);
		// cout << "H: " << h/1e3 << endl;
		qBG.push_back (q);
		rpBG.push_back (rho-1.1646*exp(-1.068e-4*2000));
		// cout << x.back() << " " << y.back() << endl;
	}

	// Assign the initial background fields
	BG[0] =  vBG;
	BG[1] =  uBG;
	BG[2] =  wBG;
	BG[3] =  hBG;
	BG[4] =  qBG;
	BG[5] =  rpBG;
	
	// Set up the initial transform spline on the given background grid
	// Cutoff wavelength set to two
	double wl = 2; 
	bc = SplineBase::BC_ZERO_SECOND;
	int num_nodes = r.size();
	cerr << "Background Spline Cutoff frequency " << wl
	<< ", number of nodes " << num_nodes
	<< ", and boundary condition type " << bc << "\n";

	// Set up the background vectors in r space
	bgSpline = new SplineD(&r.front(), r.size(), &vBG[0], wl, bc, num_nodes);
	if (!bgSpline->ok())
	{
		cerr << "bgSpline setup failed." << endl;
		return false;
	}
	
	// Iterate through the new radius grid for each variable
	double rincr = 1;
	double rMin = 0.;
	double rMax = (int)r.back()+1;
	for (unsigned int i = 0; i< numVars; i++) {
		bgSpline->solve(&BG[i].front());
		for (double rad = rMin; rad< rMax; rad+= rincr)
		{
			BGsave[i].push_back(bgSpline->evaluate(rad));
		}
		BG[i] = BGsave[i];
	}

	// Get the parametric winds from the input winds if data does not extend through domain
	// ParametricVortex parmVortex (&vIn[0], &r[0],  r.size());
	
	vector<real> parmWinds;
	r.clear();
	bgSpline->solve(&vIn[0]);
	for (double rad = rMin; rad< rMax; rad+= rincr)
	{
		//double vParm = parmVortex.getWindAtRadiusKM(rad);
		double vParm = bgSpline->evaluate(rad);
		double potRad = sqrt((rad*vParm*1e3 + CoriolisF*rad*rad*1e6/2)*2/CoriolisF)/1e3;
		r.push_back(rad);
		RXform.push_back(potRad);
		cout << "R/r: " << rad << "\t" << potRad << endl;
		parmWinds.push_back(vParm);
	}
	
	num_nodes = r.size();
	cerr << "Background Spline Cutoff frequency " << wl
	<< ", number of nodes " << num_nodes
	<< ", and boundary condition type " << bc << "\n";
	// Set up the background vector in r space
	bgSpline->setDomain(&r.front(), r.size(), wl, bc, num_nodes);
	if (!bgSpline->ok())
	{
		cerr << "bgSpline setup failed." << endl;
		return false;
	}
	
	// Set up the inverse coordinate transform
	SplineD::Debug(0);
	RXformSpline = new SplineD(&RXform.front(), r.size(), &r.front(), wl, bc, num_nodes);
	if (RXformSpline->ok())
	{
		// write the curve to a file
		ofstream fcurve("RXform.out");
		EvalSpline (RXformSpline, fcurve);
	} else {
		cerr << "RXformSpline setup failed." << endl;
		return false;
	}
	
	// Get the optimal nodes for the potential radius grid
	double minR = RXform[0];
	double maxR = RXform[r.size()-1];
	cout << "Min / Max R: " << minR << " / " << maxR << endl;
	minR = 0;
	maxR = (int)maxR + 1;
	double Rwidth = maxR - minR + 1;
	int RnumGridpts = (int)Rwidth;
	cout << "# control gridpoints: " << RnumGridpts << endl;

	// Implement Gaussian mish at some point
	//int mu = 4;

	// Cutoff wavelength of 2 for this mesh
	wl = 2;
	cerr << "Control Spline Cutoff frequency " << wl
	<< ", number of nodes " << RnumGridpts
	<< ", and boundary condition type " << bc << "\n";
	
	vector<real> initCtrl;
	initCtrl.assign(RnumGridpts, 0.0);
	R.clear();
	bgSpline->solve(&BG[0].front());
	for (double potRad = minR; potRad <= maxR; potRad += Rwidth/RnumGridpts) {
		R.push_back(potRad);
		double vbg = bgSpline->evaluate(RXformSpline->evaluate(potRad));
		double physrad = (-vbg + sqrt(vbg*vbg + CoriolisF*CoriolisF*potRad*potRad*1e6))/(CoriolisF*1e3);
		cout << "r/R: " << physrad << "\t" << potRad << endl;
		rXform.push_back(physrad);
	}
	ctrlSpline = new SplineD(&R.front(), R.size(), &initCtrl.front(), wl, bc, RnumGridpts);

	if (!ctrlSpline->ok()) { 
		cerr << "ctrlSpline setup failed." << endl;
		return false;
	}
	
	// Read in the TC centers, create a time-based spline
	readTCcenters();
	
	// Read in the observations, process them into weights and positions
	// Don't forget R^-1 weighting
	processMetObs();
	
	// Synthetic obs for now
	int obCase = -1;
	Observation synthOb = Observation();
	double ob;
	switch (obCase) {
		case 0:
			// Single ob
			synthOb.setRadius(40.);
			synthOb.setWeight(1., 0);
			synthOb.setWeight(0., 1);
			synthOb.setWeight(0., 2);
			synthOb.setWeight(0., 3);
			synthOb.setWeight(0., 4);
			synthOb.setWeight(0., 5);
			ob = 70;
			synthOb.setOb(ob);
			obVector.push_back(synthOb);

			synthOb.setWeight(0., 0);
			synthOb.setWeight(1., 1);
			ob = 10;
			synthOb.setOb(ob);
			obVector.push_back(synthOb);

			synthOb.setWeight(0., 1);
			synthOb.setWeight(1., 2);
			ob = 5;
			synthOb.setOb(ob);
			obVector.push_back(synthOb);
			
			synthOb.setWeight(0., 2);
			synthOb.setWeight(1., 3);
			ob = 10000;
			synthOb.setOb(ob);
			obVector.push_back(synthOb);
			
			synthOb.setWeight(0., 3);
			synthOb.setWeight(1., 4);
			ob = .2;
			synthOb.setOb(ob);
			obVector.push_back(synthOb);
			
			synthOb.setWeight(0., 4);
			synthOb.setWeight(1., 5);
			ob = .1;
			synthOb.setOb(ob);
			obVector.push_back(synthOb);
			
			break;
		case 1:
			// Two Obs
			synthOb.setRadius(50.);
			synthOb.setWeight(1., 0);
			synthOb.setOb(41.91);
			obVector.push_back(synthOb);
			synthOb.setRadius(125.);
			synthOb.setWeight(1., 0);
			synthOb.setOb(31.30);
			obVector.push_back(synthOb);
			break;
		case 2:
			// Sparse Obs
			for (unsigned int rt = 0; rt < r.size(); rt += 2)
			{
				synthOb.setRadius(r.at(rt));
				synthOb.setWeight(1., 0);
				synthOb.setOb(vIn.at(rt));
				obVector.push_back(synthOb);
			}
			break;
	}

	cout << "Number of Observations: " << obVector.size() << endl;
	int stateSize = RnumGridpts * numVars;
	cost1d = new CostFunctionR(obVector.size(), stateSize);
	cost1d->initialize(bgSpline, &r, BG, ctrlSpline, &R, &RXform, &rXform, &obVector.front()); 

	return true;
}

void VarDriver1d::runParametricBG()
{
	double CQTOL = 0.05;
	double CQRMS = 999;
	int iter=0;
	while ((CQRMS > CQTOL) and (iter < 20)) {
		iter++;
		cout << "Outer Loop Iteration: " << iter << endl;
		cost1d->minimize();
		CQRMS = updateXforms();
		cost1d->initState();
	}
	cout << "Increment RMS Tolerance of " << CQTOL << " reached in "
		<< iter << " iterations. Writing analysis results..." << endl;
}

bool VarDriver1d::finalizeParametricBG()
{
	double wl = 2; 
	int num_nodes = r.size();
	SplineD* analysisSpline = new SplineD(&r.front(), r.size(), &BG[0].front(), wl, bc, num_nodes);
	if (analysisSpline->ok())
	{
		// Write some info to a file for plotting
		ofstream fcurve("Analysis.out");
		ostream_iterator<string> os(fcurve, "\t ");
		*os++ = "Radius";
		*os++ = "v Background";
		*os++ = "v Increment";
		*os++ = "v Analysis";
		*os++ = "u Background";
		*os++ = "u Increment";
		*os++ = "u Analysis";
		*os++ = "w Background";
		*os++ = "w Increment";
		*os++ = "w Analysis";
		*os++ = "hstar Background";
		*os++ = "hstar Increment";
		*os++ = "hstar Analysis";
		*os++ = "qv Background";
		*os++ = "qv Increment";
		*os++ = "qv Analysis";
		*os++ = "rho prime Background";
		*os++ = "rho prime Increment";
		*os++ = "rho prime Analysis";
		
		fcurve << endl;
		ostream_iterator<double> of(fcurve, "\t ");
		double x = bgSpline->Xmin();
		double xs = (bgSpline->Xmax() - x) / 2000.0;
		for (; x <= bgSpline->Xmax(); x += xs)
		{
			*of++ = x;
			for (unsigned int var = 0; var < numVars; var++) {
				bgSpline->solve(&BGsave[var].front());
				analysisSpline->solve(&BG[var].front());
				double bg = bgSpline->evaluate (x);
				double a = analysisSpline->evaluate (x);
				*of++ = bg;
				*of++ = a-bg;
				*of++ = a;
			}		   
			fcurve << endl;
		}
	} else {
		cerr << "Analysis has a problem :(" << endl;
		return false;
	}

	
	// Get the temperature, pressure, and E
	bgSpline->solve(&BG[3].front());
	analysisSpline->solve(&BG[5].front());
	vector<real> temp;
	vector<real> press;
	vector<real> h;
	vector<real> hstarfwd;
	for (unsigned int T = 235; T < 315; T++) {
		temp.push_back(T);
		hstarfwd.push_back(0.0);
	}
	SplineD* tSpline = new SplineD(&temp.front(), temp.size(), &hstarfwd.front(), wl, bc, temp.size());
	double alt =2000;
	temp.clear();
	for (unsigned int ri = 0; ri < r.size(); ri++) {
		double rad = r[ri];
		double hstar = bgSpline->evaluate(rad);
		hstarfwd.clear();
		double rhoPrime = analysisSpline->evaluate(rad);
		double rho = 1.1646*exp(-1.068e-4*alt) + rhoPrime;
		for (unsigned int T = 235; T < 315; T++) {
			double es = (6.1078 * exp (17.2694 * (1.0 - 237.3 / (T - 273.15 + 237.3))));
			double qs = (1000.0 * 0.622 * es / (2.87*T*rho - es));
			double hs = 1005.7*T +  2.5e3*qs + 9.81*alt;
			hstarfwd.push_back(hstar*1000 - hs);
		}
		tSpline->solve(&hstarfwd.front());

		// Newton's method initial guess is from hydrostatics
		//double T =  ((9.81 * 1.1646*exp(-1.068e-4*alt) / 1.068e-4) - 5460.)/(rho*287.);
		double T = 280;
		double hmin = 1e34;
		while (fabs(hmin) > 0.01) {
			double h = tSpline->evaluate(T);
			double hprime = tSpline->slope(T);
			if (hprime != 0) {
				T = T - h/hprime;
				hmin = h;
			}
		}
		
		// Now push back temp, press, and h
		temp.push_back(T);
		press.push_back(T*rho*287./100);
	}

	ofstream fcurve("Coefficients.out");
	ostream_iterator<string> os(fcurve, "\t ");
	*os++ = "r";
	*os++ = "radius";
	*os++ = "rhov";
	*os++ = "rhou";
	*os++ = "rhow";
	*os++ = "hstar";
	*os++ = "qv";
	*os++ = "rhoprime";
	*os++ = "temperature";
	*os++ = "pressure";
		
	fcurve << endl;
	ostream_iterator<double> of(fcurve, "\t ");
	for (unsigned int ri = 0; ri < r.size(); ri++) {
		double rad = r[ri];
		*of++ = ri;
		*of++ = rad;
		for (unsigned int var = 0; var < numVars; var++) {
			analysisSpline->solve(&BG[var].front());
			*of++ = analysisSpline->getCoefficient(ri);
		}
		analysisSpline->solve(&temp.front());
		*of++ = analysisSpline->getCoefficient(ri);
		analysisSpline->solve(&press.front());
		*of++ = analysisSpline->getCoefficient(ri);
		fcurve << endl;	
	}
	
	delete analysisSpline;
	delete RXformSpline;
	cost1d->finalize();
	delete cost1d;
	
	return true;
}

bool VarDriver1d::initParaboloid()
{
	
	costAnalytic = CostFunctionAnalytic();
	costAnalytic.initialize();
	return true;
}

void VarDriver1d::runParaboloid()
{
	costAnalytic.minimize();
}

bool VarDriver1d::finalizeParaboloid()
{
	
	costAnalytic.finalize();
	return true;
}

double VarDriver1d::updateXforms()
{
	double wl = 2; 
	int num_nodes = r.size();
	double *Cq = new double[num_nodes*numVars];
	cost1d->getCq(Cq);
	double RMS = 0;
	vector<real> incr; 
	// Need to update rho with perturbation
	double rhoBar = 1.1646*exp(-1.068e-4*2000);
	for (unsigned int p = 0; p < r.size(); p++) {
		incr.push_back(0.0);
	}
	
	SplineD* incrSpline = new SplineD(&r.front(), r.size(), &incr.front(), wl, bc, num_nodes);
	if (!incrSpline->ok())
	{
		cerr << "Analysis has a problem :(" << endl;
		return false;
	} 
	
	// Increment the variables
	for (unsigned int var = 0; var < numVars; var++) {
		incr.clear();
		BG[var].clear();
		for (unsigned int p = 0; p < r.size(); p++) {
			unsigned int i = p + var*r.size();
			incr.push_back(Cq[i]);
		}
		bgSpline->solve(&BG[var].front());
		incrSpline->solve(&incr.front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			double bg = bgSpline->evaluate (rad);
			double cq = incrSpline->evaluate (rad);
			// cout << "Incr:" << ri << "\t" << bg << "\t" << cq << endl;
			double newbg = bg + cq;
			BG[var].push_back(newbg);
		}
	}
	
	// Compute an updated transform
	bgSpline->solve(&BG[0].front());
	incrSpline->solve(&BG[5].front());
	RXform.clear();
	for (unsigned int ri = 0; ri < r.size(); ri++) {
		double rad = r[ri];
		double rhov = bgSpline->evaluate (rad);
		double rp = incrSpline->evaluate (rad);
		double newv = rhov/(rhoBar + rp);
		double potRad = sqrt((rad*newv*1e3 + CoriolisF*rad*rad*1e6/2)*2/CoriolisF)/1e3;
		RXform.push_back(potRad);
	}
	RXformSpline->setDomain(&RXform.front(), r.size(), wl, bc, num_nodes);
	RXformSpline->solve(&r.front());
	rXform.clear();	
	for (unsigned int Ri = 0; Ri < R.size(); Ri++) {
		double potRad = R[Ri];
		double rhov = bgSpline->evaluate(RXformSpline->evaluate(potRad));
		double rp = incrSpline->evaluate(RXformSpline->evaluate(potRad));
		double newv = rhov/(rhoBar + rp);
		double physrad = (-newv + sqrt(newv*newv + CoriolisF*CoriolisF*potRad*potRad*1e6))/(CoriolisF*1e3);
		//cout << "r/R: " << physrad << "\t" << potRad << endl;
		rXform.push_back(physrad);
	}
	
	delete incrSpline;

	// Compute the RMS of the increment vector
	for (unsigned int var = 0; var < numVars; var++) {
		RMS = 0;
		for (unsigned int p = 0; p < r.size(); p++) {
			unsigned int i = p + var*r.size();
			RMS += Cq[i]*Cq[i];
		}
		RMS = sqrt(RMS/r.size());
		cout << "Variable RMS = " << var << "\t" << RMS << endl;
	}
	
	RMS = 0;
	for (unsigned int ri = 0; ri < numVars*r.size(); ri++) {
		RMS += Cq[ri]*Cq[ri];
		//cout << "Cq:\t" << r[ri] << "\t" << Cq[ri] << endl;
	}
	RMS = sqrt(RMS/(numVars*r.size()));
	cout << "Cq RMS = " << RMS << endl;
	delete[] Cq;
	
	return RMS;
	
}	

void VarDriver1d::EvalSpline (SplineD* spline, ostream &out)
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

void VarDriver1d::bgOut() {
	// And finally write the curve to a file
	ofstream fcurve("Background.out");
	ostream_iterator<string> os(fcurve, "\t ");
	*os++ = "Radius";
	*os++ = "v Background";
	*os++ = "u Background";
	
	fcurve << endl;
	ostream_iterator<double> of(fcurve, "\t ");
	double x = bgSpline->Xmin();
	double xs = (bgSpline->Xmax() - x) / 2000.0;
	for (; x <= bgSpline->Xmax(); x += xs)
	{
		*of++ = x;
		for (unsigned int var = 0; var < numVars; var++) {
			bgSpline->solve(&BGsave[var].front());
			double bg = bgSpline->evaluate (x);
			*of++ = bg;
		}		   
		fcurve << endl;
	}
}


void VarDriver1d::processMetObs()
{

	// Use the background rho prime to correct radar velocities
	bgSpline->solve(&BG[5].front());
	
	// Check the data directory for files
	QDir dataPath("/Users/mbell/Development/tcvar/vardata");
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
			if ((obTime < startTime) and (obTime > endTime)) continue;
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
				(abs(metOb.getAltitude() - 2000) > 200))
				continue;

			varOb.setRadius(rad);
			float Um = tcVector[tci].getUmean();
			float Vm = tcVector[tci].getVmean();
			varOb.setAltitude(metOb.getAltitude());
			double rhoBar = 1.1646*exp(-1.068e-4*metOb.getAltitude());
			double rhoPrimeBG = bgSpline->evaluate(rad);

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
					varOb.setOb((rho+rhoPrimeBG)*Vdopp);
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

