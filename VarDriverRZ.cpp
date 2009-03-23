/*
 *  VarDriverRZ.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriverRZ.h"
#include "Dorade.h"
#include <iterator>
#include <fstream>
#include <cmath>
#include <QTextStream>
#include <QFile>
#include "RecursiveFilter.h"

VarDriverRZ::VarDriverRZ()
	: VarDriver()
{
	numVars = 6;
	maxHeights = 50; // Can I make this dynamic?
	rhoBase = 1.1646;
	rhoInvScaleHeight = 1.068e-4;
	zincr = 250;
	rincr = 1;
	maxIter = 1;
	CQTOL = 0.5;
}

VarDriverRZ::~VarDriverRZ()
{
	for (unsigned int zi = 0; zi < maxHeights; zi++) {
		delete[] BG[zi];
		delete[] BGsave[zi];
	}
	delete[] BG;
	delete[] BGsave;
	delete[] scalarSpline;
	delete[] vecSpline;
	//delete[] ctrlSpline;
	delete zSpline;
	delete zSplinePsi;
	//delete[] RnumGridpts;
	delete[] RXform;
	delete[] obs;
	delete[] bgB;
	delete costRZ;

}


void VarDriverRZ::preProcessMetObs()
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
			
			// Make sure the ob is within the time limits
			MetObs metOb = metData->at(i);
			QDateTime obTime = metOb.getTime();
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
			double rhoBar = rhoBase*exp(-rhoInvScaleHeight*metOb.getAltitude());
			double rhoPrimeBG = zSpline->evaluate(metOb.getAltitude());
			
			// Initialize the weights
			varOb.setWeight(0., 0);
			varOb.setWeight(0., 1);
			varOb.setWeight(0., 2);
			varOb.setWeight(0., 3);
			varOb.setWeight(0., 4);
			varOb.setWeight(0., 5);
			double u, v, w, rho, qv, energy, rhov, rhou, rhow; 
			switch (metOb.getObType()) {
				case (MetObs::dropsonde):
					varOb.setType(MetObs::dropsonde);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					w = metOb.getVerticalVelocity();
					rho = metOb.getDryDensity();
					qv = metOb.getQv();
					energy = metOb.getTotalEnergy();
					
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
					if ((energy != -999) and (rho != -999)) {
						// energy 5 kJ error
						varOb.setWeight(1., 3);
						varOb.setOb(rho*energy/1e3);
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
					varOb.setType(MetObs::flightlevel);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					w = metOb.getVerticalVelocity();
					rho = metOb.getDryDensity();
					qv = metOb.getQv();
					energy = metOb.getTotalEnergy();
					
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
					if ((energy != -999) and (rho != -999)) {
						// energy 5 kJ error
						varOb.setWeight(1., 3);
						varOb.setOb(rho*energy/1e3);
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
					varOb.setType(MetObs::radar);
					// Geometry terms
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
					
					/* density correction term (rhoo/rho)*0.45 [rho(Z)=rho_o exp-(z/H), where 
					 H is the scale height = 9.58125 from Gray's inner 2 deg composite] 
					 0.45 density correction from Beard (1985, JOAT pp 468-471) 
					 Adjusted to use Jordan hydrostatic scale height -MB */
					double DCOR=exp(0.45*metOb.getAltitude()*0.0001068);
					
					// The snow relationship (Atlas et al., 1973) --- VT=0.817*Z**0.063  (m/s) 
					double VTS=-DCOR * (0.817*pow(ZZ,(double)0.063));
					
					// The rain relationship (Joss and Waldvogel,1971) --- VT=2.6*Z**.107 (m/s) */
					double VTR=-DCOR * (2.6*pow(ZZ,(double).107));
					/* Test if height is in the transition region between SNOW and RAIN
					   defined as hlow in km < H < hhi in km
					   if in the transition region do a linear weight of VTR and VTS */
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
					
					// Theoretically, rhoPrime could be included as a prognostic variable here...
					// However, adding another unknown without an extra equation makes the problem even more underdetermined
					// so assume it is small and ignore it
					// double rhopWgt = -Vdopp;
					//varOb.setWeight(rhopWgt, 5);
					
					// Set the error according to the spectrum width and potential fall speed error (assume 2 m/s?)
					double DopplerError = metOb.getSpectrumWidth() + wWgt*2.;
					varOb.setError(DopplerError);
					varOb.setOb((rhoBar+rhoPrimeBG)*Vdopp);
					obVector.push_back(varOb);
					
					break;
			}
		} 
		cout << obVector.size() << " total observations." << endl;
	}
	
	delete metData;
	
	// Write the Obs to a summary text file
	ofstream obstream("Observations.out");
	ostream_iterator<string> os(obstream, "\t ");
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
	obstream << endl;

	ostream_iterator<double> od(obstream, "\t ");
	for (unsigned int i=0; i < obVector.size(); i++) {
		Observation ob = obVector.at(i);
		*od++ = ob.getType();
		*od++ = ob.getRadius();
		*od++ = ob.getAltitude();
		// NULL 3rd dimension
		*od++ = -999.;
		*od++ = ob.getOb();
		*od++ = ob.getInverseError();
		for (unsigned int var = 0; var < numVars; var++)
			*od++ = ob.getWeight(var);

		obstream << endl;	
	}
	
	// Load the observations into a vector
	obs = new real[obVector.size()*10];
	for (unsigned int m=0; m < obVector.size(); m++) {
		int n = m*10;
		Observation ob = obVector.at(m);
		obs[n] = ob.getOb();
		obs[n+1] = ob.getInverseError();
		for (unsigned int var = 0; var < numVars; var++) {
			obs[n+2+var] = ob.getWeight(var);
		}
		obs[n+2+numVars] = ob.getRadius();
		obs[n+3+numVars] = ob.getAltitude();
	}	
	
	// All done preprocessing
	if (!processedFiles) {
		cout << "No files processed, nothing to do :(" << endl;
		// return 0;
	} else {
		cout << "Finished preprocessing " << processedFiles << " files." << endl;
	}
	
}

bool VarDriverRZ::loadMetObs()
{

	// Our generic observation
	Observation varOb;
	double wgt[numVars];
	double radius, alt, nulldim, ob, error;
	int type;
	ifstream obstream("./Observations.out");
	while (obstream >> type >> radius >> alt >> nulldim >> ob >> error
		   >> wgt[0] >> wgt[1] >> wgt[2] >> wgt[3] >> wgt[4] >> wgt[5])
	{
		varOb.setOb(ob);
		varOb.setRadius(radius);
		varOb.setAltitude(alt);
		varOb.setType(type);
		varOb.setError(1./error);
		for (unsigned int var = 0; var < numVars; var++)
			varOb.setWeight(wgt[var],var);
		obVector.push_back(varOb);
	}

	// Load the observations into a vector
	obs = new real[obVector.size()*10];
	for (unsigned int m=0; m < obVector.size(); m++) {
		int n = m*10;
		Observation ob = obVector.at(m);
		obs[n] = ob.getOb();
		obs[n+1] = ob.getInverseError();
		for (unsigned int var = 0; var < numVars; var++) {
			obs[n+2+var] = ob.getWeight(var);
		}
		obs[n+2+numVars] = ob.getRadius();
		obs[n+3+numVars] = ob.getAltitude();
	}	
	
	return true;
}

bool VarDriverRZ::loadBGfromFile()
{
	
	// Read in the background state
	// Read the r and v pairs from a file
	double height, radius, v, psi, rhoe, q, rho;
	int zi = 0;
	vector<real> vIn, vBG, psiBG, rhoeBG, qBG, rpBG;
	vector<real>* vInit = new vector<real>[maxHeights];
	ifstream vdata("./RZbackground.txt");
	vdata.width(14);
	while (vdata >> height >> radius >> v >> psi >> rhoe >> q >> rho)
	{
		if (z.empty()) z.push_back (height);
		if (height != z.back()) {
			// Assign the initial background fields
			BG[zi][0] =  vBG;
			BG[zi][1] =  psiBG;
			BG[zi][2] =  rhoeBG;
			BG[zi][3] =  qBG;
			BG[zi][4] =  rpBG;
			BG[zi][5] =  vIn;
			r.clear(); vIn.clear();
			vBG.clear(); psiBG.clear();
			rhoeBG.clear(); qBG.clear(); rpBG.clear();
			z.push_back (height);
			zi++;
		}		
		r.push_back (radius);
		
		// No negative v for potential radius transform!!!
		// Need to handle this more gracefully
		if (v < 0)  v = 0;
		vIn.push_back (v);
		vBG.push_back (v*rho);
		psiBG.push_back (psi/1e6);
		rhoeBG.push_back (rhoe/1e3);
		qBG.push_back (q);
		rpBG.push_back (rho-rhoBase*exp(-rhoInvScaleHeight*height));
	}
	
	// Assign the final strip
	BG[zi][0] =  vBG;
	BG[zi][1] =  psiBG;
	BG[zi][2] =  rhoeBG;
	BG[zi][3] =  qBG;
	BG[zi][4] =  rpBG;
	BG[zi][5] =  vIn;
	
	delete[] vInit;
	
	// Check that z.size is not bigger than allocated array
	if (z.size() > maxHeights) {
		cerr << "Memory overflow in z direction :" << z.size() << ">\t" << maxHeights << endl;
		return false;
	} else if (!z.size()) {
		cerr << "No heights! Problem reading BG file" << endl;
	}
	
	return true;
	
}	

bool VarDriverRZ::setupMishAndRXform()
{
	// Set up a vertical spline on the original grid
	double wl = 2*zincr; 
	bc = SplineBase::BC_ZERO_SECOND;
	int numZnodes = z.size();
	vector<real> zBuffer;
	zBuffer.assign(z.size(), 0.0);	
	zSpline = new SplineD(&z.front(), z.size(), &zBuffer.front(), wl, SplineBase::BC_ZERO_SECOND, numZnodes);
	zSplinePsi = new SplineD(&z.front(), z.size(), &zBuffer.front(), wl, SplineBase::BC_LZERO_RSECOND, numZnodes);
	cout << "Vertical Spline Cutoff frequency " << wl
	<< ", number of nodes " << numZnodes
	<< ", and boundary condition type " << bc << "\n";
	if (!zSpline->ok()) {
		cout << "Vertical spline has a problem :(" << endl;
		return false;
	}
	
	// Num levels set to original grid first
	numHeights = z.size();
	
	///Realign the vertical grid on the Gaussian points
	double zMin = 0.;
	double zMax = z.back();
	z.clear();
	numZnodes = (int)((zMax - zMin)/zincr) + 1;
	int zinit = 1;
	for (unsigned int vi = 0; vi< numVars; vi++) {
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			for (unsigned int zi = 0; zi < numHeights; zi++) {
				zBuffer[zi] = BG[zi][vi].at(ri);
			}
			SplineD* bgSpline;
			if (vi == 1) {
				bgSpline = zSplinePsi;
			} else {
				bgSpline = zSpline;
			}
			bgSpline->solve(&zBuffer.front());
			unsigned int zi = 0;
			for (double height = zMin; height < zMax; height += zincr) {
				for (int mu = -1; mu <= 1; mu += 2) {
					double gheight = height + (0.5*sqrt(1./3.)* mu + 0.5) * zincr;
					// cout << gheight << "\t" << zSpline->evaluate(gheight) << endl;
					if (zi >= maxHeights) {
						cerr << "Memory overflow in z direction: " << z.size() << ">\t" << numHeights << endl;
						return false;
					}
					if (zi < numHeights) {
						BG[zi][vi].at(ri) = bgSpline->evaluate(gheight);
					} else {
						BG[zi][vi].push_back(bgSpline->evaluate(gheight));
					}
					if (zinit) z.push_back(gheight);
					zi++;
					
				}
			}
			zinit = 0;
		}
	}
	
	// Expand spline and levels for Gaussian grid
	zSpline->setDomainGQ(&z.front(), z.size(), wl, bc, numZnodes);
	numHeights = z.size();
	ja = zSpline->getQfactored();
	
	// Set up a radial spline on the original grid
	wl = 2*rincr;
	vecSpline = new SplineD[numHeights];
	scalarSpline = new SplineD[numHeights];
	int numrnodes = r.size();
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		vecSpline[zi].setDomain(&r.front(), r.size(), wl, SplineBase::BC_LZERO_RSECOND, numrnodes);
		scalarSpline[zi].setDomain(&r.front(), r.size(), wl, SplineBase::BC_ZERO_SECOND, numrnodes);
	}
	
	// Iterate through the new radius grid for each variable
	double rMin = 0.;
	double rMax = (double)r.back();
	numrnodes = (int)((rMax - rMin)/rincr) + 1;
	r.clear();
	// And set up the transform array for physical to potential radius
	vector<real> parmWinds;	
	RXform = new vector<real>[numHeights];
	for (unsigned int vi = 0; vi< numVars; vi++) {
		SplineD* bgSpline;
		if ((vi > 1) and (vi < 5)) {
			bgSpline = scalarSpline;
		} else {
			bgSpline = vecSpline;
		}		
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			bgSpline[zi].solve(&BG[zi][vi].front());
			for (double rad = rMin; rad< rMax; rad+= rincr)
			{
				// Gaussian quadrature points
				for (int mu = -1; mu <= 1; mu += 2) {
					double grad = rad + 0.5*sqrt(1./3.)* mu + 0.5;
					BGsave[zi][vi].push_back(bgSpline[zi].evaluate(grad));
					//if (vi == 2) 
					//	cout << zi << "\t" << BGsave[zi][vi].back() << endl;
					if (vi == 5) { // Raw tangential wind for r->R xform;
						double vRaw = BGsave[zi][vi].back();
						double potRad = (grad*vRaw*1e3 + CoriolisF*grad*grad*1e6/2)*2/CoriolisF;
						if (potRad < 0) {
							// Inertially unstable!
							cout << "Inertially unstable background field! Removing instability" << endl;
							cout << zi << "\t" << rad << "\t" << mu << "\t" << vRaw << "\n";
							potRad = 0;
						} else {
							potRad = sqrt(potRad)/1e3;
						}
						if (zi == 0) r.push_back(grad);
						RXform[zi].push_back(potRad);
					}
				}
			}
			BG[zi][vi].clear();
			BG[zi][vi] = BGsave[zi][vi];
		}
	}
	
	// Set up the background vector in Gaussian r space
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		vecSpline[zi].setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_LZERO_RSECOND, numrnodes);
		scalarSpline[zi].setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_ZERO_SECOND, numrnodes);
	}
	ia = scalarSpline[0].getQfactored();
	
	cout << "Radial Spline Cutoff frequency " << wl
	<< ", number of nodes " << numrnodes
	<< ", and boundary condition type " << bc << "\n";
	
	// Get the nodal representation
	bgB = new real[numrnodes*numZnodes*(numVars-1)];
	real** bgrz = new real*[numrnodes];
	for (int rb = 0; rb < numrnodes; rb++) {
		bgrz[rb] = new real[z.size()];
	}
	
	for (unsigned int vi = 0; vi< (numVars-1); vi++) {
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			const real* bgr = scalarSpline[zi].getBGQ(&BG[zi][vi].front());
			for (int rb = 0; rb < numrnodes; rb++) {
				bgrz[rb][zi] = bgr[rb];
			}
		}
		for (int rb = 0; rb < numrnodes; rb++) {
			const real* bgz = zSpline->getBGQ(bgrz[rb]);
			for (int zb = 0; zb < numZnodes; zb++) {
				bgB[(numVars-1)*numrnodes*zb +(numVars-1)*rb + vi] = bgz[zb];
			}
		}
	}	

	// Make sure the saved field is numerically consistent with the bgB's
	real* bgz = new real[numZnodes];
	real* bgr = new real[numrnodes];
	for (unsigned int vi = 0; vi< (numVars-1); vi++) {
		for (int rb = 0; rb < numrnodes; rb++) {
			for (int zb = 0; zb < numZnodes; zb++) {
				bgz[zb] = bgB[(numVars-1)*numrnodes*zb +(numVars-1)*rb + vi];
			}
			if (vi == 1) {
				zSplinePsi->solveBGQ(bgz);
				for (unsigned int zi = 0; zi < z.size(); zi++) {
					bgrz[rb][zi] = zSplinePsi->evaluate(z[zi]);
				}
			} else {
				zSpline->solveBGQ(bgz);
				for (unsigned int zi = 0; zi < z.size(); zi++) {
					bgrz[rb][zi] = zSpline->evaluate(z[zi]);
				}
			}
		}
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			for (int rb = 0; rb < numrnodes; rb++) {
				bgr[rb] = bgrz[rb][zi];
			}
			if (vi <= 1) {
				if (!vecSpline[zi].solveBGQ(bgr)) cout << "BGQ Failed!\n";
				for (unsigned int ri = 0; ri < r.size(); ri++) {
					BGsave[zi][vi].at(ri) = vecSpline[zi].evaluate(r[ri]);
				}
			} else {
				if (!scalarSpline[zi].solveBGQ(bgr)) cout << "BGQ Failed!\n";
				for (unsigned int ri = 0; ri < r.size(); ri++) {
					BGsave[zi][vi].at(ri) = scalarSpline[zi].evaluate(r[ri]);
				}
			}
		}
	}	
	
	for (int rb = 0; rb < numrnodes; rb++) {
		delete[] bgrz[rb];
	}
	delete[] bgrz;
	delete[] bgz;
	delete[] bgr;
	
	
	/* Set up Anisotropic Filter
	real* tau = new real[numrnodes*numZnodes];
	real* Rtmp = new real[z.size()];
	real* bgr = new real[numrnodes];
	
	for (int rb = 0; rb < numrnodes; rb++) {
		bgr[rb] = bgB[(numVars-1)*rb];
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			vecSpline[zi].solveGQ(&BG[zi][5].front());
			double rad = rMin + rb*rincr;
			double vRaw = vecSpline[zi].evaluate(rad);
			double potRad = (rad*vRaw*1e3 + CoriolisF*rad*rad*1e6/2)*2/CoriolisF;
			Rtmp[zi] = sqrt(potRad)/1e3;
		}
		zSpline->solve(Rtmp);
		for (int zb = 0; zb < numZnodes; zb++) {
			tau[numrnodes*zb +rb] = zSpline->evaluate(zMin + zb*zincr);
		}
	}
	
	real* rtau = new real[numrnodes];
	rtau[0] = tau[1]-tau[0];
	for (int rb = 1; rb < numrnodes; rb++) {
		rtau[rb] = tau[rb]-tau[rb-1];
		cout << rb << "\t" << tau[rb] << "\t" << rtau[rb] << endl;
	}
	RecursiveFilter* anifilter = new RecursiveFilter(5, rtau, numrnodes);
	RecursiveFilter* isofilter = new RecursiveFilter(4, 5);
	isofilter->filterArray(bgr, numrnodes);
	for (int rb = 0; rb < numrnodes; rb++) {
		cout << rb << "\t" << bgr[rb] << "\n";
		bgr[rb] = bgB[(numVars-1)*rb];
	}
	anifilter->aniFilterArray(bgr,numrnodes);
	for (int rb = 0; rb < numrnodes; rb++) {
		cout << rb << "\t" << bgr[rb] << "\n";
	}
	
	delete[] bgr;
	
	// Old formulation
	ctrlSpline = new SplineD[numZnodes];
	rXform = new vector<real>[numHeights];
	RnumGridpts = new unsigned int[numHeights];
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		scalarSpline[zi].solveGQ(&RXform[zi].front());
		for (int rb = 0; rb < numrnodes; rb++) {
			real rad = rMin	+ rb*rincr;
			rXform[zi].push_back(scalarSpline[zi].evaluate(rad));
		}
		double Rincr = 1;
		double minR = 0;
		double maxR = (int)RXform[zi].back() + 1;
		double Rwidth = (int)((maxR - minR)/Rincr) + 1;
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
		
	} */
	
	// Set the master dimensions
	imin = rMin;
	imax = rMax;
	idim = numrnodes;
	jmin = zMin;
	jmax = zMax;
	jdim = numZnodes;
	
	return true;
}

bool VarDriverRZ::initialize()
{
	// Run a RZ vortex background field
	cout << "Initializing Vortex Background" << endl;

	// Allocate memory for the BG fields
	BG = new vector<real>*[maxHeights];
	BGsave = new vector<real>*[maxHeights];
	for (unsigned int zi = 0; zi < maxHeights; zi++) {
		BG[zi] = new vector<real>[numVars];
		BGsave[zi] = new vector<real>[numVars];
	}
	
	/* Load a BG field (on a physical grid) from a file, or
	   construct a parametric field, or (eventually) get from ESMF Coupler */
	bool loadBG = true;
	if (loadBG) {
		loadBGfromFile();
	} else {
		// Not implemented yet
		// Get the parametric winds from the input winds if data does not extend through domain
		//ParametricVortex parmVortex (&BG[zi][5].front(), &r[0],  r.size());
		//double vParm = parmVortex.getWindAtRadiusKM(grad);
	}
	
	// Set up the splines on the Gaussian Grid
	setupMishAndRXform();
	
	// Read in the TC centers
	// Ideally, create a time-based spline from limited center fixes here
	// but just load 1 second centers into vector for now
	readTCcenters();
	
	// Read in the observations, process them into weights and positions
	// Either preprocess from raw observations or load an already processed Observations.out file
	bool preprocess = true;
	if (preprocess) {
		preProcessMetObs();
	} else {
		loadMetObs();
	}
	cout << "Number of Observations: " << obVector.size() << endl;
	
	// Define the sizes of the arrays we are passing to the cost function
	int stateSize = idim*jdim*(numVars-1);
	
	costRZ = new CostFunctionRZ_CPU(obVector.size(), stateSize);
	costRZ->initialize(imin, imax, idim, jmin, jmax, jdim, ia, ja, bgB, obs, vecSpline, RnumGridpts, RXform); 
	
	return true;
}

bool VarDriverRZ::run()
{

	double CQRMS = 999;
	int iter=0;
	while ((CQRMS > CQTOL) and (iter < maxIter)) {
		iter++;
		cout << "Outer Loop Iteration: " << iter << endl;
		costRZ->minimize();
		CQRMS = updateXforms();
		costRZ->initState();
	}	
	cout << "Increment RMS Tolerance of " << CQTOL << " reached in "
		<< iter << " iterations. Writing analysis results..." << endl;

	return true;

}


bool VarDriverRZ::finalize()
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

	
	// Recalculate the BG fields on the mish
	int numrnodes = vecSpline[0].nNodes();
	int numZnodes = zSpline->nNodes();
	real** bgrz = new real*[numrnodes];
	for (int rb = 0; rb < numrnodes; rb++) {
		bgrz[rb] = new real[z.size()];
	}
	real* bgz = new real[numZnodes];
	real* bgr = new real[numrnodes];
	
	for (unsigned int vi = 0; vi< (numVars-1); vi++) {
		for (int rb = 0; rb < numrnodes; rb++) {
			for (int zb = 0; zb < numZnodes; zb++) {
				bgz[zb] = bgB[(numVars-1)*numrnodes*zb +(numVars-1)*rb + vi];
			}
			if (vi == 1) {
				zSplinePsi->solveBGQ(bgz);
				for (unsigned int zi = 0; zi < z.size(); zi++) {
					bgrz[rb][zi] = zSplinePsi->evaluate(z[zi]);
				}
			} else {
				zSpline->solveBGQ(bgz);
				for (unsigned int zi = 0; zi < z.size(); zi++) {
					bgrz[rb][zi] = zSpline->evaluate(z[zi]);
				}
			}
		}
		for (unsigned int zi = 0; zi < z.size(); zi++) {
			for (int rb = 0; rb < numrnodes; rb++) {
				bgr[rb] = bgrz[rb][zi];
			}
			if (vi <= 1) {
				vecSpline[zi].solveBGQ(bgr);
				for (unsigned int ri = 0; ri < r.size(); ri++) {
					BG[zi][vi].at(ri) = vecSpline[zi].evaluate(r[ri]);
				}
			} else {
				scalarSpline[zi].solveBGQ(bgr);
				for (unsigned int ri = 0; ri < r.size(); ri++) {
					BG[zi][vi].at(ri) = scalarSpline[zi].evaluate(r[ri]);
				}
			}
		}
	}	
	
	for (int rb = 0; rb < numrnodes; rb++) {
		delete[] bgrz[rb];
	}
	delete[] bgrz;
	delete[] bgz;
	delete[] bgr;
	
	vector<real>* rhoz = new vector<real>[numHeights];
	vector<real>* uz = new vector<real>[numHeights];
	vector<real>* vz = new vector<real>[numHeights];
	vector<real>* wz = new vector<real>[numHeights];
	vector<real>* temp = new vector<real>[numHeights];
	vector<real>* press = new vector<real>[numHeights];
	
	//RhoZ
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		scalarAnalysisSpline->solveGQ(&BG[zi][4].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double alt = z[zi];
			double rad = r[ri];
			double rhoPrime = scalarAnalysisSpline->evaluate(rad);
			double rho = rhoBase*exp(-rhoInvScaleHeight*alt) + rhoPrime;
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
		
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		// Get the temperature, pressure, and E
		double alt = z[zi];
		scalarSpline[zi].solveGQ(&BG[zi][2].front());
		scalarAnalysisSpline->solveGQ(&BG[zi][3].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			double energy = scalarSpline[zi].evaluate(rad)*1000;
			double qv = scalarAnalysisSpline->evaluate(rad);
			double t = (energy/rhoz[zi].at(ri) - 2.5e3*qv - 9.81*alt - (vz[zi].at(ri)*vz[zi].at(ri) + uz[zi].at(ri)*uz[zi].at(ri)))/1005.7;
			double p = t*rhoz[zi].at(ri)*287./100;
			// Now push back temp, press, and rho
			temp[zi].push_back(t);
			press[zi].push_back(p);
			
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

	// Write to an asi file for plotting
	QString asifile("tcvar_analysis");
	writeAsi(asifile, final);
	
	
	//Increments
	//RhoZ
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		scalarAnalysisSpline->solveGQ(&BGsave[zi][4].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double alt = z[zi];
			double rad = r[ri];
			double rhoPrime = scalarAnalysisSpline->evaluate(rad);
			double rho = rhoBase*exp(-rhoInvScaleHeight*alt) + rhoPrime;
			rhoz[zi].at(ri) = rho;
		}
	}
	
	// Winds
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		vecSpline[zi].solveGQ(&BGsave[zi][0].front());
		scalarAnalysisSpline->solveGQ(&rhoz[zi].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			double vwind = vecSpline[zi].evaluate(rad) / scalarAnalysisSpline->evaluate(rad);
			vz[zi].at(ri) = vwind;
		}
		vecSpline[zi].solveGQ(&BGsave[zi][1].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			if (rad != 0) {
				double wwind = 1e6 * (vecSpline[zi].slope(rad) / 1000) / (scalarAnalysisSpline->evaluate(rad) * rad * 1000);
				wz[zi].at(ri) = wwind;
			} else {
				wz[zi].at(ri) = 0.0;
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
			scalarAnalysisSpline->solveGQ(&rhoz[zi].front());
			if (rad != 0) {
				double uwind = 1e6 * -zSplinePsi->slope(alt) / (scalarAnalysisSpline->evaluate(rad) * rad * 1000);
				uz[zi].at(ri) = uwind;
			} else {
				uz[zi].at(ri) = 0.0;
			}
		}
	}
	
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		// Get the temperature, pressure, and E
		double alt = z[zi];
		scalarSpline[zi].solveGQ(&BGsave[zi][2].front());
		scalarAnalysisSpline->solveGQ(&BGsave[zi][3].front());
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			double rad = r[ri];
			double energy = scalarSpline[zi].evaluate(rad)*1000;
			double qv = scalarAnalysisSpline->evaluate(rad);
			double t = (energy/rhoz[zi].at(ri) - 2.5e3*qv - 9.81*alt - (vz[zi].at(ri)*vz[zi].at(ri) + uz[zi].at(ri)*uz[zi].at(ri)))/1005.7;
			double p = t*rhoz[zi].at(ri)*287./100;
			// Now push back temp, press, and rho
			temp[zi].at(ri) = t;
			press[zi].at(ri) = p;
			
		}
	}
	
	// Need to change this to B printout
	/* ofstream fcurve("Coefficients.out");
	ostream_iterator<string> os(fcurve, "\t ");
	*os++ = "z";
	*os++ = "r";
	*os++ = "radius";
	*os++ = "rhov";
	*os++ = "psi";
	*os++ = "energy";
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
	} */
	
	
	// Get the increments
	vector<real>** increment = new vector<real>*[numHeights];
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		increment[zi] = new vector<real>[11];
		for (unsigned int var = 0; var < 11; var++) {
			increment[zi][var].assign(r.size(), 0.0);
		}
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			for (unsigned int var = 0; var < numVars-1; var++) {
				increment[zi][var].at(ri) = final[zi][var].at(ri) - BGsave[zi][var].at(ri);
			}	 
			increment[zi][5].at(ri) = final[zi][5].at(ri) - temp[zi].at(ri);
			increment[zi][6].at(ri) = final[zi][6].at(ri) - press[zi].at(ri);
			increment[zi][7].at(ri) = final[zi][7].at(ri) - rhoz[zi].at(ri);
			increment[zi][8].at(ri) = final[zi][8].at(ri) - uz[zi].at(ri);
			increment[zi][9].at(ri) = final[zi][9].at(ri) - vz[zi].at(ri);
			increment[zi][10].at(ri) = final[zi][10].at(ri) - wz[zi].at(ri);
		}
	}
	
	
	
	// Write to an asi file for plotting
	asifile = "tcvar_increment";
	writeAsi(asifile, increment);
	
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
	*as++ = "energy Background";
	*as++ = "energy Increment";
	*as++ = "energy Analysis";
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
	for (unsigned int zi = 0; zi < z.size(); zi++) {
		for (unsigned int ri = 0; ri < r.size(); ri++) {
			*af++ = z[zi];
			*af++ = r[ri];
			for (unsigned int var = 0; var < 11; var++) {
				*af++ = final[zi][var].at(ri) - increment[zi][var].at(ri);
				*af++ = increment[zi][var].at(ri);
				*af++ = final[zi][var].at(ri);
			}			
			acurve << endl;
		}
	}
	
	delete scalarAnalysisSpline;
	delete vecAnalysisSpline;
	delete[] temp;
	delete[] press;
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
	costRZ->finalize();
	
	return true;
}

double VarDriverRZ::updateXforms()
{
	int wl = 2*int(rincr); 
	int numrnodes = vecSpline[0].nNodes();
	int numZnodes = zSpline->nNodes();
	double *Cq = new double[numrnodes*numZnodes*(numVars-1)];

	double RMS = 0;
	vector<real> incr(r.size(), 0.); 
	// Need to update rho with perturbation
	SplineD* vecIncrSpline = new SplineD(&r.front(), r.size(), &incr.front(), wl, SplineBase::BC_LZERO_RSECOND, numrnodes);
	vecIncrSpline->setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_LZERO_RSECOND, numrnodes);
	SplineD* scalarIncrSpline = new SplineD(&r.front(), r.size(), &incr.front(), wl, SplineBase::BC_ZERO_SECOND, numrnodes);
	scalarIncrSpline->setDomainGQ(&r.front(), r.size(), wl, SplineBase::BC_ZERO_SECOND, numrnodes);

	if (!vecIncrSpline->ok())
	{
		cerr << "Analysis has a problem :(" << endl;
		return false;
	} 
		
	/* Compute an updated transform
	double RxRMS = 0;
	for (unsigned int zi = 0; zi < z.size(); zi++) {	
		double rhoBar = rhoBase*exp(-rhoInvScaleHeight*z[zi]);
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
	RxRMS = sqrt(RxRMS/(r.size()*z.size())); */
	
	delete vecIncrSpline;
	delete scalarIncrSpline;
	
	// Compute the RMS of the increment vector
	double vMax = 0;
	double bgRMS = 0;
	
	for (unsigned int vi = 0; vi< (numVars-1); vi++) {
		for (int rb = 0; rb < numrnodes; rb++) {
			for (int zb = 0; zb < numZnodes; zb++) {
				Cq[(numVars-1)*numrnodes*zb +(numVars-1)*rb + vi] = bgB[(numVars-1)*numrnodes*zb +(numVars-1)*rb + vi];
			}
		}
	}
	
	// Increment the variables
	costRZ->updateBG();
	
	for (unsigned int vi = 0; vi < numVars-1; vi++) {
		RMS = 0;
		bgRMS = 0;
		double maxIncr = 0;
		double maxPercent = 0;
		for (int rb = 0; rb < numrnodes; rb++) {
			for (int zb = 0; zb < numZnodes; zb++) {
				unsigned int i = (numVars-1)*numrnodes*zb +(numVars-1)*rb + vi;
				double incr = bgB[i]- Cq[i];
				if (fabs(incr) > maxIncr) maxIncr = fabs(incr);
				
				if (bgB[i] != 0) {
					double percent = 100*fabs(incr/bgB[i]);
					if (percent > maxPercent) maxPercent = percent;
				}
				RMS += incr*incr;
				//cout << i << "\t" << p << "\t" << var << "\t" << Cq[i] << endl;
				bgRMS += bgB[i]*bgB[i];
			}
		}
		RMS = sqrt(RMS/(numrnodes*numZnodes));
		bgRMS = sqrt(bgRMS/(numrnodes*numZnodes));
		cout << "Variable " << vi << " RMS = " << bgRMS << endl;
		cout << "Max Increment = " << maxIncr << endl;
		cout << "Max % Increment = " << maxPercent << " %" << endl;
		cout << "RMS Increment = " << RMS << "\t( " << 100*RMS/bgRMS << " %)" << endl << endl;
		if (vi == 0) vMax = 100*RMS/bgRMS;
	}
	
	RMS = 0;
	for (unsigned int vi = 0; vi< (numVars-1); vi++) {
		for (int rb = 0; rb < numrnodes; rb++) {
			for (int zb = 0; zb < numZnodes; zb++) {
				unsigned int i = (numVars-1)*numrnodes*zb +(numVars-1)*rb + vi;
				double incr = bgB[i]- Cq[i];
				RMS += incr*incr;
			}
		}
	}
	RMS = sqrt(RMS/((numVars-1)*numrnodes*numZnodes));
	cout << "Cq (delta) RMS = " << RMS << endl;
	delete[] Cq;
	
	// Use Rx criteria
	//cout << "R transform RMS difference: " << RxRMS << endl << endl;
	//return RxRMS;
	// Use the maximum Vt as the criterion
	//return vMax;
	// Use total RMS (not very good measure I think)
	return RMS;
	
}	

void VarDriverRZ::EvalSpline (SplineD* spline, ostream &out)
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



bool VarDriverRZ::writeAsi(const QString& fileName, vector<real>** fields)
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
	fieldNames << "RV" << "SF" << "PE" << "QV" << "RP" << "T" << "P" << "RO" << "U" << "V" << "W";
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
	id[160] = (int)imin;
	id[161] = (int)imax*100;
	id[162] = (int)idim;
	id[163] = (int)(rincr * 1000);
	id[164] = 1;
	
	// Y Header
	double zincr = 0.25;
	id[165] = (int)(jmin);
	id[166] = (int)(jmax/10);
	id[167] = (int)jdim;
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
	id[310] = (int)(1);
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
	
	// Get the data on the nodes
	real* jTemp = new real[z.size()];
	real*** fieldNodes = new real**[fieldNames.size()];
	for(int n = 0; n < fieldNames.size(); n++) {
		fieldNodes[n] = new real*[idim];
		for (int i=0; i < idim; i++) {
			fieldNodes[n][i] = new real[jdim];
		}
	}
	
	for(int n = 0; n < fieldNames.size(); n++) {
		for (int i = 0; i < idim; i++) {
			for (unsigned int zi = 0; zi < z.size(); zi++) {
				const real* curve;
				if ((n > 1) and (n < 8)) {
					scalarSpline[zi].solveGQ(&fields[zi][n].front());
					curve = scalarSpline[zi].curve();
					jTemp[zi] = curve[i];
				} else {
					vecSpline[zi].solveGQ(&fields[zi][n].front());
					curve = vecSpline[zi].curve();
					jTemp[zi] = curve[i];
				}
			}
			const real* curve;
			if ((n > 1) and (n < 8)) {
				zSpline->solveGQ(jTemp);
				curve = zSpline->curve();
				for(int j = 0; j < jdim; j++) {
					fieldNodes[n][i][j] = curve[j];
				}
			} else {
				zSplinePsi->solveGQ(jTemp);
				curve = zSplinePsi->curve();
				for(int j = 0; j < jdim; j++) {
					fieldNodes[n][i][j] = curve[j];
				}
			}
		}
	}
	
	// Write data
	for(int k = 0; k < 1; k++) {
		out << reset << "level" << qSetFieldWidth(2) << k+1 << endl;
		for(int j = 0; j < jdim; j++) {
			out << reset << "azimuth" << qSetFieldWidth(3) << j+1 << endl;
			for(int n = 0; n < fieldNames.size(); n++) {
				out << reset << left << fieldNames.at(n) << endl;
				int line = 0;
				for (int i = 0; i < idim;  i++){
					out << reset << qSetRealNumberPrecision(3) << scientific << qSetFieldWidth(10) << fieldNodes[n][i][j];
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
	
	
	delete[] jTemp;
	for(int n = 0; n < fieldNames.size(); n++) {
		for (int i=0; i < idim; i++) {
			delete[] fieldNodes[n][i];
		}
		delete[] fieldNodes[n];
	}
	delete[] fieldNodes;
	
	return true;
	
}     

