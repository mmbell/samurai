/*
 *  VarDriverRZ.cpp
 *  samurai
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
	maxHeights = 80; // Can I make this dynamic?
	rhoBase = 1.156;
	rhoInvScaleHeight = 9.9504e-5;
	zincr = 100;
	rincr = 1;
	maxIter = 2;
	CQTOL = 0.5;
}

VarDriverRZ::~VarDriverRZ()
{

}

bool VarDriverRZ::finalize()
{
	for (unsigned int zi = 0; zi < maxHeights; zi++) {
		delete[] BG[zi];
		delete[] BGsave[zi];
	}
	delete[] BG;
	delete[] BGsave;
	//delete[] scalarSpline;
	//delete[] vecSpline;
	//delete[] ctrlSpline;
	//delete zSpline;
	//delete zSplinePsi;
	//delete[] RnumGridpts;
	//delete[] RXform;
	delete[] obs;
	//delete[] bgB;
	//delete[] ia;
	//delete[] ja;
	delete costRZ;
	
	return true;
}

void VarDriverRZ::preProcessMetObs()
{
	
	vector<real> rhoP;
		
	// Check the data directory for files
	QDir dataPath("./vardata");
	dataPath.setFilter(QDir::Files);
	dataPath.setSorting(QDir::Name);
	QStringList filenames = dataPath.entryList();

	QDateTime startTime = tcVector.front().getTime();
	QDateTime endTime = tcVector.back().getTime();
	QString start = startTime.toString(Qt::ISODate);
	QString end = endTime.toString(Qt::ISODate);
	cout << "TC centers found from " << start.toAscii().data() << " to " << end.toAscii().data() << endl;	
	
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
			case (sfmr):
				if (!read_sfmr(metFile, metData))
					cout << "Error reading sfmr file" << endl;
				break;
			case (wwind):
				if (!read_wwind(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (qscat):
				if (!read_qscat(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (ascat):
				if (!read_ascat(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (cen):
				continue;
			default:
				cout << "Unknown data type, skipping..." << endl;
				continue;
		}
		
		processedFiles++;
		
		// Process the metObs into Observations
		for (int i = 0; i < metData->size(); ++i) {
			
			// Make sure the ob is within the time limits
			MetObs metOb = metData->at(i);
			QDateTime obTime = metOb.getTime();
			//cout << "Obtime: " << obTime.toString(Qt::ISODate).toAscii().data() << endl;
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
			if ((rad < r.front()) or (rad > r.back()) or (rad == 0.0) or
				(metOb.getAltitude() < 0) or
				(metOb.getAltitude() > z.back()))
				continue;
			
			varOb.setRadius(rad);
			real Um = tcVector[tci].getUmean();
			real Vm = tcVector[tci].getVmean();

			real height = metOb.getAltitude();
			varOb.setAltitude(height);
			
			//Set the theta, but not used for this analysis
			real theta = 180.0 * atan2(y, x) / Pi;
			varOb.setTheta(theta);
			
			// Reference states			
			real rhoBar = rhoBase*exp(-rhoInvScaleHeight*height);
			//real qBar = 19.562 - 0.004066*height + 7.8168e-7*height*height;
			real qBar = 18.189 - 0.0060609*height + 8.37e-7*height*height - 4.6009e-11*height*height*height;
			if (qBar < 0) qBar = 0;
			//real hBar = 3.5e5;
			real hBar = 3.4633e5 - 10.751*height + 0.0021949*height*height - 1.7019e-7*height*height*height
				+ 4.858e-12*height*height*height*height;
			// Use bilinear interpolation here too for now, eventually probably a spline
			//real rhoaBG = bilinearField(rad, height, 4)/100. + rhoBar;
			//real qBG = bilinearField(rad, height, 3) + qBar;
			//real rhoBG = rhoaBG*(1+qBG/1000.);
			
			// Initialize the weights
			varOb.setWeight(0., 0);
			varOb.setWeight(0., 1);
			varOb.setWeight(0., 2);
			varOb.setWeight(0., 3);
			varOb.setWeight(0., 4);
			varOb.setWeight(0., 5);
			double u, v, w, rho, rhoa, qv, energy, rhov, rhou, rhow, wspd, vBG, uBG; 
			switch (metOb.getObType()) {
				case (MetObs::dropsonde):
					varOb.setType(MetObs::dropsonde);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					w = metOb.getVerticalVelocity();
					rho = metOb.getMoistDensity();
					rhoa = metOb.getAirDensity();
					qv = metOb.getQv();
					energy = metOb.getMoistStaticEnergy();
					
					// Separate obs for each measurement
					// rho v 1 m/s error
					if ((u != -999) and (rho != -999)) {
						varOb.setWeight(1., 0);
						rhov = rho*(-(u - Um)*y + (v-Vm)*x)/rad;
						varOb.setOb(rhov);
						varOb.setError(2.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						// rho u 1 m/s error
						varOb.setWeight(1., 1);
						rhou = rho*((u - Um)*x + (v-Vm)*y)/rad;
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(2.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);
					}
					if ((w != -999) and (rho != -999)) {
						// rho w 2 m/s error
						varOb.setWeight(1., 2);
						rhow = rho*w;
						varOb.setOb(rhow);
						varOb.setError(4.);
						obVector.push_back(varOb);
						varOb.setWeight(0., 2);
					}
					if (energy != -999) {
						// energy 5 kJ error
						varOb.setWeight(1., 3);
						varOb.setOb((energy - hBar)*1.e-3);
						varOb.setError(5.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 3);
					}
					if (qv != -999) {
						// Qv 2 g/kg error
						varOb.setWeight(1., 4);
						varOb.setOb(qv-qBar);
						varOb.setError(2.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 4);
					}
					if (rhoa != -999) {
						// Rho prime .1 kg/m^3 error
						varOb.setWeight(1., 5);
						varOb.setOb((rhoa-rhoBar)*100);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 5);
					}
					
					break;
					
				case (MetObs::flightlevel):
					varOb.setType(MetObs::flightlevel);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					w = metOb.getVerticalVelocity();
					rho = metOb.getMoistDensity();
					rhoa = metOb.getAirDensity();
					qv = metOb.getQv();
					energy = metOb.getMoistStaticEnergy();
					
					// Separate obs for each measurement
					// rho v 1 m/s error
					if ((u != -999) and (rho != -999)) {
						varOb.setWeight(1., 0);
						rhov = rho*(-(u - Um)*y + (v-Vm)*x)/rad;
						varOb.setOb(rhov);
						varOb.setError(2.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						// rho u 1 m/s error
						varOb.setWeight(1., 1);
						rhou = rho*((u - Um)*x + (v-Vm)*y)/rad;
						varOb.setOb(rhou);
						varOb.setError(2.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);
					}
					if ((w != -999) and (rho != -999)) {
						// rho w 1 dm/s error
						varOb.setWeight(1., 2);
						rhow = rho*w;
						varOb.setOb(rhow);
						varOb.setError(2.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 2);
					}
					if (energy != -999) {
						// energy 5 kJ error
						varOb.setWeight(1., 3);
						varOb.setOb((energy - hBar)*1.e-3);
						varOb.setError(5.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 3);
					}
					if (qv != -999) {
						// Qv 2 g/kg error
						varOb.setWeight(1., 4);
						varOb.setOb(qv-qBar);
						varOb.setError(2.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 4);
					}
					if (rhoa != -999) {
						// Rho prime .1 kg/m^3 error
						varOb.setWeight(1., 5);
						varOb.setOb((rhoa-rhoBar)*100);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 5);
					}
					
					break;

				case (MetObs::sfmr):
					varOb.setType(MetObs::sfmr);
					wspd = metOb.getWindSpeed();
					vBG = 1.e3*bilinearField(rad, height, 0)/rad;
					uBG = -1.e5*bilinearField(rad, 20., 1)/(rad*20.);
					varOb.setWeight(1., 0);
					//varOb.setWeight(1., 1);
					varOb.setOb(wspd);
					varOb.setError(10.0);
					obVector.push_back(varOb);
					break;
					
				case (MetObs::qscat):
					varOb.setType(MetObs::qscat);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					if (u != -999) {
						varOb.setWeight(1., 0);
						rhov = (-(u - Um)*y + (v-Vm)*x)/rad;
						varOb.setOb(rhov);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						// rho u 1 m/s error
						varOb.setWeight(1., 1);
						rhou = ((u - Um)*x + (v-Vm)*y)/rad;
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);
					}
					break;
					
				case (MetObs::ascat):
					varOb.setType(MetObs::ascat);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					if (u != -999) {
						varOb.setWeight(1., 0);
						rhov = (-(u - Um)*y + (v-Vm)*x)/rad;
						varOb.setOb(rhov);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						// rho u 1 m/s error
						varOb.setWeight(1., 1);
						rhou = ((u - Um)*x + (v-Vm)*y)/rad;
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);
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
					double DopplerError = metOb.getSpectrumWidth() + fabs(wWgt)*2.;
					if (DopplerError < 2.0) DopplerError = 2.0;
					varOb.setError(DopplerError);
					varOb.setOb(Vdopp);
					obVector.push_back(varOb);
					
					break;
					
			}

		} 
		cout << obVector.size() << " total observations." << endl;
	}
	
	delete metData;
	
	// Write the Obs to a summary text file
	ofstream obstream("Observations.in");
	// Header messes up reload
	/*ostream_iterator<string> os(obstream, "\t ");
	*os++ = "Type";
	*os++ = "r";
	*os++ = "z";
	*os++ = "theta";
	*os++ = "Observation";
	*os++ = "Inverse Error";
	*os++ = "Weight 1";
	*os++ = "Weight 2";
	*os++ = "Weight 3";
	*os++ = "Weight 4";
	*os++ = "Weight 5";
	*os++ = "Weight 6";
	obstream << endl; */

	ostream_iterator<double> od(obstream, "\t ");
	for (unsigned int i=0; i < obVector.size(); i++) {
		Observation ob = obVector.at(i);
		*od++ = ob.getType();
		*od++ = ob.getRadius();
		*od++ = ob.getAltitude();
		// Unused 3rd dimension
		*od++ = ob.getTheta();
		*od++ = ob.getOb();
		*od++ = ob.getInverseError();
		for (unsigned int var = 0; var < numVars; var++)
			*od++ = ob.getWeight(var);

		obstream << endl;	
	}
		
	// Load the observations into a vector
	obs = new real[obVector.size()*11];
	for (unsigned int m=0; m < obVector.size(); m++) {
		int n = m*11;
		Observation ob = obVector.at(m);
		obs[n] = ob.getOb();
		obs[n+1] = ob.getInverseError();
		for (unsigned int var = 0; var < numVars; var++) {
			obs[n+2+var] = ob.getWeight(var);
		}
		obs[n+2+numVars] = ob.getRadius();
		obs[n+3+numVars] = ob.getAltitude();
		obs[n+4+numVars] = ob.getType();
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
	ifstream obstream("./Observations.in");
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
	obs = new real[obVector.size()*11];
	for (unsigned int m=0; m < obVector.size(); m++) {
		int n = m*11;
		Observation ob = obVector.at(m);
		obs[n] = ob.getOb();
		obs[n+1] = ob.getInverseError();
		for (unsigned int var = 0; var < numVars; var++) {
			obs[n+2+var] = ob.getWeight(var);
		}
		obs[n+2+numVars] = ob.getRadius();
		obs[n+3+numVars] = ob.getAltitude();
		obs[n+4+numVars] = ob.getType();
	}	
	
	return true;
}

bool VarDriverRZ::loadBGfromTextfile()
{
	
	// Read in the background state from a formatted textfile
	// Read the r and v pairs from a file
	double height, radius, v, psi, h, q, rho;
	int zi = 0;
	vector<real> vIn, rvBG, psiBG, hBG, qBG, rpBG;
	vector<real>* vInit = new vector<real>[maxHeights];
	ifstream vdata("./RZbackground.in");
	vdata.width(14);
	while (vdata >> height >> radius >> v >> psi >> h >> q >> rho)
	{
		//cout << height << "\t" << radius << "\t" << v << "\t" << psi << "\t"  << h << "\t" << q << "\t" << rho << endl;
		if (z.empty()) z.push_back (height);
		if (height != z.back()) {
			// Assign the initial background fields
			BG[zi][0] =  rvBG;
			BG[zi][1] =  psiBG;
			BG[zi][2] =  hBG;
			BG[zi][3] =  qBG;
			BG[zi][4] =  rpBG;
			BG[zi][5] =  vIn;
			r.clear(); vIn.clear();
			rvBG.clear(); psiBG.clear();
			hBG.clear(); qBG.clear(); rpBG.clear();
			z.push_back (height);
			zi++;
		}		
		r.push_back (radius);
		
		// No negative v for potential radius transform!!!
		// Need to handle this more gracefully
		if (v < 0)  v = 0;
		real rhoBar = rhoBase*exp(-rhoInvScaleHeight*height);
		//real qBar = 19.562 - 0.004066*height + 7.8168e-7*height*height;
		real qBar = 18.189 - 0.0060609*height + 8.37e-7*height*height - 4.6009e-11*height*height*height;
		if (qBar < 0) qBar = 0;

		//real hBar = 3.5e5;
		real hBar = 3.4633e5 - 10.751*height + 0.0021949*height*height - 1.7019e-7*height*height*height
			+ 4.858e-12*height*height*height*height;
		vIn.push_back (v);
		rvBG.push_back (rho*radius*v*1.e-3);
		psiBG.push_back (psi*1.e-8);
		hBG.push_back ((h-hBar)*1.e-3);
		qBG.push_back (q-qBar);
		rpBG.push_back ((rho/(1 + q/1000.) - rhoBar)*100.);
	}
	
	// Assign the final strip
	BG[zi][0] =  rvBG;
	BG[zi][1] =  psiBG;
	BG[zi][2] =  hBG;
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
		exit(1);
	}
	
	return true;
	
}	

bool VarDriverRZ::loadBGfromFNL()
{
	
	/* Read in data from a GFS analysis
	GRIB_FILE gf;
	GRIB_MESSAGE m;
	const char* gribname = "fnl_080914_00_00.grb";
	int ret = gf.OpenRead(gribname);
	if (ret != 0) return -1;
	
	ret = gf.ReadMessage(m);
	if (ret != 0) return -1;
	
	int npix = m.grid.nx;
	int nlin = m.grid.ny;
	
	size_t total_size = m.grid.nxny;
	printf( "%s_%4d%02d%02d_%02d%02d", m.field.VarName( ),
			m.gtime.year, m.gtime.month, m.gtime.day,
			m.gtime.hour, m.gtime.minute);
	*/
	return false;
}

bool VarDriverRZ::bilinearMish()
{
	
	// Do a simple bilinear interpolation to the Mish
	// Num levels set to original grid first
	numHeights = z.size();
	
	///Realign the vertical grid on the Gaussian points
	double zMin = 0.;
	double zMax = z.back();
	int numZnodes = (int)((zMax - zMin)/zincr) + 1;
	double rMin = 0.;
	double rMax = (double)r.back();
	int numrnodes = (int)((rMax - rMin)/rincr) + 1;
	// Load the BG into a vector
	bgU = new real[4*(numrnodes-1)*(numZnodes-1)*(numVars-1)];
	//bgU[(numVars-1)*r.size()*zi +(numVars-1)*ri + vi] = BG[zi][vi].at(ri);

	// Set the master dimensions
	imin = rMin;
	imax = rMax;
	idim = numrnodes;
	jmin = zMin;
	jmax = zMax;
	jdim = numZnodes;
	
	for (int ri = 0; ri < (idim-1); ri++) {
		for (int rmu = -1; rmu <= 1; rmu += 2) {
			real rad = rMin + rincr * (ri + (0.5*sqrt(1./3.) * rmu + 0.5));
			for (int zi = 0; zi < (jdim-1); zi++) {
				for (int zmu = -1; zmu <= 1; zmu += 2) {
					real height = zMin + zincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
					int ii = -1;
					int jj = -1;
					// Find the brackets
					for (unsigned int i=0; i < r.size()-1; i++) {
						if ((r.at(i) <= rad) and (r.at(i+1) > rad)) {
							ii = i; 
							break;
						}
					}
					for (unsigned int j=0; j < z.size()-1; j++) {
						if ((z.at(j) <= height) and (z.at(j+1) > height)) {
							jj = j; 
							break;
						}
					}
					if ((ii < 0) or (jj < 0)) { cout << "Problem in bilinear interpolation!\n"; break; }
					real rmid = (rad-r.at(ii))/(r.at(ii+1)-r.at(ii));
					real zmid = (height-z.at(jj))/(z.at(jj+1)-z.at(jj));
					int bgZ = zi*2 + (zmu+1)/2;
					int bgR = ri*2 + (rmu+1)/2;
					for (unsigned int vi = 0; vi < (numVars-1); vi++) {
						bgU[(numVars-1)*(numrnodes-1)*2*bgZ +(numVars-1)*bgR + vi] =
						(1-rmid)*(1-zmid)*BG[jj][vi].at(ii) + rmid*(1-zmid)*BG[jj][vi].at(ii+1) +
						rmid*zmid*BG[jj+1][vi].at(ii+1) + (1-rmid)*zmid*BG[jj+1][vi].at(ii);
						//if (vi == 0) {
						//	cout << bgU[(numVars-1)*(numrnodes-1)*2*bgZ +(numVars-1)*bgR + vi] <<
						//	"\t" << BG[jj][vi].at(ii) << "\t" << BG[jj][vi].at(ii+1) <<
						//	"\t" << BG[jj+1][vi].at(ii+1) << "\t" << BG[jj+1][vi].at(ii) << endl;
						//}
					}
				}
			}
		}
	}
		
	return true;
	
}

real VarDriverRZ::bilinearField(real radius, real height, int var)
{
	
	int ii = -1;
	int jj = -1;
	// Find the brackets
	for (unsigned int i=0; i < r.size()-1; i++) {
		if ((r.at(i) <= radius) and (r.at(i+1) > radius)) {
			ii = i; 
			break;
		}
	}
	for (unsigned int j=0; j < z.size()-1; j++) {
		if ((z.at(j) <= height) and (z.at(j+1) > height)) {
			jj = j; 
			break;
		}
	}
	if ((ii < 0) or (jj < 0)) { 
		cout << "Problem in bilinear interpolation!\n"; 
		return -999.; 
	}
	real rmid = (radius-r.at(ii))/(r.at(ii+1)-r.at(ii));
	real zmid = (height-z.at(jj))/(z.at(jj+1)-z.at(jj));
	real field= 
		(1-rmid)*(1-zmid)*BG[jj][var].at(ii) + rmid*(1-zmid)*BG[jj][var].at(ii+1) +
		rmid*zmid*BG[jj+1][var].at(ii+1) + (1-rmid)*zmid*BG[jj+1][var].at(ii);
	
	return field;
}
	

bool VarDriverRZ::setupMishAndRXform()
{
	
	
	/* Set up Anisotropic Filter
	real* tau = new real[numrnodes*numZnodes];
	real* Rtmp = new real[z.size()];
	bgr = new real[numrnodes];
	
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
	//rtau[0] = tau[1]-tau[0];
	for (int rb = 0; rb < numrnodes; rb++) {
		rtau[rb] = tau[rb];
		cout << rb << "\t" << tau[rb] << "\t" << rtau[rb] << endl;
	}
	RecursiveFilter* anifilter = new RecursiveFilter(5, rtau, numrnodes);
	RecursiveFilter* isofilter = new RecursiveFilter(4, 5);
	isofilter->filterArray(bgr, numrnodes);
	cout << "Isotropic\n";
	for (int rb = 0; rb < numrnodes; rb++) {
		cout << rb << "\t" << bgr[rb] << "\n";
		bgr[rb] = bgB[(numVars-1)*rb];
	}
	cout << "Ansotropic\n";
	anifilter->aniFilterArray(bgr,numrnodes);
	for (int rb = 0; rb < numrnodes; rb++) {
		cout << rb << "\t" << bgr[rb] << "\n";
	} 
	delete[] Rtmp;
	delete[] tau;
	delete[] bgr;
	
	*/
	
	return true;
}

bool VarDriverRZ::initialize()
{
	// Run a RZ vortex background field
	cout << "Initializing RZ Background" << endl;

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
		loadBGfromTextfile();
	} else {
		loadBGfromFNL();
	}
	// Get the parametric winds from the input winds if data does not extend through domain
	//ParametricVortex parmVortex (&BG[zi][5].front(), &r[0],  r.size());
	//double vParm = parmVortex.getWindAtRadiusKM(grad);
		
	// Set up the data on the Gaussian Grid
	bilinearMish();
	
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
	int stateSize = 4*(idim-1)*(jdim-1)*(numVars-1);
	
	costRZ = new CostFunctionRZ_CPU(obVector.size(), stateSize);
	costRZ->initialize(imin, imax, idim, jmin, jmax, jdim, ia, ja, bgU, obs, RnumGridpts, RXform); 
	
	return true;
}

bool VarDriverRZ::run()
{

	double CQRMS = 999;
	int iter=0;
	while ((CQRMS > CQTOL) and (iter < maxIter)) {
		iter++;
		cout << "Outer Loop Iteration: " << iter << endl;
		costRZ->initState();
		costRZ->minimize();
		// Increment the variables
		costRZ->updateBG();
	}	
	cout << "Increment RMS Tolerance of " << CQTOL << " reached in "
		<< iter << " iterations. Writing analysis results..." << endl;

	return true;

}


