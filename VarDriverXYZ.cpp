/*
 *  VarDriverXYZ.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriverXYZ.h"
#include "Dorade.h"
#include <iterator>
#include <fstream>
#include <cmath>
#include <QTextStream>
#include <QFile>
#include <QVector>
#include <iomanip>
#include "RecursiveFilter.h"
#include <GeographicLib/TransverseMercatorExact.hpp>

// Constructor
VarDriverXYZ::VarDriverXYZ()
	: VarDriver()
{
	numVars = 7;
}

// Destructor
VarDriverXYZ::~VarDriverXYZ()
{
}

/* This routine is the main initializer of the analysis */

bool VarDriverXYZ::initialize(const QDomElement& configuration)
{
	// Run a XYZ vortex background field
	cout << "Initializing SAMURAI XYZ" << endl;
	
	// Parse the XML configuration file
	if (!parseXMLconfig(configuration)) return false;
	
	// Validate the XYZ specific parameters
	if (!validateXMLconfig()) return false;
	
	// Define the grid dimensions
	imin = configHash.value("xmin").toFloat();
	imax = configHash.value("xmax").toFloat();
	iincr = configHash.value("xincr").toFloat();
	idim = (int)((imax - imin)/iincr) + 1;
	
	jmin = configHash.value("ymin").toFloat();
	jmax = configHash.value("ymax").toFloat();
	jincr = configHash.value("yincr").toFloat();
	jdim = (int)((jmax - jmin)/jincr) + 1;
	
	kmin = configHash.value("zmin").toFloat();
	kmax = configHash.value("zmax").toFloat();
	kincr = configHash.value("zincr").toFloat();
	kdim = (int)((kmax - kmin)/kincr) + 1;

	// The recursive filter uses a fourth order stencil to spread the observations, so less than 4 gridpoints will cause a memory fault
	if (idim < 4) {
		cout << "X dimension is less than 4 gridpoints and recursive filter will fail. Aborting...\n";
		return false;
	}
	if (jdim < 4) {
		cout << "Y dimension is less than 4 gridpoints and recursive filter will fail. Aborting...\n";
		return false;
	}	
	if (kdim < 4) {
		cout << "Z dimension is less than 4 gridpoints and recursive filter will fail. Aborting...\n";
		return false;
	}
	
	// Define the sizes of the arrays we are passing to the cost function
	cout << "xMin\txMax\txIncr\tyMin\tyMax\tyIncr\tzMin\tzMax\tzIncr\n";
	cout << imin << "\t" <<  imax << "\t" <<  iincr << "\t";
	cout << jmin << "\t" <<  jmax << "\t" <<  jincr << "\t";
	cout << kmin << "\t" <<  kmax << "\t" <<  kincr << "\n\n";

	// Increase the "internal" size of the grid for the zero BC condition
	if (configHash.value("horizontalbc") == "R0") {
		imin -= iincr;
		imax += iincr;
		idim	+= 2;
		jmin -= jincr;
		jmax += jincr;
		jdim += 2;
	}
	
	if (configHash.value("verticalbc") == "R0") {
		kmin -= kincr;
		kmax += kincr;
		kdim += 2;
	}

	int uStateSize = 8*(idim-1)*(jdim-1)*(kdim-1)*(numVars);
	int bStateSize = idim*jdim*kdim*numVars;
	cout << "Physical (mish) State size = " << uStateSize << "\n";
	cout << "Nodal State size = " << bStateSize << ", Grid dimensions:\n";
	
	// Load the BG into a empty vector
	bgU = new real[uStateSize];
	bgWeights = new real[uStateSize];
	for (int i=0; i < uStateSize; i++) {
		bgU[i] = 0.;
		bgWeights[i] = 0.;
	}		
	
	// Define the Reference state
	if (configHash.value("refstate") == "dunion_mt") {
		referenceState = dunion_mt;
	} else {
		cout << "Reference state not defined!\n";
		exit(-1);
	}
	
	cout << "Reference profile: Z\t\tQv\tRhoa\tRho\tH\tTemp\tPressure\n";
	for (real k = kmin; k < kmax+kincr; k+= kincr) {
		cout << "                   " << k << "\t";
		for (int i = 0; i < 6; i++) {
			real var = getReferenceVariable(i, k*1000);
			if (i == 0) var = bhypInvTransform(var);
			cout << setw(9) << setprecision(4)  << var << "\t";
		}
		cout << "\n";
	}
	cout << setprecision(9);
	
	// Read in the Frame centers
	// Ideally, create a time-based spline from limited center fixes here
	// but just load 1 second centers into vector for now
	readFrameCenters();
	
	// Get the reference center
	QTime reftime = QTime::fromString(configHash.value("reftime"), "hh:mm:ss");
	QString refstring = reftime.toString();
	bool foundref = false;
	for (unsigned int fi = 0; fi < frameVector.size(); fi++) {
		QDateTime frametime = frameVector[fi].getTime();
		if (reftime == frametime.time()) {
			QString tempstring;
			QDate refdate = frametime.date();
			QDateTime unixtime(refdate, reftime, Qt::UTC);
			configHash.insert("reflat", tempstring.setNum(frameVector[fi].getLat()));
			configHash.insert("reflon", tempstring.setNum(frameVector[fi].getLon()));
			configHash.insert("reftime", tempstring.setNum(unixtime.toTime_t()));
			cout << "Found matching reference time " << refstring.toStdString()
			<< " at " << frameVector[fi].getLat() << ", " << frameVector[fi].getLon() << "\n";
			foundref = true;
			break;
		}
	}
	if (!foundref) {
		cout << "Error finding reference time, please check date and time in XML file\n";
		return false;
	}
	
	/* Set the maximum number of iterations to the multipass reduction factor
	Multiple outer loops will reduce the cutoff wavelengths and background error variance */
	maxIter = configHash.value("num_iterations").toFloat();

	/* Optionally load a set of background estimates and interpolate to the Gaussian mish */
	bool loadBG = configHash.value("load_background").toInt();
	int numbgObs = 0;
	if (loadBG) {
		numbgObs = loadBackgroundObs();
		if (numbgObs < 0) {
			cout << "Error loading background file\n";
			return false;
		}
	}
	
	/* Optionally adjust the interpolated background to satisfy mass continuity
	 and match the supplied points exactly. In essence, do a SAMURAI analysis using
	 the background estimates as "observations" */
	bool adjustBG = configHash.value("adjust_background").toInt();
	if (adjustBG and numbgObs) {
		if (!adjustBackground(bStateSize)) {
			cout << "Error adjusting background\n";
			return false;
		}
	}
	
	// Read in the observations, process them into weights and positions
	// Either preprocess from raw observations or load an already processed Observations.in file
	bool preprocess = configHash.value("preprocess_obs").toInt();
	if (preprocess) {
		if (!preProcessMetObs()) {
			cout << "Error pre-processing observations\n";
			return false;
		}
	} else {
		if (!loadMetObs()) {
			cout << "Error loading observations\n";
			return false;
		}
	}
	cout << "Number of New Observations: " << obVector.size() << endl;		
	
	// We are done with the bgWeights, so free up that memory
	delete[] bgWeights;	
	
	obCostXYZ = new CostFunctionXYZ(obVector.size(), bStateSize);
	obCostXYZ->initialize(&configHash, bgU, obs);
	
	// If we got here, then everything probably went OK!
	return true;
}

/* This routine drives the CostFunction minimization
 There is support for an outer loop to change the background
 error covariance or update non-linear observation operators */

bool VarDriverXYZ::run()
{
	int iter=1;
	while (iter <= maxIter) {
		cout << "Outer Loop Iteration: " << iter << endl;
		obCostXYZ->initState(iter);
		obCostXYZ->minimize();
		obCostXYZ->updateBG();
		iter++;
		
		// Optionally update the analysis parameters for an additional iteration
		updateAnalysisParams(iter);
	}	
	
	return true;
	
}

/* Clean up all that allocated memory */

bool VarDriverXYZ::finalize()
{
	obCostXYZ->finalize();
	delete[] obs;
	delete[] bgU;
	delete obCostXYZ;	
	return true;
}

/* Pre-process the observations into a single vector
 On the wishlist is some integrated QC here other than just spatial thresholding */

bool VarDriverXYZ::preProcessMetObs()
{
	
	vector<real> rhoP;

	// Geographic functions
	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM;
	real referenceLon = configHash.value("reflon").toFloat();
	
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
					cout << "Error reading cls file" << endl;
				break;
			case (sec):
				if (!read_sec(metFile, metData))
					cout << "Error reading sec file" << endl;
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
			case (eol):
				if (!read_eol(metFile, metData))
					cout << "Error reading eol file" << endl;
				break;
			case (qscat):
				if (!read_qscat(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (ascat):
				if (!read_ascat(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (nopp):
				if (!read_nopp(metFile, metData))
					cout << "Error reading wwind file" << endl;
				break;
			case (cimss):
				if (!read_cimss(metFile, metData))
					cout << "Error reading cimss file" << endl;
				break;
			case (dwl):
				if (!read_dwl(metFile, metData))
					cout << "Error reading dwl file" << endl;
				break;
			case (cen):
				continue;				
			default:
				cout << "Unknown data type, skipping..." << endl;
				continue;
		}
		
		processedFiles++;
		
		// Process the metObs into Observations
		QDateTime startTime = frameVector.front().getTime();
		QDateTime endTime = frameVector.back().getTime();
		for (int i = 0; i < metData->size(); ++i) {
			
			// Make sure the ob is within the time limits
			MetObs metOb = metData->at(i);
			QDateTime obTime = metOb.getTime();
			QString obstring = obTime.toString(Qt::ISODate);
			QString tcstart = startTime.toString(Qt::ISODate);
			QString tcend = endTime.toString(Qt::ISODate);		
			if ((obTime < startTime) or (obTime > endTime)) continue;
			int fi = startTime.secsTo(obTime);
			if ((fi < 0) or (fi > (int)frameVector.size())) {
				cout << "Time problem with observation " << fi << endl;
				continue;
			}
			real Um = frameVector[fi].getUmean();
			real Vm = frameVector[fi].getVmean();
						
			// Get the X, Y & Z
			real tcX, tcY, metX, metY;
			tm.Forward(referenceLon, frameVector[fi].getLat() , frameVector[fi].getLon() , tcX, tcY);
			tm.Forward(referenceLon, metOb.getLat() , metOb.getLon() , metX, metY);
			real obX = (metX - tcX)/1000.;
			real obY = (metY - tcY)/1000.;
			real heightm = metOb.getAltitude();
			real obZ = heightm/1000.;
			
			// Make sure the ob is in the domain
			if ((obX < imin) or (obX > imax) or
				(obY < jmin) or (obY > jmax) or
				(obZ < kmin) or (obZ > kmax))
				continue;
			
			// Restrict the horizontal domain if we are using the R0 BC
			if (configHash.value("horizontalbc") == "R0") {
				if ((obX < (imin+iincr)) or (obX > (imax-iincr)) or
					(obY < (jmin+jincr)) or (obY > (jmax-jincr))) 
					continue;
			}

			if (configHash.value("verticalbc") == "R0") {
				if ((obZ < (kmin+kincr)) or (obZ > (kmax-kincr)))
					continue;
			}
			
			
			// Create an observation and set its basic info
			Observation varOb;
			varOb.setCartesianX(obX);
			varOb.setCartesianY(obY);
			varOb.setAltitude(obZ);
			varOb.setTime(obTime.toTime_t());
			
			// Reference states			
			real rhoBar = getReferenceVariable(rhoaref, heightm);
			real qBar = getReferenceVariable(qvbhypref, heightm);
			real tBar = getReferenceVariable(tempref, heightm);
			
			// Initialize the weights
			varOb.setWeight(0., 0);
			varOb.setWeight(0., 1);
			varOb.setWeight(0., 2);
			varOb.setWeight(0., 3);
			varOb.setWeight(0., 4);
			varOb.setWeight(0., 5);
			varOb.setWeight(0., 6);
			real u, v, w, rho, rhoa, qv, tempk, rhov, rhou, rhow, wspd; 
			switch (metOb.getObType()) {
				case (MetObs::dropsonde):
					varOb.setType(MetObs::dropsonde);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					w = metOb.getVerticalVelocity();
					rho = metOb.getMoistDensity();
					rhoa = metOb.getAirDensity();
					qv = metOb.getQv();
					tempk = metOb.getTemperature();
					
					// Separate obs for each measurement
					// rho v 1 m/s error
					if ((u != -999) and (rho != -999)) {
						// rho u 1 m/s error
						varOb.setWeight(1., 0);
						rhou = rho*(u - Um);
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						rhov = rho*(v - Vm);
						varOb.setOb(rhov);
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
					if (tempk != -999) {
						// temperature 1 K error
						varOb.setWeight(1., 3);
						varOb.setOb(tempk - tBar);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 3);
					}
					if (qv != -999) {
						// Qv 0.5 g/kg error
						varOb.setWeight(1., 4);
						qv = bhypTransform(qv);
						varOb.setOb(qv-qBar);
						varOb.setError(0.5);
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
					tempk = metOb.getTemperature();
					
					// Separate obs for each measurement
					// rho v 1 m/s error
					if ((u != -999) and (rho != -999)) {						
						// rho u 1 m/s error
						varOb.setWeight(1., 0);
						rhou = rho*(u - Um);
						varOb.setOb(rhou);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						rhov = rho*(v - Vm);
						varOb.setOb(rhov);
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
					if (tempk != -999) {
						// temperature 1 K error
						varOb.setWeight(1., 3);
						varOb.setOb(tempk - tBar);
						varOb.setError(1.0);
						obVector.push_back(varOb);
						varOb.setWeight(0., 3);
					}
					if (qv != -999) {
						// Qv 0.5 g/kg error
						varOb.setWeight(1., 4);
						qv = bhypTransform(qv);
						varOb.setOb(qv-qBar);
						varOb.setError(0.5);
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
					// This needs to be redone for the Cartesian case
					//vBG = 1.e3*bilinearField(obX, obY, 0);
					//uBG = -1.e5*bilinearField(obX, 20., 1)/(rad*20.);
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
						// rho u 1 m/s error
						varOb.setWeight(1., 0);
						rhou = (u - Um);
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						// Multiply by rho later from grid values
						rhov = (v - Vm);
						varOb.setOb(rhov);
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
						// rho u 1 m/s error
						varOb.setWeight(1., 0);
						rhou = (u - Um);
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						// Multiply by rho later from grid values
						rhov = (v - Vm);
						varOb.setOb(rhov);
						varOb.setError(2.5);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);						
					}
					break;
					
				case (MetObs::AMV):
					varOb.setType(MetObs::AMV);
					u = metOb.getCartesianUwind();
					v = metOb.getCartesianVwind();
					if (u != -999) {
						// rho u 10 m/s error
						varOb.setWeight(1., 0);
						rhou = (u - Um);
						//cout << "RhoU: " << rhou << endl;
						varOb.setOb(rhou);
						varOb.setError(10.);
						obVector.push_back(varOb);
						varOb.setWeight(0., 0);
						
						varOb.setWeight(1., 1);
						// Multiply by rho later from grid values
						rhov = (v - Vm);
						varOb.setOb(rhov);
						varOb.setError(10.);
						obVector.push_back(varOb);
						varOb.setWeight(0., 1);						
					}
					break;
					
				case (MetObs::lidar):
				{
					varOb.setType(MetObs::lidar);
					// Geometry terms
					real az = metOb.getAzimuth()*Pi/180.;
					real el = metOb.getElevation()*Pi/180.;
					real uWgt = sin(az)*cos(el);
					real vWgt = cos(az)*cos(el);
					real wWgt = sin(el);
					
					// Fall speed is assumed zero since we are dealing with aerosols
					real db = metOb.getReflectivity();
					real vr = metOb.getRadialVelocity();
					real w_term = 0.0;  
					real Vdopp = vr - w_term*sin(el) - Um*sin(az)*cos(el) - Vm*cos(az)*cos(el);
					
					varOb.setWeight(uWgt, 0);
					varOb.setWeight(vWgt, 1);
					varOb.setWeight(wWgt, 2);
										
					// Set the error according to the spectrum width and power
					real DopplerError = metOb.getSpectrumWidth() + log(50/db);
					if (DopplerError < 1.0) DopplerError = 1.0;
					varOb.setError(DopplerError);
					varOb.setOb(Vdopp);
					obVector.push_back(varOb);
					varOb.setWeight(0., 0);	
					varOb.setWeight(0., 1);	
					varOb.setWeight(0., 2);
					
					break;
				}	
				case (MetObs::radar):
				{
					varOb.setType(MetObs::radar);
					// Geometry terms
					real az = metOb.getAzimuth()*Pi/180.;
					real el = metOb.getElevation()*Pi/180.;
					real uWgt = sin(az)*cos(el);
					real vWgt = cos(az)*cos(el);
					real wWgt = sin(el);
					
					// Fall speed
					real Z = metOb.getReflectivity();
					real H = metOb.getAltitude();
					real ZZ=pow(10.0,(Z*0.1));
					real zeroC = 4800.;
					real hlow= zeroC; 
					real hhi= hlow + 1000;
					
					/* density correction term (rhoo/rho)*0.45 
					 0.45 density correction from Beard (1985, JOAT pp 468-471) */
					real rho = getReferenceVariable(rhoref, H);
					real rhosfc = getReferenceVariable(rhoref, 0.);
					real DCOR = pow((rhosfc/rho),(real)0.45);
					
					// The snow relationship (Atlas et al., 1973) --- VT=0.817*Z**0.063  (m/s) 
					real VTS=-DCOR * (0.817*pow(ZZ,(real)0.063));
					
					// The rain relationship (Joss and Waldvogel,1971) --- VT=2.6*Z**.107 (m/s) */
					real VTR=-DCOR * (2.6*pow(ZZ,(real).107));
					
					/* Test if height is in the transition region between SNOW and RAIN
					   defined as hlow in km < H < hhi in km
					   if in the transition region do a linear weight of VTR and VTS */
					if ((Z > 20) and (Z <= 30)) {
						real WEIGHTR=(Z-20)/(10);
						real WEIGHTS=1.-WEIGHTR;
						VTS=(VTR*WEIGHTR+VTS*WEIGHTS)/(WEIGHTR+WEIGHTS);
					} else if (Z > 30) {
						VTS=VTR;
					}
					real w_term=VTR*(hhi-H)/1000 + VTS*(H-hlow)/1000;  
					if (H < hlow) w_term=VTR; 
					if (H > hhi) w_term=VTS;
					real Vdopp = metOb.getRadialVelocity() - w_term*sin(el) - Um*sin(az)*cos(el) - Vm*cos(az)*cos(el);
					
					varOb.setWeight(uWgt, 0);
					varOb.setWeight(vWgt, 1);
					varOb.setWeight(wWgt, 2);
					
					/* Theoretically, rhoPrime could be included as a prognostic variable here...
					   However, adding another unknown without an extra equation makes the problem even more underdetermined
					   so assume it is small and ignore it
					   real rhopWgt = -Vdopp;
					  varOb.setWeight(rhopWgt, 5); */
					
					// Set the error according to the spectrum width and potential fall speed error (assume 2 m/s?)
					real DopplerError = metOb.getSpectrumWidth() + fabs(wWgt)*2.;
					if (DopplerError < 1.0) DopplerError = 1.0;
					varOb.setError(DopplerError);
					varOb.setOb(Vdopp);
					obVector.push_back(varOb);
					varOb.setWeight(0., 0);	
					varOb.setWeight(0., 1);	
					varOb.setWeight(0., 2);
					
					// Reflectivity observations
					QString gridref = configHash.value("qrvariable");
					real qr = 0.;
					if (gridref == "qr") {
						// Do the gridding as part of the variational synthesis using Z-M relationships
						// Z-M relationships from Gamache et al (1993) JAS
						real rainmass = pow(ZZ/14630.,(real)0.6905);
						real icemass = pow(ZZ/670.,(real)0.5587);
						if ((Z > 20) and (Z <= 30)) {
							real WEIGHTR=(Z-20)/(10);
							real WEIGHTS=1.-WEIGHTR;
							icemass=(rainmass*WEIGHTR+icemass*WEIGHTS)/(WEIGHTR+WEIGHTS);
						} else if (Z > 30) {
							icemass=rainmass;
						}
						
						real precipmass = rainmass*(hhi-H)/1000 + icemass*(H-hlow)/1000;
						if (H < hlow) precipmass = rainmass;
						if (H > hhi) precipmass = icemass;
						qr = bhypTransform(precipmass/rhoBar);
						
						/* Include an observation of this quantity in the variational synthesis
						varOb.setOb(qr);
						varOb.setWeight(1., 6);
						varOb.setError(1.0);
						obVector.push_back(varOb); */
						
					} else if (gridref == "dbz") {
						qr = (Z+35.)*0.1;
						/* Include an observation of this quantity in the variational synthesis
						 varOb.setOb(qr);
						 varOb.setWeight(1., 6);
						 varOb.setError(1.0);
						 obVector.push_back(varOb); */
						
					}

					// Do a Exponential & power weighted interpolation of the reflectivity/qr in a grid box
					real ROI = configHash.value("reflectivityroi").toFloat();
					real Rsquare = (iincr*ROI)*(iincr*ROI) + (jincr*ROI)*(jincr*ROI) + (kincr*ROI)*(kincr*ROI);
					for (int zi = 0; zi < (kdim-1); zi++) {	
						for (int zmu = -1; zmu <= 1; zmu += 2) {
							real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
							if (fabs(zPos-obZ) > kincr*ROI*2.) continue;
							for (int xi = 0; xi < (idim-1); xi++) {
								for (int xmu = -1; xmu <= 1; xmu += 2) {
									real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
									if (fabs(xPos-obX) > iincr*ROI*2.) continue;
									
									for (int yi = 0; yi < (jdim-1); yi++) {
										for (int ymu = -1; ymu <= 1; ymu += 2) {
											real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
											if (fabs(yPos-obY) > jincr*ROI*2.) continue;
											real rSquare = (obX-xPos)*(obX-xPos) + (obY-yPos)*(obY-yPos) + (obZ-zPos)*(obZ-zPos); 
											int bgI = xi*2 + (xmu+1)/2;
											int bgJ = yi*2 + (ymu+1)/2;
											int bgK = zi*2 + (zmu+1)/2;
											int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
											if (rSquare < Rsquare) {
												real weight = exp(-2.302585092994045*rSquare/Rsquare);
												//real weight = (Rsquare - rSquare)/(Rsquare + rSquare);
												//if (qr > bgU[bIndex +6]) bgU[bIndex +6] = qr;
												bgU[bIndex +6] += weight*qr;
												bgWeights[bIndex] += weight;
											}
										}
									}
								}
							}
						}
					}
					
					break;
				}	
							
			}

		} 
		cout << obVector.size() << " total observations." << endl;
	}
	
	delete metData;
	
	// Finish reflectivity interpolation
	Observation varOb;
	varOb.setTime(configHash.value("reftime").toFloat());	
	real pseudow_weight = configHash.value("use_dbz_pseudow").toFloat();
	varOb.setWeight(0., 0);
	varOb.setWeight(0., 1);
	varOb.setWeight(1., 2);
	varOb.setWeight(0., 3);
	varOb.setWeight(0., 4);
	varOb.setWeight(0., 5);
	varOb.setWeight(0., 6);	
	varOb.setError(pseudow_weight);
	varOb.setOb(0.);
	for (int xi = 0; xi < (idim-1); xi++) {
		for (int xmu = -1; xmu <= 1; xmu += 2) {
			real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
			
			for (int yi = 0; yi < (jdim-1); yi++) {
				for (int ymu = -1; ymu <= 1; ymu += 2) {
					real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
					real maxrefHeight = -1;
					for (int zi = 0; zi < (kdim-1); zi++) {	
						for (int zmu = -1; zmu <= 1; zmu += 2) {
							real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
							
							int bgI = xi*2 + (xmu+1)/2;
							int bgJ = yi*2 + (ymu+1)/2;
							int bgK = zi*2 + (zmu+1)/2;
							int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
							if (bgWeights[bIndex] != 0) {
								bgU[bIndex +6] /= bgWeights[bIndex];
							}
							if (bgU[bIndex +6] > 0) {
								maxrefHeight = zPos;
							}
						}
					}

					// Set an upper boundary condition for W
					if ((maxrefHeight > 0) and (maxrefHeight < kmax)
						and (pseudow_weight > 0.0)) {
						varOb.setCartesianX(xPos);
						varOb.setCartesianY(yPos);
						varOb.setAltitude(maxrefHeight);					
						obVector.push_back(varOb);
					}
					
					// Set a lower boundary condition for W
					// Ideally use a terrain map here, but just use Z=0 for now
					if (pseudow_weight > 0.0) {
						varOb.setCartesianX(xPos);
						varOb.setCartesianY(yPos);
						varOb.setAltitude(0.0);					
						obVector.push_back(varOb);
					}
					
					
				}
			}
		}
	}	
	cout << obVector.size() << " total observations including pseudo W obs" << endl;
	
	// Write the Obs to a summary text file
	ofstream obstream("samurai_Observations.in");
	// Header messes up reload
	/*ostream_iterator<string> os(obstream, "\t ");
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
	obstream << endl; */

	ostream_iterator<real> od(obstream, "\t ");
	ostream_iterator<int> oi(obstream, "\t ");
	for (int i=0; i < obVector.size(); i++) {
		Observation ob = obVector.at(i);
		*od++ = ob.getOb();
		*od++ = ob.getInverseError();
		*od++ = ob.getCartesianX();
		*od++ = ob.getCartesianY();
		*od++ = ob.getAltitude();		
		*oi++ = ob.getType();
		*oi++ = ob.getTime();
		for (unsigned int var = 0; var < numVars; var++)
			*od++ = ob.getWeight(var);

		obstream << endl;	
	}
	
	// Load the observations into a vector
	obs = new real[obVector.size()*14];
	for (int m=0; m < obVector.size(); m++) {
		int n = m*14;
		Observation ob = obVector.at(m);
		obs[n] = ob.getOb();
		obs[n+1] = ob.getInverseError();
		obs[n+2] = ob.getCartesianX();
		obs[n+3] = ob.getCartesianY();
		obs[n+4] = ob.getAltitude();
		obs[n+5] = ob.getType();
		obs[n+6] = ob.getTime();
		for (unsigned int var = 0; var < numVars; var++) {
			obs[n+7+var] = ob.getWeight(var);
		}
	}	
	
	// All done preprocessing
	if (!processedFiles) {
		cout << "No files processed, nothing to do :(" << endl;
		// return 0;
	} else {
		cout << "Finished preprocessing " << processedFiles << " files." << endl;
	}
	
	return true;
}

/* Load the meteorological observations from a file into a vector */

bool VarDriverXYZ::loadMetObs()
{
	
	Observation varOb;
	real wgt[numVars];
	real xPos, yPos, zPos, ob, error;
	int type;
	int time;
	cout << "Loading preprocessed observations from samurai_Observations.in" << endl;
	
	// Open and read the file
	ifstream obstream("samurai_Observations.in");
	while (obstream >> ob >> error >> xPos >> yPos >> zPos >> type >> time
		   >> wgt[0] >> wgt[1] >> wgt[2] >> wgt[3] >> wgt[4] >> wgt[5] >> wgt[6])
	{
		varOb.setOb(ob);
		varOb.setCartesianX(xPos);
		varOb.setCartesianY(yPos);
		varOb.setAltitude(zPos);
		varOb.setType(type);
		varOb.setTime(time);
		varOb.setError(1./error);
		for (unsigned int var = 0; var < numVars; var++)
			varOb.setWeight(wgt[var],var);
		obVector.push_back(varOb);
	}
	
	// Load the observations into the vector
	obs = new real[obVector.size()*14];
	for (int m=0; m < obVector.size(); m++) {
		int n = m*14;
		Observation ob = obVector.at(m);
		obs[n] = ob.getOb();
		obs[n+1] = ob.getInverseError();
		obs[n+2] = ob.getCartesianX();
		obs[n+3] = ob.getCartesianY();
		obs[n+4] = ob.getAltitude();
		obs[n+5] = ob.getType();
		obs[n+6] = ob.getTime();
		for (unsigned int var = 0; var < numVars; var++) {
			obs[n+7+var] = ob.getWeight(var);
		}
		
	}	
	
	return true;
}

/* Load the background estimates from a file */

int VarDriverXYZ::loadBackgroundObs()
{
	// Turn Debug on if there are problems with the vertical spline interpolation,
	// Eventually this should be replaced with the internal spline code
	// SplineD::Debug(1);

	// Geographic functions
	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM;
	real referenceLon = configHash.value("reflon").toFloat();
	
	QVector<real> logheights, uBG, vBG, wBG, tBG, qBG, rBG;
	SplineD* bgSpline;
	int time;
	QString bgTimestring, tcstart, tcend;
	real lat, lon, alt, u, v, w, t, qv, rhoa;
	real bgX, bgY, bgZ;
	// backgroundroi is in km, ROI is gridpoints
	real ROI = configHash.value("backgroundroi").toFloat() / iincr;
	real Rsquare = (iincr*ROI)*(iincr*ROI) + (jincr*ROI)*(jincr*ROI);
	ifstream bgstream("./samurai_Background.in");
	if (!bgstream.good()) {
		cout << "Error opening samurai_Background.in for reading.\n";
		exit(1);
	}
	cout << "Loading background onto Gaussian mish with " << ROI << " grid length radius of influence" << endl;
	
	while (bgstream >> time >> lat >> lon >> alt >> u >> v >> w >> t >> qv >> rhoa)
	{
	
		// Process the metObs into Observations
		QDateTime startTime = frameVector.front().getTime();
		QDateTime endTime = frameVector.back().getTime();
		
		// Make sure the bg is within the time limits
		QDateTime bgTime;
		bgTime.setTimeSpec(Qt::UTC);
		bgTime.setTime_t(time);
		bgTimestring = bgTime.toString(Qt::ISODate);
		tcstart = startTime.toString(Qt::ISODate);
		tcend = endTime.toString(Qt::ISODate);		
		if ((bgTime < startTime) or (bgTime > endTime)) continue;
		int tci = startTime.secsTo(bgTime);
		if ((tci < 0) or (tci > (int)frameVector.size())) {
			cout << "Time problem with observation " << tci << "secs more than center entries" << endl;
			continue;
		}
		
		real Um = frameVector[tci].getUmean();
		real Vm = frameVector[tci].getVmean();

		// Get the X, Y & Z
		real tcX, tcY, metX, metY;
		tm.Forward(referenceLon, frameVector[tci].getLat() , frameVector[tci].getLon() , tcX, tcY);
		tm.Forward(referenceLon, lat, lon , metX, metY);
		bgX = (metX - tcX)/1000.;
		bgY = (metY - tcY)/1000.;
		real heightm = alt;
		bgZ = heightm/1000.;
				
		// Make sure the ob is in the Interpolation domain
		if ((bgX < (imin-(ROI*iincr*2))) or (bgX > (imax+(ROI*iincr*2))) or
			(bgY < (jmin-(ROI*jincr*2))) or (bgY > (jmax+(ROI*jincr*2)))
			or (bgZ < kmin)) //Allow for higher values for interpolation purposes
			continue;

		// Reference states			
		real rhoBar = getReferenceVariable(rhoaref, heightm);
		real qBar = getReferenceVariable(qvbhypref, heightm);
		real tBar = getReferenceVariable(tempref, heightm);

		real rho = rhoa + rhoa*qv/1000.;
		real rhou = rho*(u - Um);
		real rhov = rho*(v - Vm);
		real rhow = rho*w;
		real tprime = t - tBar;
		qv = bhypTransform(qv);
		real qvprime = qv-qBar;
		real rhoprime = (rhoa-rhoBar)*100;
		real logZ = log(bgZ);
		// We assume here that the background precipitation field is always zero
		real qr = 0.;
		bgIn << bgX << bgY << logZ << time << rhou << rhov << rhow << tprime << qvprime << rhoprime << qr ;
		if (logheights.size() == 0) {
			// First column
			logheights.push_back(logZ);
			uBG.push_back(rhou);
			vBG.push_back(rhov);
			wBG.push_back(rhow);
			tBG.push_back(tprime);
			qBG.push_back(qvprime);
			rBG.push_back(rhoprime);
		} else if (logZ > logheights.back()) {
			// Same column
			logheights.push_back(logZ);
			uBG.push_back(rhou);
			vBG.push_back(rhov);
			wBG.push_back(rhow);
			tBG.push_back(tprime);
			qBG.push_back(qvprime);
			rBG.push_back(rhoprime);
		} else {
			// Solve for the spline
			if (logheights.size() == 1) {
				cerr << "Only one level found in background spline setup. Please check Background.in to ensure sorting by Z coordinate and re-run." << endl;
				return -1;
			}
			bgSpline = new SplineD(&logheights.front(), logheights.size(), uBG.data(), 0, SplineBase::BC_ZERO_SECOND);
			if (!bgSpline->ok()) {
				cerr << "bgSpline setup failed." << endl;
				return -1;
			}
			
			// Exponential interpolation in horizontal, b-Spline interpolation on log height in vertical
			for (int zi = 0; zi < (kdim-1); zi++) {	
				for (int zmu = -1; zmu <= 1; zmu += 2) {
					real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
					real logzPos = log(zPos);
					
					for (int xi = 0; xi < (idim-1); xi++) {
						for (int xmu = -1; xmu <= 1; xmu += 2) {
							real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
							if (fabs(xPos-bgX) > iincr*ROI*2.) continue;
							
							for (int yi = 0; yi < (jdim-1); yi++) {
								for (int ymu = -1; ymu <= 1; ymu += 2) {
									real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
									if (fabs(yPos-bgY) > jincr*ROI*2.) continue;
									
									real rSquare = (bgX-xPos)*(bgX-xPos) + (bgY-yPos)*(bgY-yPos);
									int bgI = xi*2 + (xmu+1)/2;
									int bgJ = yi*2 + (ymu+1)/2;
									int bgK = zi*2 + (zmu+1)/2;
									int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
									if (rSquare < Rsquare) {
										real weight = exp(-2.302585092994045*rSquare/Rsquare);
										if (logzPos > logheights.front()) {
											bgSpline->solve(uBG.data());
											bgU[bIndex] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(vBG.data());
											bgU[bIndex +1] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(wBG.data());
											bgU[bIndex +2] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(tBG.data());
											bgU[bIndex +3] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(qBG.data());
											bgU[bIndex +4] += weight*(bgSpline->evaluate(logzPos));
											bgSpline->solve(rBG.data());
											bgU[bIndex +5] += weight*(bgSpline->evaluate(logzPos));
											bgWeights[bIndex] += weight;
										} else {
											// Below the spline interpolation
											bgU[bIndex] += weight*uBG.front();
											bgU[bIndex +1] += weight*vBG.front();
											bgU[bIndex +2] += weight*wBG.front();
											bgU[bIndex +3] += weight*tBG.front();
											bgU[bIndex +4] += weight*qBG.front();
											bgU[bIndex +5] += weight*rBG.front();
											bgWeights[bIndex] += weight;
										}											
									}
								}
							}
						}
					}
				}
			}
			
			delete bgSpline;
			logheights.clear();
			uBG.clear();
			vBG.clear();
			wBG.clear();
			tBG.clear();
			qBG.clear();
			rBG.clear();
			
			logheights.push_back(log(bgZ));
			uBG.push_back(rhou);
			vBG.push_back(rhov);
			wBG.push_back(rhow);
			tBG.push_back(tprime);
			qBG.push_back(qvprime);
			rBG.push_back(rhoprime);			
		}
	}				

	if (!logheights.size()) {
		// Error reading in the background field
		cout << "No background estimates read in. Please check the time and location of your background field.\n";
		cout << "Observation window: " << tcstart.toStdString() << " to " << tcend.toStdString() << "\n";
		cout << "Background time: " << bgTimestring.toStdString() << "\n";
		return -1;
	}
	
	// Solve for the last spline
	bgSpline = new SplineD(&logheights.front(), logheights.size(), &uBG[0], 2, SplineBase::BC_ZERO_SECOND);
	if (!bgSpline->ok())
	{
		cerr << "bgSpline setup failed." << endl;
		return -1;
	}
	
	// Exponential interpolation in horizontal, b-Spline interpolation on log height in vertical
	for (int zi = 0; zi < (kdim-1); zi++) {	
		for (int zmu = -1; zmu <= 1; zmu += 2) {
			real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
			real logzPos = log(zPos);
			
			for (int xi = 0; xi < (idim-1); xi++) {
				for (int xmu = -1; xmu <= 1; xmu += 2) {
					real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
					if (fabs(xPos-bgX) > ROI*iincr*2.) continue;
					
					for (int yi = 0; yi < (jdim-1); yi++) {
						for (int ymu = -1; ymu <= 1; ymu += 2) {
							real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
							if (fabs(yPos-bgY) > ROI*jincr*2.) continue;
							
							real rSquare = (bgX-xPos)*(bgX-xPos)/(iincr*iincr)+ (bgY-yPos)*(bgY-yPos)/(jincr*jincr);
							int bgI = xi*2 + (xmu+1)/2;
							int bgJ = yi*2 + (ymu+1)/2;
							int bgK = zi*2 + (zmu+1)/2;
							int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
							if (rSquare < Rsquare) {
								real weight = exp(-2.302585092994045*rSquare/Rsquare);
								if (logzPos > logheights.front()) {
									bgSpline->solve(uBG.data());
									bgU[bIndex] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(vBG.data());
									bgU[bIndex +1] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(wBG.data());
									bgU[bIndex +2] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(tBG.data());
									bgU[bIndex +3] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(qBG.data());
									bgU[bIndex +4] += weight*(bgSpline->evaluate(logzPos));
									bgSpline->solve(rBG.data());
									bgU[bIndex +5] += weight*(bgSpline->evaluate(logzPos));
									bgWeights[bIndex] += weight;
								} else {
									// Below the spline interpolation
									bgU[bIndex] += weight*uBG.front();
									bgU[bIndex +1] += weight*vBG.front();
									bgU[bIndex +2] += weight*wBG.front();
									bgU[bIndex +3] += weight*tBG.front();
									bgU[bIndex +4] += weight*qBG.front();
									bgU[bIndex +5] += weight*rBG.front();
									bgWeights[bIndex] += weight;
								}								
							}
						}
					}
				}
			}
		}
	}
	
	delete bgSpline;
	logheights.clear();
	uBG.clear();
	vBG.clear();
	wBG.clear();
	tBG.clear();
	qBG.clear();
	rBG.clear();
	
	
	int numbgObs = bgIn.size()*7/11;
	if (numbgObs > 0) {
		// Check interpolation
		for (int zi = 0; zi < (kdim-1); zi++) {	
			for (int zmu = -1; zmu <= 1; zmu += 2) {
				real zPos = kmin + kincr * (zi + (0.5*sqrt(1./3.) * zmu + 0.5));
				
				for (int xi = 0; xi < (idim-1); xi++) {
					for (int xmu = -1; xmu <= 1; xmu += 2) {
						real xPos = imin + iincr * (xi + (0.5*sqrt(1./3.) * xmu + 0.5));
						
						for (int yi = 0; yi < (jdim-1); yi++) {
							for (int ymu = -1; ymu <= 1; ymu += 2) {
								real yPos = jmin + jincr * (yi + (0.5*sqrt(1./3.) * ymu + 0.5));
								int bgI = xi*2 + (xmu+1)/2;
								int bgJ = yi*2 + (ymu+1)/2;
								int bgK = zi*2 + (zmu+1)/2;
								int bIndex = numVars*(idim-1)*2*(jdim-1)*2*bgK + numVars*(idim-1)*2*bgJ +numVars*bgI;
								for (unsigned int var = 0; var < numVars; var++) {
									if (bgWeights[bIndex] != 0) {
										bgU[bIndex +var] /= bgWeights[bIndex];
									} else {
										cout << "Empty background mish at " << xPos << ", " << yPos << ", " << zPos << endl;
									}
								}				
								bgWeights[bIndex] = 0.;
							}
						}
					}
				}
			}
		}	
	} else {
		cout << "No background observations loaded" << endl;
		return 0;
	}
	
	cout << numbgObs << " background observations loaded" << endl;	
	return numbgObs;
}

bool VarDriverXYZ::adjustBackground(const int& bStateSize)
{
	/* Set the minimum filter length to the background resolution, not the analysis resolution
	 to avoid artifacts when running interpolating to small mesoscale grids */
	
	// Load the observations into a vector
	int numbgObs = bgIn.size()*7/11;
	bgObs = new real[numbgObs*14];
	for (int m=0; m < numbgObs*14; m++) bgObs[m] = 0.;
	
	int p = 0;
	for (int m=0; m < bgIn.size(); m+=11) {
		real bgX = bgIn[m];
		real bgY = bgIn[m+1];
		real bgZ = exp(bgIn[m+2]);
		real bgTime = bgIn[m+3];
		if ((bgX < imin) or (bgX > imax) or
			(bgY < jmin) or (bgY > jmax) or
			(bgZ < kmin) or (bgZ > kmax)) {
			numbgObs -= 7;
			continue;
		}
		// Restrict the horizontal domain if we are using the R0 BC
		if (configHash.value("horizontalbc") == "R0") {
			if ((bgX < (imin+iincr)) or (bgX > (imax-iincr)) or
				(bgY < (jmin+jincr)) or (bgY > (jmax-jincr))) 
				continue;
		}
		
		if (configHash.value("verticalbc") == "R0") {
			if ((bgZ < (kmin+kincr)) or (bgZ > (kmax-kincr)))
				continue;
		}
		
		for (unsigned int n = 0; n < numVars; n++) {
			bgObs[p] = bgIn[m+4+n];
			// Error of background = 1
			bgObs[p+1] = 1.;
			bgObs[p+2] = bgX;
			bgObs[p+3] = bgY;
			bgObs[p+4] = bgZ;
			// Null type
			bgObs[p+5] = -1;
			bgObs[p+6] = bgTime;
			bgObs[p+7+n] = 1.;
			p += 14;
		}
	}	
	
	// Adjust the background field to the spline mish
	bgCostXYZ = new CostFunctionXYZ(numbgObs, bStateSize);
	bgCostXYZ->initialize(&configHash, bgU, bgObs);
	/* Set the iteration to zero -- this will prevent writing the background file until after the adjustment
	    which is presumably what you want most of the time. Otherwise, you would not be here */
	int bgIter = 1;
	bgCostXYZ->initState(bgIter);
	bgCostXYZ->minimize();
	
	// Increment the variables
	bgCostXYZ->updateBG();
	bgCostXYZ->finalize();
	
	delete bgCostXYZ;
	delete[] bgObs;
		
	return true;
}


/* Any updates needed for additional analysis iterations go here */
 
void VarDriverXYZ::updateAnalysisParams(const int& iteration)
{
	QString iter;
	iter.setNum(iteration);
	
	QString key = "uerror_" + iter;
	QString val = configHash.value(key);
	configHash.insert("uerror", val);
	
	key = "verror_" + iter;
	val = configHash.value(key);
	configHash.insert("verror", val);
	
	key = "werror_" + iter;
	val = configHash.value(key);
	configHash.insert("werror", val);
	
	key = "terror_" + iter;
	val = configHash.value(key);
	configHash.insert("terror", val);
	
	key = "qverror_" + iter;
	val = configHash.value(key);
	configHash.insert("qverror", val);
	
	key = "rhoerror_" + iter;
	val = configHash.value(key);
	configHash.insert("rhoerror", val);
	
	key = "qrerror_" + iter;
	val = configHash.value(key);
	configHash.insert("qrerror", val);	
	
	key = "mcweight_" + iter;
	val = configHash.value(key);
	configHash.insert("mcweight", val);
	
	key = "xfilter_" + iter;
	val = configHash.value(key);
	configHash.insert("xfilter", val);
	
	key = "yfilter_" + iter;
	val = configHash.value(key);
	configHash.insert("yfilter", val);
	
	key = "zfilter_" + iter;
	val = configHash.value(key);
	configHash.insert("zfilter", val);
	
	key = "x_spline_cutoff_" + iter;
	val = configHash.value(key);
	configHash.insert("x_spline_cutoff", val);
	
	key = "y_spline_cutoff_" + iter;
	val = configHash.value(key);
	configHash.insert("y_spline_cutoff", val);
	
	key = "z_spline_cutoff_" + iter;
	val = configHash.value(key);
	configHash.insert("z_spline_cutoff", val);
	
}

/* This routine validates that all required parameters are present
It currently does not check the validity of a particular parameter, just that it exists */

bool VarDriverXYZ::validateXMLconfig()
{
	
	// Validate the hash -- multiple passes are not validated currently
	QStringList configKeys;
	configKeys << "xmin" << "xmax" << "xincr" <<
	"ymin" << "ymax" << "yincr" <<
	"zmin" << "zmax" << "zincr" <<
	"xfilter" << "yfilter" << "zfilter" <<
	"uerror" << "verror" << "werror" << "terror" << 
	"qverror" << "rhoerror" << "qrerror" << "mcweight" << 
	"radardbz" << "radarvel" << "radarsw" << "radarskip" << "radarstride" << "dynamicstride" <<
	"horizontalbc" << "verticalbc" << "use_dbz_pseudow" <<
	"x_spline_cutoff" << "y_spline_cutoff" << "z_spline_cutoff";
	
	for (int i = 0; i < configKeys.count(); i++) {
		if (!configHash.contains(configKeys.at(i))) {
			cout <<	"No configuration found for <" << configKeys.at(i).toStdString() << "> aborting..." << endl;
			return false;
		}
	}
	return true;
}	
	
	
	
