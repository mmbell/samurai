/*
 *  VarDriverThermo.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriverThermo.h"
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
VarDriverThermo::VarDriverThermo()
: VarDriver3D()
{
	numVars = 7;
    numDerivatives = 4;
    obMetaSize = 7;
}

// Destructor
VarDriverThermo::~VarDriverThermo()
{
}

/* This routine is the main initializer of the analysis */

bool VarDriverThermo::initialize(const QDomElement& configuration)
{
	// Run a 3D vortex background field
	cout << "Initializing SAMURAI Thermo" << endl;
	
	// Parse the XML configuration file
	if (!parseXMLconfig(configuration)) return false;
	
	// Validate the 3D specific parameters
	if (!validateXMLconfig()) return false;
    
    // Validate the run geometry
    if (configHash.value("mode") == "XYZ") {
        runMode = XYZ;
    } else if (configHash.value("mode") == "RTZ") {
        runMode = RTZ;
	} else {
        cout << "Unrecognized run mode " << configHash.value("mode").toStdString() << ", Aborting...\n";
        return false;
    }
    
	// Define the grid dimensions
	imin = configHash.value("i_min").toFloat();
	imax = configHash.value("i_max").toFloat();
	iincr = configHash.value("i_incr").toFloat();
	idim = (int)((imax - imin)/iincr) + 1;
	
	jmin = configHash.value("j_min").toFloat();
	jmax = configHash.value("j_max").toFloat();
	jincr = configHash.value("j_incr").toFloat();
	jdim = (int)((jmax - jmin)/jincr) + 1;
	
	kmin = configHash.value("k_min").toFloat();
	kmax = configHash.value("k_max").toFloat();
	kincr = configHash.value("k_incr").toFloat();
	kdim = (int)((kmax - kmin)/kincr) + 1;
    
	// The recursive filter uses a fourth order stencil to spread the observations, so less than 4 gridpoints will cause a memory fault
	if (idim < 4) {
		cout << "i dimension is less than 4 gridpoints and recursive filter will fail. Aborting...\n";
		return false;
	}
	if (jdim < 4) {
		cout << "j dimension is less than 4 gridpoints and recursive filter will fail. Aborting...\n";
		return false;
	}	
	if (kdim < 4) {
		cout << "k dimension is less than 4 gridpoints and recursive filter will fail. Aborting...\n";
		return false;
	}
	
	// Define the sizes of the arrays we are passing to the cost function
	cout << "iMin\tiMax\tiIncr\tjMin\tjMax\tjIncr\tkMin\tkMax\tkIncr\n";
	cout << imin << "\t" <<  imax << "\t" <<  iincr << "\t";
	cout << jmin << "\t" <<  jmax << "\t" <<  jincr << "\t";
	cout << kmin << "\t" <<  kmax << "\t" <<  kincr << "\n\n";

	int uStateSize = 8*(idim+1)*(jdim+1)*(kdim+1)*(numVars);
	int bStateSize = (idim+2)*(jdim+2)*(kdim+2)*numVars;
	cout << "Physical (mish) State size = " << uStateSize << "\n";
	cout << "Nodal State size = " << bStateSize << ", Grid dimensions:\n";
	
	// Load the BG into a empty vector
	bgU = new real[uStateSize];
	bgWeights = new real[uStateSize];
	for (int i=0; i < uStateSize; i++) {
		bgU[i] = 0.;
		bgWeights[i] = 0.;
	}		
	
	/* Set the data path */
	dataPath.setPath(configHash.value("data_directory"));
	if (dataPath.exists()) {
		dataPath.setFilter(QDir::Files);
		dataPath.setSorting(QDir::Name);
	} else {
		cout << "Can't find data directory: " << configHash.value("data_directory").toStdString() << endl;
		return false;
	}
	
	/* Check to make sure the output path exists */
	QDir outputPath(configHash.value("output_directory"));
	if (!outputPath.exists()) {
		cout << "Can't find output directory: " << configHash.value("output_directory").toStdString() << endl;
		return false;
	}
	
	// Define the Reference state
	QString refSounding = dataPath.absoluteFilePath(configHash.value("ref_state"));
    refstate = new ReferenceState(refSounding);	
	cout << "Reference profile: Z\t\tQv\tRhoa\tRho\tH\tTemp\tPressure\n";
	for (real k = kmin; k < kmax+kincr; k+= kincr) {
		cout << "                   " << k << "\t";
		for (int i = 0; i < 6; i++) {
			real var = refstate->getReferenceVariable(i, k*1000);
			if (i == 0) var = refstate->bhypInvTransform(var);
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
	QTime reftime = QTime::fromString(configHash.value("ref_time"), "hh:mm:ss");
	QString refstring = reftime.toString();
	bool foundref = false;
	for (unsigned int fi = 0; fi < frameVector.size(); fi++) {
		QDateTime frametime = frameVector[fi].getTime();
		if (reftime == frametime.time()) {
			QString tempstring;
			QDate refdate = frametime.date();
			QDateTime unixtime(refdate, reftime, Qt::UTC);
			configHash.insert("ref_lat", tempstring.setNum(frameVector[fi].getLat()));
			configHash.insert("ref_lon", tempstring.setNum(frameVector[fi].getLon()));
			configHash.insert("ref_time", tempstring.setNum(unixtime.toTime_t()));
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
	maxIter = configHash.value("num_iterations").toInt();
    		
	/* Optionally load a set of background estimates and interpolate to the Gaussian mish */
	
	if(!this->readNcFile()) {
				cout << "Reading Nc-File failed ...Exit." << endl;
				return EXIT_FAILURE;
			}
	
	if(!this->loadObsVector()) {
				cout << "Loading ObsVector failed ...Exit." << endl;
				return EXIT_FAILURE;
			}			
	
// AF for now set the background to zero and don't allow any additional observations, i.e. skip loadBackgroundObs, adjustBackground, preProcessMetObs and loadMetObs

	// We are done with the bgWeights, so free up that memory
	delete[] bgWeights;	
	   
	//obCost3D = new CostFunctionThermo(obVector.size(), bStateSize);
	//obCost3D->initialize(&configHash, bgU, obs, refstate);
	
	// If we got here, then everything probably went OK!
	return true;
}

/* This routine drives the CostFunction minimization
 There is support for an outer loop to change the background
 error covariance or update non-linear observation operators */

bool VarDriverThermo::run()
{
	//int iter=1;
	//while (iter <= maxIter) {
		//cout << "Outer Loop Iteration: " << iter << endl;
		//obCost3D->initState(iter);
		//obCost3D->minimize();
		//obCost3D->updateBG();
		//iter++;
		
		// Optionally update the analysis parameters for an additional iteration
		//updateAnalysisParams(iter);
	//}	
	
	return true;
	
}

/* Clean up all that allocated memory */

bool VarDriverThermo::finalize()
{
	//obCost3D->finalize();
	delete[] obs;
	delete[] bgU;
	//delete obCost3D;
	delete refstate;
	return true;
}


bool VarDriverThermo::readNcFile()
{
  cout << "Read in NetCDF File " << endl;
  QString ncFileName = "/bora/rita2005/Samurai/wrf/20050920_19.nc";
  if (ncFile.readNetCDF(ncFileName.toAscii().data()) != 0) {
	  cout << "Error reading NetCDF file\n";
	  exit(1);
	}
  
  double a = ncFile.calc_A(4,5,6);
  std::cout << "A is: " << a << "\n";

  double b = ncFile.calc_B(4,5,6);
  std::cout << "B is: " << b << "\n";
  
  double c = ncFile.calc_C(4,5,6);
  std::cout << "C is: " << c << "\n";

  return true;
  
}

bool VarDriverThermo::loadObsVector()
{
  cout << "Load Obs Vector " << endl;
  return true;
  
}
