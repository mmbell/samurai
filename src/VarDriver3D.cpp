/*
 *  VarDriver3D.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriver3D.h"
#include "Dorade.h"
#include <iterator>
#include <fstream>
#include <cmath>
#include <QTextStream>
#include <QFile>
#include <QVector>
#include <iomanip>
#include "RecursiveFilter.h"
#include "samurai.h"

// TODO debug

// Constructor
VarDriver3D::VarDriver3D()
: VarDriver()
{
  numVars = 7;
  numDerivatives = 4;
  bkgdAdapter = NULL;
}

// Destructor
VarDriver3D::~VarDriver3D()
{
}

// Initialize from parsed xml

bool VarDriver3D::initialize(const QDomElement& configuration)
{
  // Run a 3D vortex background field
  cout << "Initializing SAMURAI 3D" << endl;

  // Parse the XML configuration file
  if (!parseXMLconfig(configuration)) return false;
  
  return validateDriver();
}

// Initialize from a passed structure

bool VarDriver3D::initialize(const samurai_config &configSam)
{
  // Run a 3D vortex background field
  cout << "Initializing SAMURAI 3D" << endl;

  // Parse the Samurai configuration structure passed from COAMPS
  if (! parseSamuraiConfig(configSam)) return false;
  
  return validateDriver();
}

// Validate and finish initializing the driver

bool VarDriver3D::validateDriver()
{
  // Make sure the config has all the keys we need
  
  if ( ! validateConfig())
    return false;
  
  // std::cout << "==== Content of configHash" << std::endl;
  // dump_hash(configHash);
  
  // Validate the run geometry
  if (configHash.value("mode") == "XYZ") {
    runMode = XYZ;
  } else if (configHash.value("mode") == "RTZ") {
    runMode = RTZ;
  } else {
    cout << "Unrecognized run mode " << configHash.value("mode").toStdString() << ", Aborting...\n";
    return false;
  }

  // Set the projection (to be used by the cost functions)
  projection.setProjection(projectionFromConfig());

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

  // TODO
  // Centers and Met Observations are tightly coupled.
  //
  // So if we allow different centers between each runs, preProcessMetObs() or loadMetObs()
  // have to be called again.
  //
  // What else is dependent on the data structures created by readFrameCenters() and the *MetObs() ?
  //
  // Are all the Met Obs still available, or only the ones that were matched to the center time frames?

  // With coamps, grid dimensions come from the run arguments, not the fixed confid
  // So only do validation that depends on grid dimensions if fixedGrid is set.

  if (fixedGrid) {   // If runGrid, this will be done in the run(......) call
    if ( ! validateFixedGrid() )
      return false;

    // Set the background Obs adapter to get data from a file (if config says so)
    QString loadBG = configHash.value("load_background");
    if (loadBG == "true")
      bkgdAdapter = new BkgdStream(dataPath.absoluteFilePath("samurai_Background.in").toLatin1().data());
    
    return gridDependentInit();
  }
  
  // If we got here, then everything probably went OK!
  return true;
}

// This used to be in validateDriver.
// All the grid dependent initialization is performed here.
// Prerequisit: Grid dimensions have been set and verified.

bool VarDriver3D::gridDependentInit()
{
  // Define the Reference state

  QString refSounding = dataPath.absoluteFilePath(configHash.value("ref_state"));
  refstate = new ReferenceState(refSounding);
  // cout << "Reference profile: Z\t\tQv\tRhoa\tRho\tH\tTemp\tPressure\n";
  for (real k = kmin; k < kmax + kincr; k += kincr) {
    // cout << "                   " << k << "\t";
    for (int i = 0; i < 6; i++) {
      real var = refstate->getReferenceVariable(i, k*1000);
      if (i == 0) var = refstate->bhypInvTransform(var);
      // cout << setw(9) << setprecision(4)  << var << "\t";
    }
    // cout << "\n";
  }
  cout << setprecision(9);
  
  // Set the maximum number of iterations to the multipass reduction factor
  //  Multiple outer loops will reduce the cutoff wavelengths and background error variance
  
  maxIter = configHash.value("num_iterations").toInt();
    
  // Read in the Frame centers (if runGrid, user is responsible for managing center
  // structure instead)
  //
  // Ideally, create a time-based spline from limited center fixes here
  // but just load 1 second centers into vector for now

  if ( fixedGrid) {
    readFrameCenters();
    if ( ! findReferenceCenter() ) {
      cout << "Error finding reference time, please check date and time in XML file\n";
      return false;
    }
  }
  
  // Optionally load a set of background coefficients directly
  QString loadBGcoeffs = configHash.value("load_bg_coefficients");

  if (loadBGcoeffs == "true")
    if (! loadBackgroundCoeffs() )
      return false;

  // Optionally load a set of background estimates and interpolate to the Gaussian mish
  
  int numbgObs = 0;
  if(bkgdAdapter != NULL) {
    numbgObs = loadBackgroundObs();
    if (numbgObs < 0) {
      cout << "Error loading background Obs\n";
      return false;
    }
  }
  
  // Optionally adjust the interpolated background to satisfy mass continuity
  // and match the supplied points exactly. In essence, do a SAMURAI analysis using
  // the background estimates as "observations"

  QString adjustBG = configHash.value("adjust_background");
  if ((adjustBG == "true") and numbgObs) {
    if ( ! adjustBackground()) {
      cout << "Error adjusting background\n";
      return false;
    }
  }

  if ( ! loadMetObs() )
    return false;
  initObCost3D();
  return true;
}

// Find the center that matches the ref_time

bool VarDriver3D::findReferenceCenter()
{
  QTime reftime = QTime::fromString(configHash.value("ref_time"), "hh:mm:ss");
  QString refstring = reftime.toString();

  for (unsigned int fi = 0; fi < frameVector.size(); fi++) {
    QDateTime frametime = frameVector[fi].getTime();
    if (reftime == frametime.time()) {
      QString tempstring;
      QDate refdate = frametime.date();
      QDateTime unixtime(refdate, reftime, Qt::UTC);
      configHash.insert("ref_lat", tempstring.setNum(frameVector[fi].getLat()));
      configHash.insert("ref_lon", tempstring.setNum(frameVector[fi].getLon()));
      // Overwrite hh:mm:ss string with epoch
      configHash.insert("ref_time", tempstring.setNum(unixtime.toTime_t()));
      cout << "Found matching reference time " << refstring.toStdString()
	   << " at " << frameVector[fi].getLat() << ", " << frameVector[fi].getLon() << "\n";

      return true;
    }
  }
  return true;
}

//
// Initialize background observations and run
//

bool VarDriver3D::run(int nx, int ny, int nsigma,

		      // ----- new -----
		      char cdtg[10],	// "12Z oct 4 2015 -> "2015100412"
		      int delta,	// delta * iter past cdtg
		      int iter,
		      float imin, float imax, float iincr, // used to come from config
		      float jmin, float jmax, float jincr,
		      // ----- new -----
		      
		      float *sigmas,
		      float *latitude, // 2D arrays
		      float *longitude,
		      float *u1,	// 3D array (nx, ny, nsigma)  
		      float *v1,
		      float *w1,
		      float *th1,
		      float *p1,
		      // These are output values
		      float *usam,	// 3D array
		      float *vsam,
		      float *wsam,
		      float *thsam,
		      float *psam)
{
  fillRunCenters(cdtg, delta, iter, *latitude, *longitude);
  
  if (! validateRunGrid(nx, ny, nsigma,
			imin, imax, iincr,
			jmin, jmax, jincr, sigmas) )
    return false;


  // Clean up from previous run if needed
  
  if(bkgdAdapter != NULL)
    delete bkgdAdapter;
  
  // Set the background obs adapter to get data from the passed arrays
  
  bkgdAdapter = new BkgdArray(nx, ny, nsigma,
			      cdtg, delta, iter,
			      sigmas,
			      latitude, longitude,
			      u1, v1, w1, th1, p1);

  if (! gridDependentInit() )
    return false;
#if 0  
  // convert the time info
  // cdtg:  current date-time group. oc. 12th 2015 at 12Z -> "2015100412"
  // delta: model time step on coarse grid
  // iter: current iteration
  // delta * iter: model time past cdtg
  
  if ( ! loadRunBackgroundObs(nx, ny, nsigma,
			      cdtg, delta, iter,
			      sigmas, latitude, longitude,
			      u1, v1, w1, th1, p1))
    return false;
#endif
  
  if( run() )
    return obCost3D->copyResults(nx, ny, nsigma, usam, vsam, wsam, thsam, psam);
  return false;
}


/* This routine drives the CostFunction minimization
 There is support for an outer loop to change the background
 error covariance or update non-linear observation operators */

bool VarDriver3D::run()
{
	int iter=1;
	while (iter <= maxIter) {
		if (iter < maxIter) {
			configHash["save_mish"] = "true";
		} else {
			configHash["save_mish"] = "false";
		}
		cout << "Outer Loop Iteration: " << iter << endl;
		obCost3D->initState(iter);
		obCost3D->minimize();
		obCost3D->updateBG();
		iter++;

		// Optionally update the analysis parameters for an additional iteration
		updateAnalysisParams(iter);
	}

	return true;

}

/* Clean up all that allocated memory */

bool VarDriver3D::finalize()
{
	obCost3D->finalize();
	delete[] obs;
	delete[] bgU;
	delete obCost3D;
	delete refstate;
	return true;
}

/* Pre-process the observations into a single vector
 On the wishlist is some integrated QC here other than just spatial thresholding */

bool VarDriver3D::preProcessMetObs()
{

  vector<real> rhoP;

  // Convert the bg dBZ back to Z for further processing with real radar data

  if ((configHash.value("qr_variable") == "dbz") and
      (configHash.value("load_background") == "true") and
      (configHash.value("adjust_background") == "false")) {
    for (int ki = -1; ki < (kdim); ki++) {
      for (int kmu = -1; kmu <= 1; kmu += 2) {
	for (int ii = -1; ii < (idim); ii++) {
	  for (int imu = -1; imu <= 1; imu += 2) {
	    for (int ji = -1; ji < (jdim); ji++) {
	      for (int jmu = -1; jmu <= 1; jmu += 2) {
		int bgI = (ii+1)*2 + (imu+1)/2;
		int bgJ = (ji+1)*2 + (jmu+1)/2;
		int bgK = (ki+1)*2 + (kmu+1)/2;
		int bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
		real dBZ = bgU[bIndex +6]*10. - 35.;
		real ZZ = pow(10.0,(dBZ*0.1));
		bgU[bIndex +6] = ZZ;
		bgWeights[bIndex] = 1.0;
	      }
	    }
	  }
	}
      }
    }
  }

  // Geographic functions
  //GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();

  real referenceLon = configHash.value("ref_lon").toFloat();

  // Find the zero C line using Newton's method

  real zeroClevel = 273.15;
  real height = 5000;
  real tmin = 1e34;
  int iter = 0;

  while ((fabs(tmin) > 0.1) and (iter < 5000)) {
    real t = refstate->getReferenceVariable(ReferenceVariable::tempref, height) - zeroClevel;
    real tprime = (refstate->getReferenceVariable(ReferenceVariable::tempref, height+500.)
		   - refstate->getReferenceVariable(ReferenceVariable::tempref, height-500.))/1000.;
    if (tprime != 0) {
      height = height - t/tprime;
      tmin = t;
    }
    iter++;
  }
  zeroClevel = height;
  cout << "Found zero C level at " << zeroClevel << " based on reference sounding" << endl;

  QStringList filenames = dataPath.entryList();
  int processedFiles = 0;
  QList<MetObs>* metData = new QList<MetObs>;
  cout << "Found " << filenames.size() << " data files to read..." << endl;
  
  for (int i = 0; i < filenames.size(); ++i) {
    metData->clear();
    QString file = filenames.at(i);
    QStringList fileparts = file.split(".");
    if (fileparts.isEmpty()) {
      cout << "Unknown file! " << file.toLatin1().data() << endl;
      continue;
    }
    QString suffix = fileparts.last();
    QString prefix = fileparts.first();
    if (prefix == "swp") {
      // Switch it to suffix
      suffix = "swp";
    }
    cout << "Processing " << file.toLatin1().data() << " of type " << suffix.toLatin1().data() << endl;
    QFile metFile(dataPath.filePath(file));

    // Read different types of files
    if (! read_met_obs_file(dataSuffix.value(suffix), metFile, metData))
      continue;

    processedFiles++;

    int obsProblem = 0;
    int timeProblem = 0;
    int coordProblem = 0;
    int domainProblem = 0;
    int radiusProblem = 0;
    
    // Process the metObs into Observations
    QDateTime startTime = frameVector.front().getTime();
    QDateTime endTime = frameVector.back().getTime();
    int prevobs = obVector.size();
    
    for (int i = 0; i < metData->size(); ++i) {

      // Make sure the ob is within the time limits
      MetObs metOb = metData->at(i);
      QDateTime obTime = metOb.getTime();
      QString obstring = obTime.toString(Qt::ISODate);
      QString tcstart = startTime.toString(Qt::ISODate);
      QString tcend = endTime.toString(Qt::ISODate);
      if ((obTime < startTime) or (obTime > endTime)) {
	obsProblem++;
	continue;
      }
      int fi = startTime.secsTo(obTime);
      if ((fi < 0) or (fi > (int)frameVector.size())) {
	cout << "Time problem with observation " << fi << endl;
	timeProblem++;
	continue;
      }
      real Um = frameVector[fi].getUmean();
      real Vm = frameVector[fi].getVmean();

      // Get the X, Y & Z
      real tcX, tcY, metX, metY;
      if ((metOb.getLat() == -999) or (metOb.getLat() == -999)) {
	coordProblem++;
	continue;
      }
      projection.Forward(referenceLon, frameVector[fi].getLat() , frameVector[fi].getLon() , tcX, tcY);
      projection.Forward(referenceLon, metOb.getLat() , metOb.getLon() , metX, metY);
      real obX = (metX - tcX)/1000.;
      real obY = (metY - tcY)/1000.;
      real heightm = metOb.getAltitude();
      real obZ = heightm/1000.;
      real obRadius = sqrt(obX*obX + obY*obY);
      real obTheta = 180.0 * atan2(obY, obX) / Pi;
      if (configHash.value("allow_negative_angles") != "true") {
	if (obTheta < 0) obTheta += 360.0;
      }

        // Make sure the ob is in the domain
      if (runMode == XYZ) {
	if ((obX < imin) or (obX > imax) or
	    (obY < jmin) or (obY > jmax) or
	    (obZ < kmin) or (obZ > kmax)) {
	  domainProblem++;
	  continue;
	}
      } else if (runMode == RTZ) {
	if ((obRadius < imin) or (obRadius > imax) or
	    (obTheta < jmin) or (obTheta > jmax) or
	    (obZ < kmin) or (obZ > kmax)) {
	  domainProblem++;
	  continue;
	}
	
	if (obRadius == 0.0) {
	  radiusProblem++;
	  continue;
	}
      }
      // Create an observation and set its basic info
      Observation varOb;
      varOb.setCartesianX(obX);
      varOb.setCartesianY(obY);
      varOb.setRadius(obRadius);
      varOb.setTheta(obTheta);
      varOb.setAltitude(obZ);
      varOb.setTime(obTime.toTime_t());

      // Reference states
      real rhoBar = refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm);
      real qBar = refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm);
      real tBar = refstate->getReferenceVariable(ReferenceVariable::tempref, heightm);

      // Initialize the weights
      for (unsigned int var = 0; var < numVars; ++var) {
	for (unsigned int d = 0; d < numDerivatives; ++d) {
	  varOb.setWeight(0.0, var, d);
	}
      }
      
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
	  if (runMode == XYZ) {
	    rhou = rho*(u - Um);
	  } else if (runMode == RTZ) {
	    rhou = rho*((u - Um)*obX + (v - Vm)*obY)/obRadius;
	  }
	  //cout << "RhoU: " << rhou << endl;
	  varOb.setOb(rhou);
	  varOb.setError(configHash.value("dropsonde_rhou_error").toFloat());
	  obVector.push_back(varOb);
	  
	  varOb.setWeight(0., 0);

	  varOb.setWeight(1., 1);
	  if (runMode == XYZ) {
	    rhov = rho*(v - Vm);
	  } else if (runMode == RTZ) {
	    rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	  }
	  varOb.setOb(rhov);
	  varOb.setError(configHash.value("dropsonde_rhov_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 1);

	}
	if ((w != -999) and (rho != -999)) {
	  // rho w 1.5 m/s error
	  varOb.setWeight(1., 2);
	  rhow = rho*w;
	  varOb.setOb(rhow);
	  varOb.setError(configHash.value("dropsonde_rhow_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 2);
	}
	if (tempk != -999) {
	  // temperature 1 K error
	  varOb.setWeight(1., 3);
	  varOb.setOb(tempk - tBar);
	  varOb.setError(configHash.value("dropsonde_tempk_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 3);
	}
	if (qv != -999) {
	  // Qv 0.5 g/kg error
	  varOb.setWeight(1., 4);
	  qv = refstate->bhypTransform(qv);
	  varOb.setOb(qv-qBar);
	  varOb.setError(configHash.value("dropsonde_qv_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 4);
	}
	if (rhoa != -999) {
	  // Rho prime .1 kg/m^3 error
	  varOb.setWeight(1., 5);
	  varOb.setOb((rhoa-rhoBar)*100);
	  varOb.setError(configHash.value("dropsonde_rhoa_error").toFloat());
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
	  if (runMode == XYZ) {
	    rhou = rho*(u - Um);
	  } else if (runMode == RTZ) {
	    rhou = rho*((u - Um)*obX + (v - Vm)*obY)/obRadius;
	  }
	  varOb.setOb(rhou);
	  varOb.setError(configHash.value("flightlevel_rhou_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 0);

	  varOb.setWeight(1., 1);
	  if (runMode == XYZ) {
	    rhov = rho*(v - Vm);
	  } else if (runMode == RTZ) {
	    rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	  }
	  varOb.setOb(rhov);
	  varOb.setError(configHash.value("flightlevel_rhov_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 1);
	}
	if ((w != -999) and (rho != -999)) {
	  // rho w 1 dm/s error
	  varOb.setWeight(1., 2);
	  rhow = rho*w;
	  varOb.setOb(rhow);
	  varOb.setError(configHash.value("flightlevel_rhow_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 2);
	}
	if (tempk != -999) {
	  // temperature 1 K error
	  varOb.setWeight(1., 3);
	  varOb.setOb(tempk - tBar);
	  varOb.setError(configHash.value("flightlevel_tempk_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 3);
	}
	if (qv != -999) {
	  // Qv 0.5 g/kg error
	  varOb.setWeight(1., 4);
	  qv = refstate->bhypTransform(qv);
	  varOb.setOb(qv-qBar);
	  varOb.setError(configHash.value("flightlevel_qv_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 4);
	}
	if (rhoa != -999) {
	  // Rho prime .1 kg/m^3 error
	  varOb.setWeight(1., 5);
	  varOb.setOb((rhoa-rhoBar)*100);
	  varOb.setError(configHash.value("flightlevel_rhoa_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 5);
	}

	break;

      case (MetObs::insitu):
	varOb.setType(MetObs::insitu);
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
	  if (runMode == XYZ) {
	    rhou = rho*(u - Um);
	  } else if (runMode == RTZ) {
	    rhou = rho*((u - Um)*obX + (v - Vm)*obY)/obRadius;
	  }
	  //cout << "RhoU: " << rhou << endl;
	  varOb.setOb(rhou);
	  varOb.setError(configHash.value("insitu_rhou_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 0);

	  varOb.setWeight(1., 1);
	  if (runMode == XYZ) {
	    rhov = rho*(v - Vm);
	  } else if (runMode == RTZ) {
	    rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	  }
	  varOb.setOb(rhov);
	  varOb.setError(configHash.value("insitu_rhov_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 1);

	}
	if ((w != -999) and (rho != -999)) {
	  // rho w 1.5 m/s error
	  varOb.setWeight(1., 2);
	  rhow = rho*w;
	  varOb.setOb(rhow);
	  varOb.setError(configHash.value("insitu_rhow_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 2);
	}
	if (tempk != -999) {
	  // temperature 1 K error
	  varOb.setWeight(1., 3);
	  varOb.setOb(tempk - tBar);
	  varOb.setError(configHash.value("insitu_tempk_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 3);
	}
	if (qv != -999) {
	  // Qv 0.5 g/kg error
	  varOb.setWeight(1., 4);
	  qv = refstate->bhypTransform(qv);
	  varOb.setOb(qv-qBar);
	  varOb.setError(configHash.value("insitu_qv_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 4);
	}
	if (rhoa != -999) {
	  // Rho prime .1 kg/m^3 error
	  varOb.setWeight(1., 5);
	  varOb.setOb((rhoa-rhoBar)*100);
	  varOb.setError(configHash.value("insitu_rhoa_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 5);
	}

	break;

      case (MetObs::mtp):
	varOb.setType(MetObs::mtp);
	rhoa = metOb.getDryDensity(); // Pressure/density is from dry air only?
	tempk = metOb.getTemperature();
	if (tempk != -999) {
	  // temperature 1 K error
	  varOb.setWeight(1., 3);
	  varOb.setOb(tempk - tBar);
	  varOb.setError(metOb.getTemperatureError() + configHash.value("mtp_tempk_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 3);
	}
	if (rhoa != -999) {
	  // Rho prime .1 kg/m^3 error
	  varOb.setWeight(1., 5);
	  varOb.setOb((rhoa-rhoBar)*100);
	  varOb.setError(configHash.value("mtp_rhoa_error").toFloat());
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
	//varOb.setWeight(1., 0);
	varOb.setWeight(1., 1);
	varOb.setOb(wspd);
	varOb.setError(configHash.value("sfmr_windspeed_error").toFloat());
	obVector.push_back(varOb);
	break;

      case (MetObs::qscat):
	varOb.setType(MetObs::qscat);
	u = metOb.getCartesianUwind();
	v = metOb.getCartesianVwind();
	if (u != -999) {
	  // rho u 1 m/s error
	  // Multiply by rho later from grid values
	  varOb.setWeight(1., 0);
	  if (runMode == XYZ) {
	    rhou = (u - Um);
	  } else if (runMode == RTZ) {
	    rhou = ((u - Um)*obX + (v - Vm)*obY)/obRadius;
	  }
	  //cout << "RhoU: " << rhou << endl;
	  varOb.setOb(rhou);
	  varOb.setError(configHash.value("qscat_rhou_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 0);

	  varOb.setWeight(1., 1);
	  // Multiply by rho later from grid values
	  if (runMode == XYZ) {
	    rhov = (v - Vm);
	  } else if (runMode == RTZ) {
	    rhov = (-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	  }
	  varOb.setOb(rhov);
	  varOb.setError(configHash.value("qscat_rhov_error").toFloat());
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
	  // Multiply by rho later from grid values
	  varOb.setWeight(1., 0);
	  if (runMode == XYZ) {
	    rhou = (u - Um);
	  } else if (runMode == RTZ) {
	    rhou = ((u - Um)*obX + (v - Vm)*obY)/obRadius;
	  }
	  //cout << "RhoU: " << rhou << endl;
	  varOb.setOb(rhou);
	  varOb.setError(configHash.value("ascat_rhou_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 0);

	  varOb.setWeight(1., 1);
	  // Multiply by rho later from grid values
	  if (runMode == XYZ) {
	    rhov = (v - Vm);
	  } else if (runMode == RTZ) {
	    rhov = (-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	  }
	  varOb.setOb(rhov);
	  varOb.setError(configHash.value("ascat_rhov_error").toFloat());
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
	  // Multiply by rho later from grid values
	  varOb.setWeight(1., 0);
	  if (runMode == XYZ) {
	    rhou = (u - Um);
	  } else if (runMode == RTZ) {
	    rhou = ((u - Um)*obX + (v - Vm)*obY)/obRadius;
	  }
	  //cout << "RhoU: " << rhou << endl;
	  varOb.setOb(rhou);
	  varOb.setError(configHash.value("amv_rhou_error").toFloat());
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 0);

	  varOb.setWeight(1., 1);
	  // Multiply by rho later from grid values
	  if (runMode == XYZ) {
	    rhov = (v - Vm);
	  } else if (runMode == RTZ) {
	    rhov = (-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	  }
	  varOb.setOb(rhov);
	  varOb.setError(configHash.value("amv_rhov_error").toFloat());
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
	  real uWgt, vWgt;
	  if (runMode == XYZ) {
	    uWgt = sin(az)*cos(el);
	    vWgt = cos(az)*cos(el);
	  } else if (runMode == RTZ) {
	    uWgt = (obX*sin(az)*cos(el) + obY*cos(az)*cos(el))/obRadius;
	    vWgt = (obX*cos(az)*cos(el) - obY*sin(az)*cos(el))/obRadius;
	  }
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
	  real DopplerError = metOb.getSpectrumWidth()*configHash.value("lidar_sw_error").toFloat()
	    + log(configHash.value("lidar_power_error").toFloat()/db);
	  if (DopplerError < configHash.value("lidar_min_error").toFloat())
	    DopplerError = configHash.value("lidar_min_error").toFloat();
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
	  real uWgt, vWgt;
	  if (runMode == XYZ) {
	    uWgt = sin(az)*cos(el);
	    vWgt = cos(az)*cos(el);
	  } else if (runMode == RTZ) {
	    uWgt = (obX*sin(az)*cos(el) + obY*cos(az)*cos(el))/obRadius;
	    vWgt = (obX*cos(az)*cos(el) - obY*sin(az)*cos(el))/obRadius;
	  }
	  real wWgt = sin(el);
	  // Restrict to horizontal component only
	  if (configHash.value("horizontal_radar_appx") == "true") wWgt = 0;

	  // Fall speed
	  real Z = metOb.getReflectivity();
	  real w_term = 0.0;
	  real ZZ = -999.0;
	  if (Z > -999.0) {
	    real H = metOb.getAltitude();
	    ZZ=pow(10.0,(Z*0.1));
	    real melting_zone = 1000 * configHash.value("melting_zone_width").toFloat();
	    real hlow= zeroClevel;
	    real hhi= hlow + melting_zone;

	    /* density correction term (rhoo/rho)*0.45
	       0.45 density correction from Beard (1985, JOAT pp 468-471) */
	    real rho = refstate->getReferenceVariable(ReferenceVariable::rhoref, H);
	    real rhosfc = refstate->getReferenceVariable(ReferenceVariable::rhoref, 0.);
	    real DCOR = pow((rhosfc/rho),(real)0.45);

	    // The snow relationship (Atlas et al., 1973) --- VT=0.817*Z**0.063  (m/s)
	    real VTS=-DCOR * (0.817*pow(ZZ,(real)0.063));

	    // The rain relationship (Joss and Waldvogel,1971) --- VT=2.6*Z**.107 (m/s) */
	    real VTR=-DCOR * (2.6*pow(ZZ,(real).107));

	    /* Test if height is in the transition region between SNOW and RAIN
	       defined as hlow in km < H < hhi in km
	       if in the transition region do a linear weight of VTR and VTS */
	    real mixed_dbz = configHash.value("mixed_phase_dbz").toFloat();
	    real rain_dbz = configHash.value("rain_dbz").toFloat();
	    if ((Z > mixed_dbz) and
		(Z <= rain_dbz)) {
	      real WEIGHTR=(Z-mixed_dbz)/(rain_dbz - mixed_dbz);
	      real WEIGHTS=1.-WEIGHTR;
	      VTS=(VTR*WEIGHTR+VTS*WEIGHTS)/(WEIGHTR+WEIGHTS);
	    } else if (Z > rain_dbz) {
	      VTS=VTR;
	    }
	    w_term=VTR*(hhi-H)/melting_zone + VTS*(H-hlow)/melting_zone;
	    if (H < hlow) w_term=VTR;
	    if (H > hhi) w_term=VTS;
	  }
	  real VR = metOb.getRadialVelocity();
	  if (VR != -999.0) {
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
	    real DopplerError = fabs(wWgt)*configHash.value("radar_fallspeed_error").toFloat();
	    if (metOb.getSpectrumWidth() != -999.0) {
	      DopplerError += metOb.getSpectrumWidth()*configHash.value("radar_sw_error").toFloat();
	    }
	    if (DopplerError < configHash.value("radar_min_error").toFloat())
	      DopplerError = configHash.value("radar_min_error").toFloat();
	    varOb.setError(DopplerError);
	    varOb.setOb(Vdopp);

	    real maxel = configHash.value("max_radar_elevation").toFloat();
	    if (!maxel) maxel = 90.0;
	    if (fabs(metOb.getElevation() <= maxel))
	      obVector.push_back(varOb);

	    varOb.setWeight(0., 0);
	    varOb.setWeight(0., 1);
	    varOb.setWeight(0., 2);
	  }

	  // Reflectivity observations
	  QString gridref = configHash.value("qr_variable");
	  real qr = 0.;
	  if (ZZ > 0) {
	    if (gridref == "qr") {
	      // Do the gridding as part of the variational synthesis using Z-M relationships
	      // Z-M relationships from Gamache et al (1993) JAS
	      real H = metOb.getAltitude();
	      real melting_zone = 1000 * configHash.value("melting_zone_width").toFloat();
	      real hlow= zeroClevel;
	      real hhi= hlow + melting_zone;
	      real rainmass = pow(ZZ/14630.,(real)0.6905);
	      real icemass = pow(ZZ/670.,(real)0.5587);
	      real mixed_dbz = configHash.value("mixed_phase_dbz").toFloat();
	      real rain_dbz = configHash.value("rain_dbz").toFloat();
	      if ((Z > mixed_dbz) and
		  (Z <= rain_dbz)) {
		real WEIGHTR=(Z-mixed_dbz)/(rain_dbz - mixed_dbz);
		real WEIGHTS=1.-WEIGHTR;
		icemass=(rainmass*WEIGHTR+icemass*WEIGHTS)/(WEIGHTR+WEIGHTS);
	      } else if (Z > 30) {
		icemass=rainmass;
	      }

	      real precipmass = rainmass*(hhi-H)/melting_zone + icemass*(H-hlow)/melting_zone;
	      if (H < hlow) precipmass = rainmass;
	      if (H > hhi) precipmass = icemass;
	      qr = refstate->bhypTransform(precipmass/rhoBar);

	      //Include an observation of this quantity in the variational synthesis
	      varOb.setOb(qr);
	      varOb.setWeight(1., 6);
	      varOb.setError(1.0);
	      obVector.push_back(varOb);

	    } else if (gridref == "dbz") {
	      qr = ZZ;
	      /* Include an observation of this quantity in the variational synthesis
		 varOb.setOb(qr);
		 varOb.setWeight(1., 6);
		 varOb.setError(1.0);
		 obVector.push_back(varOb); */

	    }

	    // Do a Exponential & power weighted interpolation of the reflectivity/qr in a grid box
	    real iROI = configHash.value("i_reflectivity_roi").toFloat() / iincr;
	    real jROI = configHash.value("j_reflectivity_roi").toFloat() / jincr;
	    real kROI = configHash.value("k_reflectivity_roi").toFloat() / kincr;
	    real Rsquare = (iincr*iROI)*(iincr*iROI) + (jincr*jROI)*(jincr*jROI) + (kincr*kROI)*(kincr*kROI);
	    for (int ki = 0; ki < (kdim-1); ki++) {
	      for (int kmu = -1; kmu <= 1; kmu += 2) {
		real kPos = kmin + kincr * (ki + (0.5*sqrt(1./3.) * kmu + 0.5));
		if (fabs(kPos-obZ) > kincr*kROI*2.) continue;
		for (int ii = 0; ii < (idim-1); ii++) {
		  for (int imu = -1; imu <= 1; imu += 2) {
		    real iPos = imin + iincr * (ii + (0.5*sqrt(1./3.) * imu + 0.5));
		    if (runMode == XYZ) {
		      if (fabs(iPos-obX) > iincr*iROI*2.) continue;
		    } else if (runMode == RTZ) {
		      if (fabs(iPos-obRadius) > iincr*iROI*2.) continue;
		    }
		    for (int ji = 0; ji < (jdim-1); ji++) {
		      for (int jmu = -1; jmu <= 1; jmu += 2) {
			real jPos = jmin + jincr * (ji + (0.5*sqrt(1./3.) * jmu + 0.5));
			real rSquare = 0.0;
			if (runMode == XYZ) {
			  if (fabs(jPos-obY) > jincr*jROI*2.) continue;
			  rSquare = (obX-iPos)*(obX-iPos) + (obY-jPos)*(obY-jPos) + (obZ-kPos)*(obZ-kPos);
			} else if (runMode == RTZ) {
			  real dTheta = fabs(jPos-obTheta);
			  if (dTheta > 360.) dTheta -= 360.;
			  if (dTheta > jincr*jROI*2.) continue;
			  rSquare = (obRadius-iPos)*(obRadius-iPos) + (dTheta)*(dTheta) + (obZ-kPos)*(obZ-kPos);
			}
			// Add one extra index to account for buffer zone in analysis
			int bgI = (ii+1)*2 + (imu+1)/2;
			int bgJ = (ji+1)*2 + (jmu+1)/2;
			int bgK = (ki+1)*2 + (kmu+1)/2;
			int bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
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
	  }
	  break;
	}

      case(MetObs::mesonet):
	{
	  varOb.setType(MetObs::mesonet);
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
	    if (runMode == XYZ) {
	      rhou = rho*(u - Um);
	    } else if (runMode == RTZ) {
	      rhou = rho*((u - Um)*obX + (v - Vm)*obY)/obRadius;
	    }
	    //cout << "RhoU: " << rhou << endl;
	    varOb.setOb(rhou);
	    varOb.setError(configHash.value("mesonet_rhou_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 0);

	    varOb.setWeight(1., 1);
	    if (runMode == XYZ) {
	      rhov = rho*(v - Vm);
	    } else if (runMode == RTZ) {
	      rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	    }
	    varOb.setOb(rhov);
	    varOb.setError(configHash.value("mesonet_rhov_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 1);

	  }
	  if ((w != -999) and (rho != -999)) {
	    // rho w 1.5 m/s error
	    varOb.setWeight(1., 2);
	    rhow = rho*w;
	    varOb.setOb(rhow);
	    varOb.setError(configHash.value("mesonet_rhow_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 2);
	  }
	  if (tempk != -999) {
	    // temperature 1 K error
	    varOb.setWeight(1., 3);
	    varOb.setOb(tempk - tBar);
	    varOb.setError(configHash.value("mesonet_tempk_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 3);
	  }
	  if (qv != -999) {
	    // Qv 0.5 g/kg error
	    varOb.setWeight(1., 4);
	    qv = refstate->bhypTransform(qv);
	    varOb.setOb(qv-qBar);
	    varOb.setError(configHash.value("mesonet_qv_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 4);
	  }
	  if (rhoa != -999) {
	    // Rho prime .1 kg/m^3 error
	    varOb.setWeight(1., 5);
	    varOb.setOb((rhoa-rhoBar)*100);
	    varOb.setError(configHash.value("mesonet_rhoa_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 5);
	  }
	  break;
	}

      case(MetObs::aeri):
	{
	  varOb.setType(MetObs::aeri);
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
	    if (runMode == XYZ) {
	      rhou = rho*(u - Um);
	    } else if (runMode == RTZ) {
	      rhou = rho*((u - Um)*obX + (v - Vm)*obY)/obRadius;
	    }
	    //cout << "RhoU: " << rhou << endl;
	    varOb.setOb(rhou);
	    varOb.setError(configHash.value("aeri_rhou_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 0);

	    varOb.setWeight(1., 1);
	    if (runMode == XYZ) {
	      rhov = rho*(v - Vm);
	    } else if (runMode == RTZ) {
	      rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	    }
	    varOb.setOb(rhov);
	    varOb.setError(configHash.value("aeri_rhov_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 1);

	  }
	  if ((w != -999) and (rho != -999)) {
	    // rho w 1.5 m/s error
	    varOb.setWeight(1., 2);
	    rhow = rho*w;
	    varOb.setOb(rhow);
	    varOb.setError(configHash.value("aeri_rhow_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 2);
	  }
	  if (tempk != -999) {
	    // temperature 1 K error
	    varOb.setWeight(1., 3);
	    varOb.setOb(tempk - tBar);
	    varOb.setError(configHash.value("aeri_tempk_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 3);
	  }
	  if (qv != -999) {
	    // Qv 0.5 g/kg error
	    varOb.setWeight(1., 4);
	    qv = refstate->bhypTransform(qv);
	    varOb.setOb(qv-qBar);
	    varOb.setError(configHash.value("aeri_qv_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 4);
	  }
	  if (rhoa != -999) {
	    // Rho prime .1 kg/m^3 error
	    varOb.setWeight(1., 5);
	    varOb.setOb((rhoa-rhoBar)*100);
	    varOb.setError(configHash.value("aeri_rhoa_error").toFloat());
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 5);
	  }
	  break;
	}
      }

    } // for everything in metData

    // Show a summary of what got tossed out
    
    cout << "Observation problem: " << obsProblem << ", Time problem: " << timeProblem
	 << ", Coordinate problem: " <<  coordProblem << ", Domain problem: " << domainProblem
	 << ", Radius problem: " << radiusProblem << endl;
    
    int newobs = obVector.size() - prevobs;
    // if (metData->size() > 0) {
    if (newobs > 0) {
      cout << "Processed " << newobs << " observations from " << metData->size() << " entries ("
	   << 100.0*(float)newobs/(6*(float)metData->size()) << "%)\n";
    } else {
      cout << "No valid observations in file\n";
    }
    cout << obVector.size() << " total observations." << endl;
  }

  delete metData;

  // Finish reflectivity interpolation
  
  Observation varOb;
  varOb.setTime(configHash.value("ref_time").toInt());
  real pseudow_weight = configHash.value("dbz_pseudow_weight").toFloat();
  real mc_weight = configHash.value("mc_weight").toFloat();
  // Initialize the weights
  for (unsigned int var = 0; var < numVars; ++var) {
    for (unsigned int d = 0; d < numDerivatives; ++d) {
      varOb.setWeight(0.0, var, d);
    }
  }
  real gausspoint = 0.5*sqrt(1./3.);
  for (int iIndex = -1; iIndex < idim; iIndex++) {
    for (int ihalf = 0; ihalf <= 1; ihalf++) {
      for (int imu = -ihalf; imu <= ihalf; imu++) {
	real i = imin + iincr * (iIndex + (gausspoint * imu + 0.5*ihalf));
	for (int jIndex = -1; jIndex < jdim; jIndex++) {
	  for (int jhalf =0; jhalf <= 1; jhalf++) {
	    for (int jmu = -jhalf; jmu <= jhalf; jmu++) {
	      real j = jmin + jincr * (jIndex + (gausspoint * jmu + 0.5*jhalf));
	      real maxrefHeight = -1;
	      for (int kIndex = -1; kIndex < kdim; kIndex++) {
		for (int khalf =0; khalf <= 1; khalf++) {
		  for (int kmu = -khalf; kmu <= khalf; kmu++) {
		    real k = kmin + kincr * (kIndex + (gausspoint * kmu + 0.5*khalf));
		    // On the mish
		    if (ihalf and jhalf and khalf and
			(imu != 0) and (jmu != 0) and (kmu != 0)){
		      int bgI = (iIndex+1)*2 + (imu+1)/2;
		      int bgJ = (jIndex+1)*2 + (jmu+1)/2;
		      int bgK = (kIndex+1)*2 + (kmu+1)/2;
		      int bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
		      if (bgWeights[bIndex] != 0) {
			bgU[bIndex +6] /= bgWeights[bIndex];
		      }
		      if (configHash.value("qr_variable") == "dbz") {
			if (bgU[bIndex +6] > 0) {
			  real dbzavg = 10* log10(bgU[bIndex +6]);
			  bgU[bIndex +6] = (dbzavg+35.)*0.1;
			} else {
			  bgU[bIndex +6] = 0.0;
			}
			if (bgU[bIndex +6] > 3.5) {
			  maxrefHeight = k;
			}
		      }
		    }
		    // On the nodes for mass continuity
		    if ((mc_weight > 0.0) and
			!ihalf and !jhalf and !khalf and
			(i >= imin) and (i <= imax) and
			(j >= jmin) and (j <= jmax) and
			(k >= kmin) and (k <= kmax)) {
		      if (runMode == XYZ) {
			varOb.setCartesianX(i);
			varOb.setCartesianY(j);
			varOb.setWeight(1.0, 0, 1);
			varOb.setWeight(1.0, 1, 2);
			varOb.setWeight(1.0, 2, 3);
		      } else if (runMode == RTZ) {
			if (i > 0) {
			  varOb.setRadius(i);
			  varOb.setTheta(j);
			  real rInverse = 180.0/(i*Pi);
			  varOb.setWeight((1.0/i), 0, 0);
			  varOb.setWeight(1.0, 0, 1);
			  varOb.setWeight(rInverse, 1, 2);
			  varOb.setWeight(1.0, 2, 3);
			}
		      }
		      varOb.setAltitude(k);
		      varOb.setError(mc_weight);
		      varOb.setOb(0.);
		      obVector.push_back(varOb);
		    }
		  }
		}
	      }
	      varOb.setWeight(0.0, 0, 1);
	      varOb.setWeight(0.0, 1, 2);
	      varOb.setWeight(0.0, 2, 3);
	      varOb.setWeight(1., 2);
	      if (runMode == XYZ) {
		varOb.setCartesianX(i);
		varOb.setCartesianY(j);
	      } else if (runMode == RTZ) {
		varOb.setRadius(i);
		varOb.setTheta(j);
	      }
	      varOb.setError(pseudow_weight);
	      varOb.setOb(0.);
	      if (!ihalf and !jhalf){
		// Set an upper boundary condition for W
		if ((maxrefHeight > 0) and (maxrefHeight < kmax)
		    and (pseudow_weight > 0.0)) {
		  varOb.setAltitude(maxrefHeight);
		  obVector.push_back(varOb);
		}

		// Set a lower boundary condition for W
		// Ideally use a terrain map here, but just use Z=0 for now
		if (pseudow_weight > 0.0) {
		  varOb.setAltitude(0.0);
		  varOb.setError(pseudow_weight);
		  obVector.push_back(varOb);
		}
	      }
	      varOb.setWeight(0., 2);
	    }
	  }
	}
      }
    }
  }

  cout << obVector.size() << " total observations including pseudo-obs for W and mass continuity" << endl;

  // Write the Obs to a summary text file
  QString obFilename = dataPath.absoluteFilePath("samurai_Observations.in");
  ofstream obstream(obFilename.toLatin1().data());
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
    real invError = ob.getInverseError();
    if (!invError) {
      cout << "Undefined instrument error specification for " << ob.getType() << "instrument type!\n";
      return false;
    }
    *od++ = invError;
    if (runMode == XYZ) {
      *od++ = ob.getCartesianX();
      *od++ = ob.getCartesianY();
    } else if (runMode == RTZ) {
      *od++ = ob.getRadius();
      *od++ = ob.getTheta();
    }
    *od++ = ob.getAltitude();
    *oi++ = ob.getType();
    *oi++ = ob.getTime();
    for (unsigned int var = 0; var < numVars; var++) {
      for (unsigned int d = 0; d < numDerivatives; ++d) {
	*od++ = ob.getWeight(var, d);
      }
    }
    obstream << endl;
  }

  // Load the observations into a vector
  obs = new real[obVector.size()*(7+numVars*numDerivatives)];
  for (int m=0; m < obVector.size(); m++) {
    int n = m*(7+numVars*numDerivatives);
    Observation ob = obVector.at(m);
    obs[n] = ob.getOb();
    real invError = ob.getInverseError();
    if (!invError) {
      cout << "Undefined instrument error specification for " << ob.getType() << "instrument type!\n";
      return false;
    }
    obs[n+1] = invError;
    if (runMode == XYZ) {
      obs[n+2] = ob.getCartesianX();
      obs[n+3] = ob.getCartesianY();
    } else if (runMode == RTZ) {
      obs[n+2] = ob.getRadius();
      obs[n+3] = ob.getTheta();
    }
    obs[n+4] = ob.getAltitude();
    obs[n+5] = ob.getType();
    obs[n+6] = ob.getTime();
    for (unsigned int var = 0; var < numVars; var++) {
      for (unsigned int d = 0; d < numDerivatives; ++d) {
	int wgt_index = n + (7*(d+1)) + var;
	obs[wgt_index] = ob.getWeight(var, d);
      }
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

bool VarDriver3D::loadPreProcessMetObs()
{

    Observation varOb;
    real wgt[numVars][4];
    real iPos, jPos, kPos, ob, error;
    int type;
    int time;
    cout << "Loading preprocessed observations from samurai_Observations.in" << endl;

    // Open and read the file
    QString obFilename = dataPath.absoluteFilePath("samurai_Observations.in");
    ifstream obstream(obFilename.toLatin1().data());
    while (obstream >> ob >> error >> iPos >> jPos >> kPos >> type >> time
			>> wgt[0][0] >> wgt[0][1] >> wgt[0][2] >> wgt[0][3]
			>> wgt[1][0] >> wgt[1][1] >> wgt[1][2] >> wgt[1][3]
   			>> wgt[2][0] >> wgt[2][1] >> wgt[2][2] >> wgt[2][3]
  			>> wgt[3][0] >> wgt[3][1] >> wgt[3][2] >> wgt[3][3]
   			>> wgt[4][0] >> wgt[4][1] >> wgt[4][2] >> wgt[4][3]
   			>> wgt[5][0] >> wgt[5][1] >> wgt[5][2] >> wgt[5][3]
   			>> wgt[6][0] >> wgt[6][1] >> wgt[6][2] >> wgt[6][3])
    {
        varOb.setOb(ob);
        if (runMode == XYZ) {
            varOb.setCartesianX(iPos);
            varOb.setCartesianY(jPos);
        } else if (runMode == RTZ) {
            varOb.setRadius(iPos);
            varOb.setTheta(jPos);
        }
        varOb.setAltitude(kPos);
        varOb.setType(type);
        varOb.setTime(time);
        varOb.setError(1./error);
        for (unsigned int var = 0; var < numVars; var++) {
            for (unsigned int d = 0; d < numDerivatives; ++d) {
                varOb.setWeight(wgt[var][d],var, d);
            }
        }
        obVector.push_back(varOb);
    }

    // Load the observations into the vector
    obs = new real[obVector.size()*(7+numVars*numDerivatives)];
    for (int m=0; m < obVector.size(); m++) {
        int n = m*(7+numVars*numDerivatives);
        Observation ob = obVector.at(m);
        obs[n] = ob.getOb();
        real invError = ob.getInverseError();
        if (!invError) {
            cout << "Undefined instrument error specification for " << ob.getType() << "instrument type!\n";
            return false;
        }
        obs[n+1] = invError;
        if (runMode == XYZ) {
            obs[n+2] = ob.getCartesianX();
            obs[n+3] = ob.getCartesianY();
        } else if (runMode == RTZ) {
            obs[n+2] = ob.getRadius();
            obs[n+3] = ob.getTheta();
        }
        obs[n+4] = ob.getAltitude();
        obs[n+5] = ob.getType();
        obs[n+6] = ob.getTime();
        for (unsigned int var = 0; var < numVars; var++) {
            for (unsigned int d = 0; d < numDerivatives; ++d) {
                int wgt_index = n + (7*(d+1)) + var;
                obs[wgt_index] = ob.getWeight(var, d);
            }
        }
    }

    return true;
}

#if 0
// Load the background estimates from the given array arguments

int VarDriver3D::loadBackgroundObs(int nx, int ny, int nsigma,
				   char *ctdg, int delta, int iter, // time elements
				   float *sigmas,
				   float *latitude, float *longitude,
				   float *u1,
				   float *v1,
				   float *w1,
				   float *th1,
				   float *p1)
{
#if 0
  // TODO debug
    for(int i = 0; i < nx; i++)
      for(int j = 0; j < ny; j++)
	for(int k = 0; k < nsigma; k++) {
	  float f = *(u1 + i + nx * (j + ny * k));
	  std::cout << "i: " << i
		    << ", j: " << j
		    << ", k: " << k
		    << ", val: " << f
		    << std::endl;
	}
#endif
    
  return loadBackgroundObs();  
}

#endif

// Load the background estimates from a file

int VarDriver3D::loadBackgroundObs(const char *background_fname)
{
  bkgdAdapter = new BkgdStream(background_fname);
  return loadBackgroundObs();
}

// Background Observations can come either from a file, or from passed arguments to the run() function.

int VarDriver3D::loadBackgroundObs()
{
  if(bkgdAdapter == NULL) {
    std::cerr << "Error: VarDriver3D::loadBackgroundObs() called without adapter initialization." << std::endl;
    exit(1);
  }
    
  // Turn Debug on if there are problems with the vertical spline interpolation,
  // Eventually this should be replaced with the internal spline code
  // SplineD::Debug(1);

  // Geographic functions
  //GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();

  real referenceLon = configHash.value("ref_lon").toFloat();

  QVector<real> logheights, uBG, vBG, wBG, tBG, qBG, rBG, zBG;
  //SplineD* bgSpline;
    
  int time;
  QString bgTimestring, tcstart, tcend;
  real lat, lon, alt, u, v, w, t, qv, rhoa, qr;
  real bgX, bgY, bgRadius, bgTheta;
  real bgZ = -32768.;
  // backgroundroi is in km, ROI is gridpoints
  real iROI = configHash.value("i_background_roi").toFloat() / iincr;
  real jROI = configHash.value("j_background_roi").toFloat() / jincr;
  QString interp_mode = configHash.value("bg_interpolation");
  real maxGridDist = 3.0;
  if (interp_mode == "Cressman") {
    maxGridDist = 1.0;
  }
  real Rsquare = (iincr*iROI)*(iincr*iROI) + (jincr*jROI)*(jincr*jROI);

  if (frameVector.size() == 0) {
    cout << "Frame Vector is not initialized."  << endl;
    return -1;
  }
    
  cout << "Loading background onto Gaussian mish with " << iROI
       << " grid length radius of influence in i direction" << endl;
  cout << "and " << jROI << " grid length radius of influence in j direction" << endl;

  QDateTime startTime = frameVector.front().getTime();
  QDateTime endTime = frameVector.back().getTime();

  cout << "Start time: " << startTime.toString("yyyy-MM-dd-HH:mm:ss").toLatin1().data() << endl;
  cout << "End  time:  " << endTime.toString("yyyy-MM-dd-HH:mm:ss").toLatin1().data() << endl;    

  int debug = 0;
    
  while( bkgdAdapter->next(time, lat, lon, alt, u, v, w, t, qv, rhoa, qr) )
    {
      // Process the metObs into Observations

      // Make sure the bg is within the time limits
      QDateTime bgTime;
      bgTime.setTime_t(time);
      bgTime.setTimeSpec(Qt::UTC);
      if(debug++ < 1000)
	cout << "BgObs time:       " << bgTime.toString("yyyy-MM-dd-HH:mm:ss").toLatin1().data()    
	     << " time_t(): " << bgTime.toTime_t()
	     << endl;
      bgTimestring = bgTime.toString(Qt::ISODate);
      tcstart = startTime.toString(Qt::ISODate);
      tcend = endTime.toString(Qt::ISODate);
      if ((bgTime < startTime) or (bgTime > endTime))
	continue;
      int tci = startTime.secsTo(bgTime);
      if ((tci < 0) or (tci > (int)frameVector.size())) {
	cout << "Time problem with observation " << tci << "secs more than center entries" << endl;
	continue;
      }

      real Um = frameVector[tci].getUmean();
      real Vm = frameVector[tci].getVmean();

      // Get the X, Y & Z
      real tcX, tcY, metX, metY;
      projection.Forward(referenceLon, frameVector[tci].getLat() , frameVector[tci].getLon(),
			 tcX, tcY);
      projection.Forward(referenceLon, lat, lon , metX, metY);
      bgX = (metX - tcX) / 1000.;
      bgY = (metY - tcY) / 1000.;

      real heightm = alt;

      bgZ = heightm/1000.;
      bgRadius = sqrt(bgX*bgX + bgY*bgY);
      bgTheta = 180.0 * atan2(bgY, bgX) / Pi;
      if (configHash.value("allow_negative_angles") != "true") {
	if (bgTheta < 0) bgTheta += 360.0;
      }

      // Make sure the ob is in the Interpolation domain
	
      if (runMode == XYZ) {
	if ((bgX < (imin-iincr-(iROI*iincr*maxGridDist))) or (bgX > (imax+iincr+(iROI*iincr*maxGridDist))) or
	    (bgY < (jmin-jincr-(jROI*jincr*maxGridDist))) or (bgY > (jmax+jincr+(jROI*jincr*maxGridDist)))
	    or (bgZ < kmin)) //Allow for higher values for interpolation purposes
	  continue;
      } else if (runMode == RTZ) {
	if ((bgRadius < (imin-iincr-(iROI*iincr*maxGridDist))) or (bgRadius > (imax+iincr+(iROI*iincr*maxGridDist))) or
	    (bgTheta < jmin-jincr-(jROI*jincr*maxGridDist)) or (bgTheta > jmax+jincr+(jROI*jincr*maxGridDist)) or
	    (bgZ < kmin)) //Exceeding the Theta domain only makes sense for sectors
	  continue;
	real cylUm = (Um*bgX + Vm*bgY)/bgRadius;
	real cylVm = (-Um*bgY + Vm*bgX)/bgRadius;
	Um = cylUm;
	Vm = cylVm;
      }

      // Reference states
      real rhoBar = refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm);
      real qBar = refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm);
      real tBar = refstate->getReferenceVariable(ReferenceVariable::tempref, heightm);

      real rho = rhoa + rhoa*qv/1000.;
      real rhou = rho*(u - Um);
      real rhov = rho*(v - Vm);
      real rhow = rho*w;
      real tprime = t - tBar;
      qv = refstate->bhypTransform(qv);
      real qvprime = qv-qBar;
      real rhoprime = (rhoa-rhoBar)*100;
      real logZ = log(bgZ);
      if (configHash.value("qr_variable") == "qr") {
	qr = refstate->bhypTransform(qr);
      }
      // We assume here that the background precipitation field is always zero
      // real qr = 0.;

      if (runMode == XYZ) {
	bgIn << bgX << bgY << logZ << time << rhou << rhov << rhow << tprime << qvprime << rhoprime << qr ;
      } else if (runMode == RTZ) {
	bgIn << bgRadius << bgTheta << logZ << time << rhou << rhov << rhow << tprime << qvprime << rhoprime << qr ;
      }
      if (logheights.size() == 0) {
	// First column
	logheights.push_back(logZ);
	uBG.push_back(rhou);
	vBG.push_back(rhov);
	wBG.push_back(rhow);
	tBG.push_back(tprime);
	qBG.push_back(qvprime);
	rBG.push_back(rhoprime);
	zBG.push_back(qr);
      } else if (logZ > logheights.back()) {
	// Same column
	logheights.push_back(logZ);
	uBG.push_back(rhou);
	vBG.push_back(rhov);
	wBG.push_back(rhow);
	tBG.push_back(tprime);
	qBG.push_back(qvprime);
	rBG.push_back(rhoprime);
	zBG.push_back(qr);
      } else {
	// Solve for the spline
	if (logheights.size() == 1) {
	  cerr << "Error at " << lat << ", " << lon << "\n";
	  cerr << "Only one level found in background spline setup. "
	       << "Please check Background.in to ensure sorting by Z coordinate and re-run." << endl;
	  return -1;
	}

	// Exponential interpolation in horizontal, b-Spline interpolation on log height in vertical
#pragma omp parallel for
	for (int ki = -1; ki < (kdim); ki++) {
	  SplineD* bgSpline = new SplineD(&logheights.front(), logheights.size(),
					  uBG.data(), 0, SplineBase::BC_ZERO_SECOND);
	  if (!bgSpline->ok()) {
	    cerr << "bgSpline setup failed." << endl;
	    continue; //return -1;
	  }
	  for (int kmu = -1; kmu <= 1; kmu += 2) {
	    real kPos = kmin + kincr * (ki + (0.5*sqrt(1./3.) * kmu + 0.5));
	    if (kPos < 0) kPos = 0.001;
	    real logzPos = log(kPos);
	    if (logzPos < logheights[0]) logzPos = logheights[0];
	    //if (fabs(kPos-obZ) > kincr*ROI*2.) continue;
	    for (int ii = -1; ii < (idim); ii++) {
	      for (int imu = -1; imu <= 1; imu += 2) {
		real iPos = imin + iincr * (ii + (0.5*sqrt(1./3.) * imu + 0.5));
		if (runMode == XYZ) {
		  if (fabs(iPos-bgX) > iincr*iROI*maxGridDist) continue;
		} else if (runMode == RTZ) {
		  if (fabs(iPos-bgRadius) > iincr*iROI*maxGridDist) continue;
		}
		for (int ji = -1; ji < (jdim); ji++) {
		  for (int jmu = -1; jmu <= 1; jmu += 2) {
		    real jPos = jmin + jincr * (ji + (0.5*sqrt(1./3.) * jmu + 0.5));
		    real rSquare = 0.0;
		    if (runMode == XYZ) {
		      if (fabs(jPos-bgY) > jincr*jROI*maxGridDist) continue;
		      rSquare = (bgX-iPos)*(bgX-iPos) + (bgY-jPos)*(bgY-jPos);
		    } else if (runMode == RTZ) {
		      real dTheta = fabs(jPos-bgTheta);
		      if (dTheta > 360.) dTheta -= 360.;
		      if (dTheta > jincr*jROI*2.) continue;
		      rSquare = (bgRadius-iPos)*(bgRadius-iPos) + (dTheta)*(dTheta);
		    }
		    // Add one extra index to account for buffer zone in analysis
		    int bgI = (ii+1)*2 + (imu+1)/2;
		    int bgJ = (ji+1)*2 + (jmu+1)/2;
		    int bgK = (ki+1)*2 + (kmu+1)/2;
		    int bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
		    if (rSquare < Rsquare*maxGridDist) {
		      real weight = exp(-2.302585092994045*rSquare/Rsquare);
		      if (interp_mode == "Cressman") {
			weight = (Rsquare - rSquare)/(Rsquare + rSquare);
		      }
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
			bgSpline->solve(zBG.data());
			bgU[bIndex +6] += weight*(bgSpline->evaluate(logzPos));
			bgWeights[bIndex] += weight;
		      } else {
			// Below the spline interpolation
			bgU[bIndex] += weight*uBG.front();
			bgU[bIndex +1] += weight*vBG.front();
			bgU[bIndex +2] += weight*wBG.front();
			bgU[bIndex +3] += weight*tBG.front();
			bgU[bIndex +4] += weight*qBG.front();
			bgU[bIndex +5] += weight*rBG.front();
			bgU[bIndex +6] += weight*zBG.front();
			bgWeights[bIndex] += weight;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  delete bgSpline;
	}

	logheights.clear();
	uBG.clear();
	vBG.clear();
	wBG.clear();
	tBG.clear();
	qBG.clear();
	rBG.clear();
	zBG.clear();

	logheights.push_back(log(bgZ));
	uBG.push_back(rhou);
	vBG.push_back(rhov);
	wBG.push_back(rhow);
	tBG.push_back(tprime);
	qBG.push_back(qvprime);
	rBG.push_back(rhoprime);
	zBG.push_back(qr);
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
  if (logheights.size() == 1) {
    cerr << "Only one level found in background spline setup. Please check Background.in to ensure sorting by Z coordinate and re-run." << endl;
    return -1;
  }

  // Exponential interpolation in horizontal, b-Spline interpolation on log height in vertical
#pragma omp parallel for
  for (int ki = -1; ki < (kdim); ki++) {
    SplineD* bgSpline = new SplineD(&logheights.front(), logheights.size(), uBG.data(), 2, SplineBase::BC_ZERO_SECOND);
    if (!bgSpline->ok())
      {
	cerr << "bgSpline setup failed." << endl;
	continue; //return -1;
      }
    for (int kmu = -1; kmu <= 1; kmu += 2) {
      real kPos = kmin + kincr * (ki + (0.5*sqrt(1./3.) * kmu + 0.5));
      if (kPos < 0) kPos = 0.001;
      real logzPos = log(kPos);
      if (logzPos < logheights[0]) logzPos = logheights[0];
      //if (fabs(kPos-obZ) > kincr*ROI*2.) continue;
      for (int ii = -1; ii < (idim); ii++) {
	for (int imu = -1; imu <= 1; imu += 2) {
	  real iPos = imin + iincr * (ii + (0.5*sqrt(1./3.) * imu + 0.5));
	  if (runMode == XYZ) {
	    if (fabs(iPos-bgX) > iincr*iROI*maxGridDist)
	      continue;
	  } else if (runMode == RTZ) {
	    if (fabs(iPos-bgRadius) > iincr*iROI*maxGridDist)
	      continue;
	  }
	  for (int ji = -1; ji < (jdim); ji++) {
	    for (int jmu = -1; jmu <= 1; jmu += 2) {
	      real jPos = jmin + jincr * (ji + (0.5*sqrt(1./3.) * jmu + 0.5));
	      real rSquare = 0.0;
	      if (runMode == XYZ) {
		if (fabs(jPos-bgY) > jincr*jROI*maxGridDist)
		  continue;
		rSquare = (bgX-iPos)*(bgX-iPos) + (bgY-jPos)*(bgY-jPos);
	      } else if (runMode == RTZ) {
		real dTheta = fabs(jPos-bgTheta);
		if (dTheta > 360.) dTheta -= 360.;
		if (dTheta > jincr*jROI*2.)
		  continue;
		rSquare = (bgRadius-iPos)*(bgRadius-iPos) + dTheta*dTheta;
	      }
	      // Add one extra index to account for buffer zone in analysis
	      int bgI = (ii+1)*2 + (imu+1)/2;
	      int bgJ = (ji+1)*2 + (jmu+1)/2;
	      int bgK = (ki+1)*2 + (kmu+1)/2;
	      int bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
	      if (rSquare < Rsquare*maxGridDist) {
		real weight = exp(-2.302585092994045*rSquare/Rsquare);
		if (interp_mode == "Cressman") {
		  weight = (Rsquare - rSquare)/(Rsquare + rSquare);
		}
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
		  bgSpline->solve(zBG.data());
		  bgU[bIndex +6] += weight*(bgSpline->evaluate(logzPos));
		  bgWeights[bIndex] += weight;
		} else {
		  // Below the spline interpolation
		  bgU[bIndex] += weight*uBG.front();
		  bgU[bIndex +1] += weight*vBG.front();
		  bgU[bIndex +2] += weight*wBG.front();
		  bgU[bIndex +3] += weight*tBG.front();
		  bgU[bIndex +4] += weight*qBG.front();
		  bgU[bIndex +5] += weight*rBG.front();
		  bgU[bIndex +6] += weight*zBG.front();
		  bgWeights[bIndex] += weight;
		}
	      }
	    }
	  }
	}
      }
    }
    delete bgSpline;
  }

  logheights.clear();
  uBG.clear();
  vBG.clear();
  wBG.clear();
  tBG.clear();
  qBG.clear();
  rBG.clear();
  zBG.clear();

  int numbgObs = bgIn.size()*7/11;
  QVector<int> emptybg;
    
  if (numbgObs > 0) {
    // Check interpolation
    for (int ki = -1; ki < (kdim); ki++) {
      for (int kmu = -1; kmu <= 1; kmu += 2) {
	real kPos = kmin + kincr * (ki + (0.5*sqrt(1./3.) * kmu + 0.5));
	for (int ii = -1; ii < (idim); ii++) {
	  for (int imu = -1; imu <= 1; imu += 2) {
	    real iPos = imin + iincr * (ii + (0.5*sqrt(1./3.) * imu + 0.5));
	    for (int ji = -1; ji < (jdim); ji++) {
	      for (int jmu = -1; jmu <= 1; jmu += 2) {
		real jPos = jmin + jincr * (ji + (0.5*sqrt(1./3.) * jmu + 0.5));
		int bgI = (ii+1)*2 + (imu+1)/2;
		int bgJ = (ji+1)*2 + (jmu+1)/2;
		int bgK = (ki+1)*2 + (kmu+1)/2;
		int bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
		for (unsigned int var = 0; var < numVars; var++) {
		  if (bgWeights[bIndex] != 0) {
		    bgU[bIndex +var] /= bgWeights[bIndex];
		  } else {
		    emptybg.push_back(bIndex);
		    if (emptybg.size() < 15) {
		      cout << "Empty background mish for variable " << var << " at " << iPos << ", " << jPos << ", " << kPos << endl;
		    } else if (emptybg.size() == 15) {
		      cout << "Too many empty mish points, will no longer report.\n";
		    }
		  }
		}
		if (configHash.value("qr_variable") == "dbz") {
		  if (bgU[bIndex +6] > 0) {
		    real dbzavg = 10* log10(bgU[bIndex +6]);
		    bgU[bIndex +6] = (dbzavg+35.)*0.1;
		  } else {
		    bgU[bIndex +6] = 0.0;
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    // Attempt to fill holes everywhere
    if (emptybg.size() > 0) {
      int neighbors[6];
      real avg[numVars];
      // int uStateSize = 8*(idim+1)*(jdim+1)*(kdim+1)*(numVars);
      for (int i=0; i < emptybg.size(); i++) {
	for (unsigned int var = 0; var < numVars; var++) avg[var] = 0.0;
	// Check the neighbors
	int bIndex = emptybg[i];
	neighbors[0] = bIndex + numVars*(idim+1)*2*(jdim+1)*2;
	neighbors[1] = bIndex - numVars*(idim+1)*2*(jdim+1)*2;
	neighbors[2] = bIndex + numVars*(idim+1)*2;
	neighbors[3] = bIndex - numVars*(idim+1)*2;
	neighbors[4] = bIndex + 1;
	neighbors[5] = bIndex - 1;
	real count = 0.0;
	for (int j = 0; j < 6; j++) {
	  if ((neighbors[j] >= 0) and (neighbors[j] < uStateSize) and (bgWeights[neighbors[j]])) {
	    // A good neighbor
	    for (unsigned int var = 0; var < numVars; var++) {
	      avg[var] += bgU[neighbors[j] + var];
	      count++;
	    }
	  }
	}
	// Need at least 3 neighbors
	if (count > 2) {
	  for (unsigned int var = 0; var < numVars; var++) {
	    bgU[bIndex + var] = avg[var]/count;
	  }
	} else {
	  if (configHash.value("allow_background_missing_values") != "true") {
	    cout << "Too large a hole in the background field!\n";
	    cout << "Please check your background file or ROI values.\n";
	    cout << "If you want to allow missing (zero) values, add the following line to <options>:\n";
	    cout << "<allow_background_missing_values>true</allow_background_missing_values>\n";
	    exit(-1);
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

bool VarDriver3D::adjustBackground()
{
  /* Set the minimum filter length to the background resolution, not the analysis resolution
     to avoid artifacts when running interpolating to small mesoscale grids */

  // Load the observations into a vector
  int numbgObs = bgIn.size()*7/11;
  if (configHash.value("mc_weight").toFloat() > 0) {
    numbgObs += idim*jdim*kdim;
  }
  bgObs = new real[numbgObs*(7+numVars*numDerivatives)];
  for (unsigned int m=0; m < numbgObs*(7+numVars*numDerivatives); m++) bgObs[m] = 0.;

  int p = 0;
  real obX, obY, obRadius, obTheta;
  obX = obY = obRadius = obTheta = -32768.;
  for (int m=0; m < bgIn.size(); m+=11) {
    if (runMode == XYZ) {
      obX = bgIn[m];
      obY = bgIn[m+1];
    } else if (runMode == RTZ) {
      obRadius = bgIn[m];
      obTheta = bgIn[m+1];
    }
    real obZ = exp(bgIn[m+2]);
    real obTime = bgIn[m+3];
    // Make sure the ob is in the domain
    if (runMode == XYZ) {
      if ((obX < imin) or (obX > imax) or
	  (obY < jmin) or (obY > jmax) or
	  (obZ < kmin) or (obZ > kmax)) {
	numbgObs -= 7;
	continue;
      }
    } else if (runMode == RTZ) {
      if ((obRadius < imin) or (obRadius > imax) or
	  (obTheta < jmin) or (obTheta > jmax) or
	  (obZ < kmin) or (obZ > kmax)) {
	numbgObs -= 7;
	continue;
      }
    }

    for (unsigned int n = 0; n < numVars; n++) {
      bgObs[p] = bgIn[m+4+n];
      if ((n == 6) and (configHash.value("qr_variable") == "dbz")) {
	// Convert to dBZ control variable
	real dbzavg = 10* log10(bgIn[m+4+n]);
	bgObs[p] = (dbzavg+35.)*0.1;
      }
      // Default error of background = 0.1
      if (configHash.value("bg_obs_error").isEmpty() or
	  (configHash.value("bg_obs_error").toFloat() <= 0.0)) {
	bgObs[p+1] = 100.;
      } else {
	bgObs[p+1] = 1.0/configHash.value("bg_obs_error").toFloat();
      }
      if (runMode == XYZ) {
	bgObs[p+2] = obX;
	bgObs[p+3] = obY;
      } else if (runMode == RTZ) {
	bgObs[p+2] = obRadius;
	bgObs[p+3] = obTheta;
      }
      bgObs[p+4] = obZ;
      // Null type
      bgObs[p+5] = -1;
      bgObs[p+6] = obTime;
      bgObs[p+7+n] = 1.;
      p += (7+numVars*numDerivatives);
    }
  }

  // Add mass continuity constraint
  if (configHash.value("mc_weight").toFloat() > 0) {
    for (int iIndex = 0; iIndex < idim; iIndex++) {
      real i = imin + iincr * iIndex;
      if (i > ((idim-1)*iincr + imin)) continue;
      for (int jIndex = 0; jIndex < jdim; jIndex++) {
	real j = jmin + jincr * jIndex;
	if (j > ((jdim-1)*jincr + jmin)) continue;
	for (int kIndex = 0; kIndex < kdim; kIndex++) {
	  real k = kmin + kincr * kIndex;
	  if (k > ((kdim-1)*kincr + kmin)) continue;
	  bgObs[p] = 0.0;
	  bgObs[p+1] = configHash.value("mc_weight").toFloat();
	  bgObs[p+2] = i;
	  bgObs[p+3] = j;
	  bgObs[p+4] = k;
	  // Null type
	  bgObs[p+5] = -1;
	  bgObs[p+6] = configHash.value("ref_time").toInt();
	  if (runMode == XYZ) {
	    bgObs[p+(7*(1+1))] = 1.0;
	    bgObs[p+(7*(2+1))+1] = 1.0;
	    bgObs[p+(7*(3+1))+2] = 1.0;
	  } else if (runMode == RTZ) {
	    if (i > 0) {
	      real rInverse = 180.0/(i*Pi);
	      bgObs[p+7] = 1.0/i;
	      bgObs[p+(7*(1+1))] = 1.0;
	      bgObs[p+(7*(2+1))+1] = rInverse;
	      bgObs[p+(7*(3+1))+2] = 1.0;
	    }
	  }
	  p += (7+numVars*numDerivatives);
	}
      }
    }
  }

  // Store and set the background errors
  QString bgError[7];
  bgError[0] = configHash.value("bg_rhou_error");
  bgError[1] = configHash.value("bg_rhov_error");
  bgError[2] = configHash.value("bg_rhow_error");
  bgError[3] = configHash.value("bg_tempk_error");
  bgError[4] = configHash.value("bg_qv_error");
  bgError[5] = configHash.value("bg_rhoa_error");
  bgError[6] = configHash.value("bg_qr_error");

  QString bg_interpolation_error = "1.0";
  if (!configHash.value("bg_interpolation_error").isEmpty()) {
    bg_interpolation_error = configHash.value("bg_interpolation_error");
    cout << "Setting background interpolation error to " << bg_interpolation_error.toStdString() << "\n";
  } else {
    cout << "Using default background interpolation error of 1.0\n";
  }
  configHash["bg_rhou_error"] = bg_interpolation_error;
  configHash["bg_rhov_error"] = bg_interpolation_error;
  configHash["bg_rhow_error"] = bg_interpolation_error;
  configHash["bg_tempk_error"] = bg_interpolation_error;
  configHash["bg_qv_error"] = bg_interpolation_error;
  configHash["bg_rhoa_error"] = bg_interpolation_error;
  configHash["bg_qr_error"] = bg_interpolation_error;
  configHash["save_mish"] = "true";
  // Adjust the background field to the spline mish
  if (runMode == XYZ) {
    if (configHash.value("output_pressure_increment").toFloat() > 0) {
      bgCost3D = new CostFunctionXYP(projection, numbgObs, bStateSize);
    } else {
      bgCost3D = new CostFunctionXYZ(projection, numbgObs, bStateSize);
    }
  } else if (runMode == RTZ) {
    bgCost3D = new CostFunctionRTZ(projection, numbgObs, bStateSize);
  }
  bgCost3D->initialize(&configHash, bgU, bgObs, refstate);
  /* Set the iteration to zero -- this will prevent writing the background file until after the adjustment
     which is presumably what you want most of the time. Otherwise, you would not be here */
  int bgIter = 1;
  bgCost3D->initState(bgIter);
  bgCost3D->minimize();

  // Increment the variables
  bgCost3D->updateBG();
  bgCost3D->finalize();

  delete bgCost3D;
  delete[] bgObs;

  // Reset the background errors
  configHash["bg_rhou_error"] = bgError[0];
  configHash["bg_rhov_error"] = bgError[1];
  configHash["bg_rhow_error"] = bgError[2];
  configHash["bg_tempk_error"] = bgError[3];
  configHash["bg_qv_error"] = bgError[4];
  configHash["bg_rhoa_error"] = bgError[5];
  configHash["bg_qr_error"] = bgError[6];
  configHash["save_mish"] = "false";

  // Convert the dBZ back to Z for further processing
  if (configHash.value("qr_variable") == "dbz") {
    for (int ki = -1; ki < (kdim); ki++) {
      for (int kmu = -1; kmu <= 1; kmu += 2) {
	for (int ii = -1; ii < (idim); ii++) {
	  for (int imu = -1; imu <= 1; imu += 2) {
	    for (int ji = -1; ji < (jdim); ji++) {
	      for (int jmu = -1; jmu <= 1; jmu += 2) {
		int bgI = (ii+1)*2 + (imu+1)/2;
		int bgJ = (ji+1)*2 + (jmu+1)/2;
		int bgK = (ki+1)*2 + (kmu+1)/2;
		int bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
		real dbZ = bgU[bIndex +6]*10. - 35.;
		real ZZ = pow(10.0,(dbZ*0.1));
		bgU[bIndex +6] = ZZ;
	      }
	    }
	  }
	}
      }
    }
  }

  return true;
}


/* Any updates needed for additional analysis iterations go here */

void VarDriver3D::updateAnalysisParams(const int& iteration)
{
    QString iter;
    iter.setNum(iteration);

    QString key = "bg_rhou_error_" + iter;
    QString val = configHash.value(key);
    configHash.insert("bg_rhou_error", val);

    key = "bg_rhov_error_" + iter;
    val = configHash.value(key);
    configHash.insert("bg_rhov_error", val);

    key = "bg_rhow_error_" + iter;
    val = configHash.value(key);
    configHash.insert("bg_rhow_error", val);

    key = "bg_tempk_error_" + iter;
    val = configHash.value(key);
    configHash.insert("bg_tempk_error", val);

    key = "bg_qv_error_" + iter;
    val = configHash.value(key);
    configHash.insert("bg_qv_error", val);

    key = "bg_rhoa_error_" + iter;
    val = configHash.value(key);
    configHash.insert("bg_rhoa_error", val);

    key = "bg_qr_error_" + iter;
    val = configHash.value(key);
    configHash.insert("bg_qr_error", val);

    key = "mc_weight_" + iter;
    val = configHash.value(key);
    configHash.insert("mc_weight", val);

    key = "i_filter_length_" + iter;
    val = configHash.value(key);
    configHash.insert("i_filter_length", val);

    key = "j_filter_length_" + iter;
    val = configHash.value(key);
    configHash.insert("j_filter_length", val);

    key = "k_filter_length_" + iter;
    val = configHash.value(key);
    configHash.insert("k_filter_length", val);

    key = "i_spline_cutoff_" + iter;
    val = configHash.value(key);
    configHash.insert("i_spline_cutoff", val);

    key = "j_spline_cutoff_" + iter;
    val = configHash.value(key);
    configHash.insert("j_spline_cutoff", val);

    key = "k_spline_cutoff_" + iter;
    val = configHash.value(key);
    configHash.insert("k_spline_cutoff", val);

}

/* This routine validates that all required parameters are present
 It currently does not check the validity of a particular parameter, just that it exists */

bool VarDriver3D::validateConfig()
{
    // Validate the hash -- multiple passes are not validated currently
    QStringList configKeys;
    configKeys << "mc_weight" << "i_min" << "i_max" << "i_incr" <<
    "j_min" << "j_max" << "j_incr" <<
    "k_min" << "k_max" << "k_incr" <<
    "i_filter_length" << "j_filter_length" << "k_filter_length" <<
    "i_spline_cutoff" << "j_spline_cutoff" << "k_spline_cutoff" <<
    "i_rhou_bcL" << "i_rhou_bcR" << "j_rhou_bcL" << "j_rhou_bcR" << "k_rhou_bcL" << "k_rhou_bcR" <<
    "i_rhov_bcL" << "i_rhov_bcR" << "j_rhov_bcL" << "j_rhov_bcR" << "k_rhov_bcL" << "k_rhov_bcR" <<
    "i_rhow_bcL" << "i_rhow_bcR" << "j_rhow_bcL" << "j_rhow_bcR" << "k_rhow_bcL" << "k_rhow_bcR" <<
    "i_tempk_bcL" << "i_tempk_bcR" << "j_tempk_bcL" << "j_tempk_bcR" << "k_tempk_bcL" << "k_tempk_bcR" <<
    "i_qv_bcL" << "i_qv_bcR" << "j_qv_bcL" << "j_qv_bcR" << "k_qv_bcL" << "k_qv_bcR" <<
    "i_rhoa_bcL" << "i_rhoa_bcR" << "j_rhoa_bcL" << "j_rhoa_bcR" << "k_rhoa_bcL" << "k_rhoa_bcR" <<
    "i_qr_bcL" << "i_qr_bcR" << "j_qr_bcL" << "j_qr_bcR" << "k_qr_bcL" << "k_qr_bcR" <<
	"data_directory" << "output_directory";

    for (int i = 0; i < configKeys.count(); i++) {
        if (!configHash.contains(configKeys.at(i))) {
            cout <<	"No configuration found for <" << configKeys.at(i).toStdString() << "> aborting..." << endl;
            return false;
        }
    }
    return true;
}

bool VarDriver3D::loadBackgroundCoeffs()
{
	// Load a set of coefficients directly for the same grid
	cout << "Loading previous coefficients from samurai_Coefficients.in" << endl;

	QFile bgFile(dataPath.filePath("samurai_Coefficients.in"));
	if (!bgFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	int numbgCoeffs = 0;
	QTextStream in(&bgFile);
	while (!in.atEnd()) {
		QString line = in.readLine();
		if (line.contains("Variable")) {
			continue;
		} else {
			QStringList lineparts = line.split(QRegExp("\\s+"));
			int var = lineparts[0].toInt();
			int iIndex = lineparts[1].toInt();
			int jIndex = lineparts[2].toInt();
			int kIndex = lineparts[3].toInt();
			int bIndex = numVars*(idim+2)*(jdim+2)*kIndex + numVars*(idim+2)*jIndex +numVars*iIndex + var;
			bgU[bIndex] = lineparts[5].toFloat();
			numbgCoeffs++;
		}
	}
	cout << numbgCoeffs << " background coeffients loaded" << endl;
	if (numbgCoeffs != bStateSize) {
	  cout << "Error loading background coefficients" << endl;
	  return false;
	}
	return true;
}

// fill up variables from arguments and call the common validateGrid()

bool VarDriver3D::validateRunGrid(float n_x, float n_y, float n_z,
				  float i_min, float i_max, float i_incr,
				  float j_min, float j_max, float j_incr,
				  float *sigmas)
{
  imin = i_min;
  imax = i_max;
  jmin = j_min;
  jmax = j_max;

  iincr = i_incr;
  jincr = j_incr;

  // Do we need these params, or can we compute them from the _min, _max, _incr?
  
  idim = n_x;
  jdim = n_y;
  kdim = n_z;

  // Compute k dimension
  kmin  = sigmas[0];
  kmax  = sigmas[ (int) n_z - 1];
  kincr = (kmax - kmin) / n_z;
  
  return validateGrid();
}

//

void VarDriver3D::fillRunCenters(char *cdtg, int delta, int iter, float lat, float lon)
{
  // get rid of previous centers
  clearCenters();

  // Compute ref_time
  QDateTime cDateTime;
  QString tString;
  cDateTime.setTimeSpec(Qt::UTC);
  cDateTime = QDateTime::fromString(cdtg, "yyyyMMddHH");
  cDateTime = cDateTime.addSecs(delta * iter);

  configHash.insert("ref_time", tString.setNum(cDateTime.toTime_t()));
  configHash.insert("ref_lat",  tString.setNum(lat));
  configHash.insert("ref_lat",  tString.setNum(lon));  

  // create 3 "centers" centered around the ref time, at 2 second intervals

  float zero = 0.0;	// Framecenter constructor needs a reference... why?
  for(int delta = -2; delta <= 2; delta += 1) {
    QDateTime tTime = cDateTime.addSecs(delta);
    frameVector.push_back(FrameCenter(tTime, lat, lon, zero, zero));
  }    
}
		      
// fill up variables from the config hash  and call the common validateGrid()

bool VarDriver3D::validateFixedGrid()
{
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
  
  return validateGrid();
}

bool VarDriver3D::validateGrid()
{
  // The recursive filter uses a fourth order stencil to spread the observations,
  // so less than 4 gridpoints will cause a memory fault
  
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

  return true;
}

bool VarDriver3D::loadMetObs()
{
  // Load the BG into an empty vector
  
  uStateSize = 8*(idim+1)*(jdim+1)*(kdim+1)*(numVars);
  bStateSize = (idim+2)*(jdim+2)*(kdim+2)*numVars;
  
  cout << "Physical (mish) State size = " << uStateSize << "\n";
  cout << "Nodal State size = " << bStateSize << ", Grid dimensions:\n";

  bgU = new real[uStateSize];
  bgWeights = new real[uStateSize];
  for (int i=0; i < uStateSize; i++) {
    bgU[i] = 0.0;
    bgWeights[i] = 0.0;
  }

  // Read in the meteorological observations, process them into weights and positions
  // Either preprocess from raw observations or load an already processed Observations.in file

  QString preprocess = configHash.value("preprocess_obs");
  if (preprocess == "true") {
    if (! preProcessMetObs()) {
      cout << "Error pre-processing observations\n";
      return false;
    }
  } else {
    if (!loadPreProcessMetObs()) {
      cout << "Error loading observations\n";
      return false;
    }
  }

  // We are done with the bgWeights, so free up that memory
  delete[] bgWeights;

  if (obVector.size() == 0) {
    // No observations so quit
    cout << "No observations loaded, unable to perform analysis.\n";
    return false;
  } else {
    cout << "Number of New Observations: " << obVector.size() << endl;
  }
  return true;
}

bool VarDriver3D::initObCost3D()
{
  if (runMode == XYZ) {
    if (configHash.value("output_pressure_increment").toFloat() > 0) {
      obCost3D = new CostFunctionXYP(projection, obVector.size(), bStateSize);
    } else if (configHash.value("output_COAMPS") == "true") {
      obCost3D = new CostFunctionCOAMPS(projection, obVector.size(), bStateSize);
    } else {
      obCost3D = new CostFunctionXYZ(projection, obVector.size(), bStateSize);
    }
  } else if (runMode == RTZ) {
    obCost3D = new CostFunctionRTZ(projection, obVector.size(), bStateSize);
  }
  
  obCost3D->initialize(&configHash, bgU, obs, refstate);
  return true;
}
