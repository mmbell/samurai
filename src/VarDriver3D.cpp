/*
 *  VarDriver3D.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include <iterator>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <regex>
#include <set>

#include <iomanip>
// #include <netcdfcpp.h>
#include <Ncxx/Nc3File.hh>

#include "VarDriver3D.h"
#include "Dorade.h"
#include "RecursiveFilter.h"
#include "samurai.h"
#include "timers.h"
#include "BkgdObsLoaders.h"
#include "LineSplit.h"
#include "FileList.h"
#include "timing/gptl.h"

// Constructor
VarDriver3D::VarDriver3D() : VarDriver()
{
  numVars = 7;
  numDerivatives = 4;
  bkgdAdapter = NULL;
  bgU = NULL;
  bgWeights = NULL;
  sigmaTable = NULL;
}

// Destructor
VarDriver3D::~VarDriver3D()
{
}

// Initialize with the configHash already filled

bool VarDriver3D::initialize()
{
  // Run a 3D vortex background field
  cout << "Initializing SAMURAI 3D" << endl;
  return validateDriver();
}

// Initialize from parsed xml

bool VarDriver3D::initialize(const XMLNode& configuration)
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
  if (configHash["mode"] == "XYZ") {
    runMode = XYZ;
  } else if (configHash["mode"] == "RTZ") {
    runMode = RTZ;
  } else {
    cout << "Unrecognized run mode " << configHash["mode"] << ", Aborting...\n";
    return false;
  }

  bool fractlBkgd = ( configHash["bkgd_obs_interpolation"] == "fractl" );
  bool loadBG = ( configHash["load_background"] == "true" );

  // Warn is there are non-sensical combos of options.
  // Just 1 right now, so put it inlie

  if (fractlBkgd && ( configHash["adjust_background"] == "true"))
    std::cout << "** Warning: 'adjust_background' set to 'true' "
	      << "with 'bkgd_obs_interpolation' set to 'fractl' doesn't make sense."
	      << std::endl;

  if (fractlBkgd && ! loadBG) {
    std::cout << "** Warning: 'bkgd_obs_interpolation' set "
	      << "but 'load_background' set to false. Setting load_background to true."
	      << std::endl;
    loadBG = true;
  }

  // Set the projection (to be used by the cost functions)
  projection.setProjection(projectionFromConfig());

  /* Set the data path */
  // NOTE (NCAR): Originally we filtered these by 'files' and then sorted - this is doable with std::filesystem, but
  // that's not fully implemented in some compilers yet, so for the moment we're just assuming it's a directory
  dataPath = configHash["data_directory"];
  if (!DirectoryExists(dataPath)) {
    std::cout << "Can't find data directory: " << configHash["data_directory"] << endl;
    return false;
  }

  /* Check to make sure the output path exists */

  std::string outputPath(configHash["output_directory"]);
  if (!DirectoryExists(outputPath)) {
    std::cout << "Can't find output directory: " << configHash["output_directory"] << endl;
    return false;
  }

  // Centers and Met Observations are tightly coupled.
  //
  // So if we allow different centers between each runs, preProcessMetObs() or loadMetObs()
  // have to be called again.
  //
  // What else is dependent on the data structures created by readFrameCenters() and
  // the *MetObs() ?
  //
  // Are all the Met Obs still available, or only the ones that were matched to the center
  // time frames?

  // With coamps, grid dimensions come from the run arguments, not the fixed confid
  // So only do validation that depends on grid dimensions if fixedGrid is set.

  if (fractlBkgd) { // Grid comes from the fractl file
    if ( ! validateFractlGrid() )
      return false;

    bkgdAdapter = new BkgdFractl();
    return gridDependentInit();
  } else if (fixedGrid) {   // If runGrid, this will be done in the run(......) call
    if ( ! validateFixedGrid() )
      return false;

    // Set the background Obs adapter to get data from a file (if config says so)
    if (loadBG) {
      bkgdAdapter = new BkgdStream((dataPath + "/samurai_Background.in").c_str()); // Not cross-platform with the '/', but that's fixable later, easily.
    }
    return gridDependentInit();
  }

  // If we got here, then everything probably went OK!
  return true;
}

// This used to be in validateDriver.
// All the grid dependent initialization is performed here.
// Prerequisit: Grid dimensions have been set and verified. ***

bool VarDriver3D::gridDependentInit()
{

  START_TIMER(timei);

  // Define the Reference state

  std::string refSounding = configHash["ref_state"];
  refstate = new ReferenceState(refSounding);
  // cout << "Reference profile: Z\t\tQv\tRhoa\tRho\tH\tTemp\tPressure\n";
  for (real k = kmin; k < kmax + kincr; k += kincr) {
    // cout << "                   " << k << "\t";
    for (int i = 0; i < 6; i++) {
      real var = refstate->getReferenceVariable(i, k * 1000);
      if (i == 0) var = refstate->bhypInvTransform(var);
      // cout << setw(9) << setprecision(4)  << var << "\t";
    }
    // cout << "\n";
  }
  cout << setprecision(9);

  // Set the maximum number of iterations to the multipass reduction factor
  //  Multiple outer loops will reduce the cutoff wavelengths and background error variance

  maxIter = std::stoi(configHash["num_iterations"]);

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

  // These are used to process the obs (bkg and met)

  uStateSize = 8 * (idim + 1) * (jdim + 1) * (kdim + 1) * (numVars); // bgU size
  bStateSize = (idim + 2) * (jdim + 2) * (kdim + 2) * numVars;       // 2 mish points between nodes

  std::cout << "Physical (mish) State size = " << uStateSize << std::endl;
  std::cout << "Nodal State size = " << bStateSize << std::endl;
  std::cout << "Grid dimensions: (" << idim << ", " << jdim << ", " << kdim << ")" << std::endl;

  if (bgU != NULL)
    delete[] bgU;
  bgU = new real[uStateSize];
  std::memset(bgU, 0, uStateSize * sizeof(real)); // needed, as data not necessarily initialized to 0

  // Initialize this to zero.
  std::fill(bgU,bgU+uStateSize,0.0);
  //for (int i=0;i< uStateSize;i++) {bgU[i]=0.0;}

  // Optionally load a set of background coefficients directly

  std::string loadBGcoeffs = configHash["load_bg_coefficients"];

  if (loadBGcoeffs == "true")
    if (! loadBackgroundCoeffs() )
      return false;

  // Optionally load a set of background estimates and interpolate to the Gaussian mish

  int numbgObs = 0;
  if(bkgdAdapter != NULL) {
    START_TIMER(timeb);
    numbgObs = loadBackgroundObs();
    PRINT_TIMER("loadBackgroundObs", timeb);

    if (numbgObs < 0) {
      cout << "Error loading background Obs\n";
      return false;
    }
    //NCAR - cleanup from dangling bkgAdapter pointer, as it's not needed elsewhere:
    delete bkgdAdapter;
  }

  // Optionally adjust the interpolated background to satisfy mass continuity
  // and match the supplied points exactly. In essence, do a SAMURAI analysis using
  // the background estimates as "observations"

  std::string adjustBG = configHash["adjust_background"];
  if ((adjustBG == "true") and numbgObs) {
    if ( ! adjustBackground()) {
      cout << "Error adjusting background\n";
      return false;
    }
  }

  START_TIMER(timem);
  if ( ! loadMetObs() )
    return false;
  PRINT_TIMER("loadMetObs", timem);

  START_TIMER(timec);
  initObCost3D();
  PRINT_TIMER("initObCost3D", timec);

  PRINT_TIMER("gridDependentInit", timei);

  return true;
}

// Find the center that matches the ref_time

bool VarDriver3D::findReferenceCenter()
{
	// NOTE (NCAR) : I'm not sure whether we're JUST comparing times, or if dates matter?  Seems not, but get clarification from CSU team
  datetime reftime = ParseTime(configHash["ref_time"].c_str(), "%H:%M:%S");

  for (unsigned int fi = 0; fi < frameVector.size(); fi++) {
    datetime frametime = frameVector[fi].getTime();
    if (Time(reftime) == Time(frametime)) {
      configHash.insert("ref_lat", std::to_string(frameVector[fi].getLat()));
      configHash.insert("ref_lon", std::to_string(frameVector[fi].getLon()));
      // Note - we can't insert, since it already exists.. so we're doing the awkward 'GetMap', then assigning
      // a new value directly.  This can definitely be cleaner, but let's get it correct first, clean later.
      configHash.GetMap()->find("ref_time")->second = std::to_string(Date(frametime));
      cout << "Found matching reference time " << PrintTime(reftime) << " at " << frameVector[fi].getLat() << ", " << frameVector[fi].getLon() << "\n";
      return true;
    }
  }
  return false;
}

//
// Initialize background observations and run
//

// While this was developed to interface with coamps, nothing prevents a caller to use this interface.

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
					float *terrain_Height,
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
  if (configHash["array_order"] == "row-major")
    bkgdAdapter = new BkgdCArray(nx, ny, nsigma,
				 cdtg, delta, iter,
				 sigmas,
				 latitude, longitude, terrain_Height,
				 u1, v1, w1, th1, p1);
  else
    bkgdAdapter = new BkgdFArray(nx, ny, nsigma,
				 cdtg, delta, iter,
				 sigmas,
				 latitude, longitude, terrain_Height,
				 u1, v1, w1, th1, p1);
  if (! gridDependentInit() )
    return false;

  if( run() )
    return obCost3D->copyResults(nx, ny, nsigma, usam, vsam, wsam, thsam, psam);
  return false;
}


/* This routine drives the CostFunction minimization
 There is support for an outer loop to change the background
 error covariance or update non-linear observation operators */

bool VarDriver3D::run()
{
	int iter = 1;
	while (iter <= maxIter) {
		if (iter < maxIter) {
                       configHash.update("save_mish", "true");
		} else {
                       configHash.update("save_mish", "false");
		}
		cout << "Outer Loop Iteration: " << iter << endl;
		START_TIMER(timei);
		obCost3D->initState(iter);
		PRINT_TIMER("Cost3D Init", timei);

		START_TIMER(timem);
		obCost3D->minimize();
		PRINT_TIMER("Cost3D minimize", timem);

		START_TIMER(timeu);
		obCost3D->updateBG();
		PRINT_TIMER("Cost3d update", timeu);

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
// Set lower boundary info in the preProcessMetObs function
bool VarDriver3D::preProcessMetObs()
{
  GPTLstart("VarDriver3D::preprocessMetObs");

  vector<real> rhoP;

  // Convert the bg dBZ back to Z for further processing with real radar data

  if ((configHash["qr_variable"] == "dbz") and
      (configHash["load_background"] == "true") and
      (configHash["adjust_background"] == "false")) {
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
		real dBZ = bgU[bIndex +6] * 10. - 35.;
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

  real referenceLon = std::stof(configHash["ref_lon"]);

  // Find the zero C line using Newton's method

  real zeroClevel = 273.15;
  real height = 5000;
  real tmin = 1e34;
  int iter = 0;

  while ((fabs(tmin) > 0.1) and (iter < 5000)) {
    real t = refstate->getReferenceVariable(ReferenceVariable::tempref, height) - zeroClevel;
    real tprime = (refstate->getReferenceVariable(ReferenceVariable::tempref, height + 500.)
		   - refstate->getReferenceVariable(ReferenceVariable::tempref, height - 500.)) / 1000.;
    if (tprime != 0) {
      height = height - t / tprime;
      tmin = t;
    }
    iter++;
  }
  zeroClevel = height;
  cout << "Found zero C level at " << zeroClevel << " based on reference sounding" << endl;

  // Load Met. Observations
  auto filenames = FileList(dataPath);

  int processedFiles = 0;
  int attemptedFiles = 0;
  std::vector<MetObs>* metData = new std::vector<MetObs>;
  cout << "Found " << filenames.size() << " data files to read..." << endl;

  int totalFiles = filenames.size();
  for (std::size_t i = 0; i < filenames.size(); ++i) {
    metData->clear();

		std::string file = filenames[i];
    std::vector<std::string> parts = LineSplit(file, '.');
    std::string suffix = parts[parts.size()-1];
    std::string prefix = parts[0];


    if (prefix == "swp") {
      // Switch it to suffix
      suffix = "swp";
    }

    if (suffix == "nc") {	// cfrad file?
			if (std::regex_match(file, std::regex(".*cfrad.*\\.nc")))
				suffix = "cfrad";
    }

    cout << "Processing " << file << " of type " << suffix << endl;
    attemptedFiles++;
    // Read different types of files
    std::string fullpath = dataPath + "/" + file;
    if (! read_met_obs_file(dataSuffix[suffix], fullpath, metData))
      continue;

    processedFiles++;

    int obsProblem = 0;
    int coordProblem = 0;
    int timeProblem = 0;
    int domainProblem = 0;
    int radiusProblem = 0;

    // Process the metObs into Observations

    if (frameVector.size() == 0) {
      std::cout << "No centerfile, cannot process Met. Obs." << std::endl;
      return false;
    }

    datetime startTime_ob = frameVector.front().getTime();
    auto startTime = Date(startTime_ob);
    datetime endTime_ob = frameVector.back().getTime();
    auto endTime = Date(endTime_ob);
    int prevobs = obVector.size();
		std::cout << "i  = " << i << endl;
    for (std::size_t i = 0; i < metData->size(); ++i) {

      // Make sure the ob is within the time limits
      MetObs metOb = metData->at(i);
      datetime obTime_ob = metOb.getTime();
      auto obTime = Date(obTime_ob);


      // NOTE: Changing below line to account for msec vs. sec differences (discussion with M Bell, ongoing)
      // This makes the DesRosier case similar for now, but may need to change later
      //if ((obTime < startTime) or (obTime > endTime)) {
      if ((obTime < startTime) or (obTime >= endTime)) {
				timeProblem++;

				// if (timeProblem < 10) testing
				//   std::cout << "tcstart: " << PrintDate(startTime_ob) << ", tcend: " << PrintDate(endTime_ob) << ", obTime: " << PrintDate(obTime_ob) << std::endl;
				continue;
      }
			// std::cout << "timeProblem  = " << timeProblem << endl;
      int fi = std::chrono::duration_cast<std::chrono::seconds>(obTime_ob - startTime_ob).count();
			// std::cout << "fi = " << fi << endl;
      // if ((fi < 0) or (fi > (int)frameVector.size())) {
			// 	cout << "**Time problem with observation " << fi << ", " << startTime << ", " << obTime << endl;
			// 	timeProblem++;
			// 	continue;
      // } testing
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
      real obX = (metX - tcX) / 1000.;
      real obY = (metY - tcY) / 1000.;

      real heightm = metOb.getAltitude();

      real obZ = heightm/1000.;
      real obRadius = sqrt(obX*obX + obY*obY);
      real obTheta = 180.0 * atan2(obY, obX) / Pi;
      if (configHash["allow_negative_angles"] != "true")
	if (obTheta < 0)
	  obTheta += 360.0;

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
      varOb.setTime(obTime);  // Check this, given Qt's handling of time vs. ours (NCAR)

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
	  varOb.setError(std::stof(configHash["dropsonde_rhou_error"]));
	  obVector.push_back(varOb);

	  varOb.setWeight(0., 0);

	  varOb.setWeight(1., 1);
	  if (runMode == XYZ) {
	    rhov = rho*(v - Vm);
	  } else if (runMode == RTZ) {
	    rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	  }
	  varOb.setOb(rhov);
	  varOb.setError(std::stof(configHash["dropsonde_rhov_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 1);

	}
	if ((w != -999) and (rho != -999)) {
	  // rho w 1.5 m/s error
	  varOb.setWeight(1., 2);
	  rhow = rho*w;
	  varOb.setOb(rhow);
	  varOb.setError(std::stof(configHash["dropsonde_rhow_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 2);
	}
	if (tempk != -999) {
	  // temperature 1 K error
	  varOb.setWeight(1., 3);
	  varOb.setOb(tempk - tBar);
	  varOb.setError(std::stof(configHash["dropsonde_tempk_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 3);
	}
	if (qv != -999) {
	  // Qv 0.5 g/kg error
	  varOb.setWeight(1., 4);
	  qv = refstate->bhypTransform(qv);
	  varOb.setOb(qv-qBar);
	  varOb.setError(std::stof(configHash["dropsonde_qv_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 4);
	}
	if (rhoa != -999) {
	  // Rho prime .1 kg/m^3 error
	  varOb.setWeight(1., 5);
	  varOb.setOb((rhoa-rhoBar)*100);
	  varOb.setError(std::stof(configHash["dropsonde_rhoa_error"]));
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
	  varOb.setError(std::stof(configHash["flightlevel_rhou_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 0);

	  varOb.setWeight(1., 1);
	  if (runMode == XYZ) {
	    rhov = rho*(v - Vm);
	  } else if (runMode == RTZ) {
	    rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	  }
	  varOb.setOb(rhov);
	  varOb.setError(std::stof(configHash["flightlevel_rhov_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 1);
	}
	if ((w != -999) and (rho != -999)) {
	  // rho w 1 dm/s error
	  varOb.setWeight(1., 2);
	  rhow = rho*w;
	  varOb.setOb(rhow);
	  varOb.setError(std::stof(configHash["flightlevel_rhow_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 2);
	}
	if (tempk != -999) {
	  // temperature 1 K error
	  varOb.setWeight(1., 3);
	  varOb.setOb(tempk - tBar);
	  varOb.setError(std::stof(configHash["flightlevel_tempk_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 3);
	}
	if (qv != -999) {
	  // Qv 0.5 g/kg error
	  varOb.setWeight(1., 4);
	  qv = refstate->bhypTransform(qv);
	  varOb.setOb(qv-qBar);
	  varOb.setError(std::stof(configHash["flightlevel_qv_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 4);
	}
	if (rhoa != -999) {
	  // Rho prime .1 kg/m^3 error
	  varOb.setWeight(1., 5);
	  varOb.setOb((rhoa-rhoBar)*100);
	  varOb.setError(std::stof(configHash["flightlevel_rhoa_error"]));
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
	  varOb.setError(std::stof(configHash["insitu_rhou_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 0);

	  varOb.setWeight(1., 1);
	  if (runMode == XYZ) {
	    rhov = rho*(v - Vm);
	  } else if (runMode == RTZ) {
	    rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	  }
	  varOb.setOb(rhov);
	  varOb.setError(std::stof(configHash["insitu_rhov_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 1);

	}
	if ((w != -999) and (rho != -999)) {
	  // rho w 1.5 m/s error
	  varOb.setWeight(1., 2);
	  rhow = rho*w;
	  varOb.setOb(rhow);
	  varOb.setError(std::stof(configHash["insitu_rhow_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 2);
	}
	if (tempk != -999) {
	  // temperature 1 K error
	  varOb.setWeight(1., 3);
	  varOb.setOb(tempk - tBar);
	  varOb.setError(std::stof(configHash["insitu_tempk_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 3);
	}
	if (qv != -999) {
	  // Qv 0.5 g/kg error
	  varOb.setWeight(1., 4);
	  qv = refstate->bhypTransform(qv);
	  varOb.setOb(qv-qBar);
	  varOb.setError(std::stof(configHash["insitu_qv_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 4);
	}
	if (rhoa != -999) {
	  // Rho prime .1 kg/m^3 error
	  varOb.setWeight(1., 5);
	  varOb.setOb((rhoa-rhoBar)*100);
	  varOb.setError(std::stof(configHash["insitu_rhoa_error"]));
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
	  varOb.setError(metOb.getTemperatureError() + std::stof(configHash["mtp_tempk_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 3);
	}
	if (rhoa != -999) {
	  // Rho prime .1 kg/m^3 error
	  varOb.setWeight(1., 5);
	  varOb.setOb((rhoa-rhoBar)*100);
	  varOb.setError(std::stof(configHash["mtp_rhoa_error"]));
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
	varOb.setError(std::stof(configHash["sfmr_windspeed_error"]));
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
	  varOb.setError(std::stof(configHash["qscat_rhou_error"]));
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
	  varOb.setError(std::stof(configHash["qscat_rhov_error"]));
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
	  varOb.setError(std::stof(configHash["ascat_rhou_error"]));
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
	  varOb.setError(std::stof(configHash["ascat_rhov_error"]));
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
	  varOb.setError(std::stof(configHash["amv_rhou_error"]));
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
	  varOb.setError(std::stof(configHash["amv_rhov_error"]));
	  obVector.push_back(varOb);
	  varOb.setWeight(0., 1);
	}
	break;

			case (MetObs::terrain):
	{
		varOb.setType(MetObs::terrain);
		real dhdx = metOb.getTerrainDX();
		real dhdy = metOb.getTerrainDY();
		real drhoudx_coeff = -dhdx/sqrt(1+dhdx*dhdx+dhdy*dhdy);
		real drhoudy_coeff = -dhdy/sqrt(1+dhdx*dhdx+dhdy*dhdy);
		real drhoudz_coeff = 1/sqrt(1+dhdx*dhdx+dhdy*dhdy);

		real drhovdx_coeff = -dhdx/sqrt(1+dhdx*dhdx+dhdy*dhdy);
		real drhovdy_coeff = -dhdy/sqrt(1+dhdx*dhdx+dhdy*dhdy);
		real drhovdz_coeff = 1/sqrt(1+dhdx*dhdx+dhdy*dhdy);
		if (dhdy != -999) {
		// rho u 10 m/s error

// dudn = 0
		varOb.setWeight(drhoudx_coeff, 0, 1);
		varOb.setWeight(drhoudy_coeff, 0, 2);
		varOb.setWeight(drhoudz_coeff, 0, 3);
		varOb.setOb(0.0);
		varOb.setError(std::stof(configHash["neumann_u_weight"]));
		obVector.push_back(varOb);
		varOb.setWeight(0.0, 0, 1);
		varOb.setWeight(0.0, 0, 2);
		varOb.setWeight(0.0, 0, 3);
// dvdn = 0
		varOb.setWeight(drhovdx_coeff, 1, 1);
		varOb.setWeight(drhovdy_coeff, 1, 2);
		varOb.setWeight(drhovdz_coeff, 1, 3);
		varOb.setOb(0.0);
		varOb.setError(std::stof(configHash["neumann_v_weight"]));
		obVector.push_back(varOb);
		varOb.setWeight(0.0, 1, 1);
		varOb.setWeight(0.0, 1, 2);
		varOb.setWeight(0.0, 1, 3);
		// Dirichlet Boundary
		varOb.setWeight(dhdx, 0, 0);
		varOb.setWeight(dhdy, 1, 0);
		varOb.setWeight(-1  , 2, 0);
		varOb.setOb(0);
		varOb.setError(std::stof(configHash["dirichlet_w_weight"]));
		obVector.push_back(varOb);
		varOb.setWeight(0.0, 0, 0); // does it need to set weight?
		varOb.setWeight(0.0, 1, 0);
		varOb.setWeight(0.0, 2, 0);
		}
	break;
}

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
	  real DopplerError = metOb.getSpectrumWidth()*std::stof(configHash["lidar_sw_error"])
	    + log(std::stof(configHash["lidar_power_error"])/db);
	  if (DopplerError < std::stof(configHash["lidar_min_error"]))
	    DopplerError = std::stof(configHash["lidar_min_error"]);
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
	  if (configHash["horizontal_radar_appx"] == "true")
	    wWgt = 0;

	  // Fall speed
	  real Z = metOb.getReflectivity();
	  real w_term = 0.0;
	  real ZZ = -999.0;
	  if (Z > -999.0) {
	    real H = metOb.getAltitude();
	    ZZ=pow(10.0,(Z*0.1));
	    real melting_zone = 1000 * std::stof(configHash["melting_zone_width"]);
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
	    real mixed_dbz = std::stof(configHash["mixed_phase_dbz"]);
	    real rain_dbz = std::stof(configHash["rain_dbz"]);
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
	    real DopplerError = fabs(wWgt)*std::stof(configHash["radar_fallspeed_error"]);
	    if (metOb.getSpectrumWidth() != -999.0) {
	      DopplerError += metOb.getSpectrumWidth()*std::stof(configHash["radar_sw_error"]);
	    }
	    if (DopplerError < std::stof(configHash["radar_min_error"]))
	      DopplerError = std::stof(configHash["radar_min_error"]);
	    varOb.setError(DopplerError);
	    varOb.setOb(Vdopp);

      real maxel;
      if (configHash.exists("max_radar_elevation") == false) {
         maxel = 90.0;
      } else {
         maxel = std::stof(configHash["max_radar_elevation"]);
       }
	    if (fabs(metOb.getElevation() <= maxel))
	      obVector.push_back(varOb);

	    varOb.setWeight(0., 0);
	    varOb.setWeight(0., 1);
	    varOb.setWeight(0., 2);
	  }

	  // Reflectivity observations
	  std::string gridref = configHash["qr_variable"];
	  real qr = 0.;
	  if (ZZ > 0) {
	    if (gridref == "qr") {
	      // Do the gridding as part of the variational synthesis using Z-M relationships
	      // Z-M relationships from Gamache et al (1993) JAS
	      real H = metOb.getAltitude();
	      real melting_zone = 1000 * std::stof(configHash["melting_zone_width"]);
	      real hlow= zeroClevel;
	      real hhi= hlow + melting_zone;
	      real rainmass = pow(ZZ/14630.,(real)0.6905);
	      real icemass = pow(ZZ/670.,(real)0.5587);
	      real mixed_dbz = std::stof(configHash["mixed_phase_dbz"]);
	      real rain_dbz = std::stof(configHash["rain_dbz"]);
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
	    real iROI = std::stof(configHash["i_reflectivity_roi"]) / iincr;
	    real jROI = std::stof(configHash["j_reflectivity_roi"]) / jincr;
	    real kROI = std::stof(configHash["k_reflectivity_roi"]) / kincr;
	    real Rsquare = (iincr*iROI)*(iincr*iROI) + (jincr*jROI)*(jincr*jROI) + (kincr*kROI)*(kincr*kROI);
#pragma omp parallel for
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
			int64_t bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
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
	    varOb.setError(std::stof(configHash["mesonet_rhou_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 0);

	    varOb.setWeight(1., 1);
	    if (runMode == XYZ) {
	      rhov = rho*(v - Vm);
	    } else if (runMode == RTZ) {
	      rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	    }
	    varOb.setOb(rhov);
	    varOb.setError(std::stof(configHash["mesonet_rhov_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 1);

	  }
	  if ((w != -999) and (rho != -999)) {
	    // rho w 1.5 m/s error
	    varOb.setWeight(1., 2);
	    rhow = rho*w;
	    varOb.setOb(rhow);
	    varOb.setError(std::stof(configHash["mesonet_rhow_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 2);
	  }
	  if (tempk != -999) {
	    // temperature 1 K error
	    varOb.setWeight(1., 3);
	    varOb.setOb(tempk - tBar);
	    varOb.setError(std::stof(configHash["mesonet_tempk_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 3);
	  }
	  if (qv != -999) {
	    // Qv 0.5 g/kg error
	    varOb.setWeight(1., 4);
	    qv = refstate->bhypTransform(qv);
	    varOb.setOb(qv-qBar);
	    varOb.setError(std::stof(configHash["mesonet_qv_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 4);
	  }
	  if (rhoa != -999) {
	    // Rho prime .1 kg/m^3 error
	    varOb.setWeight(1., 5);
	    varOb.setOb((rhoa-rhoBar)*100);
	    varOb.setError(std::stof(configHash["mesonet_rhoa_error"]));
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
	    varOb.setError(std::stof(configHash["aeri_rhou_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 0);

	    varOb.setWeight(1., 1);
	    if (runMode == XYZ) {
	      rhov = rho*(v - Vm);
	    } else if (runMode == RTZ) {
	      rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
	    }
	    varOb.setOb(rhov);
	    varOb.setError(std::stof(configHash["aeri_rhov_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 1);

	  }
	  if ((w != -999) and (rho != -999)) {
	    // rho w 1.5 m/s error
	    varOb.setWeight(1., 2);
	    rhow = rho*w;
	    varOb.setOb(rhow);
	    varOb.setError(std::stof(configHash["aeri_rhow_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 2);
	  }
	  if (tempk != -999) {
	    // temperature 1 K error
	    varOb.setWeight(1., 3);
	    varOb.setOb(tempk - tBar);
	    varOb.setError(std::stof(configHash["aeri_tempk_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 3);
	  }
	  if (qv != -999) {
	    // Qv 0.5 g/kg error
	    varOb.setWeight(1., 4);
	    qv = refstate->bhypTransform(qv);
	    varOb.setOb(qv-qBar);
	    varOb.setError(std::stof(configHash["aeri_qv_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 4);
	  }
	  if (rhoa != -999) {
	    // Rho prime .1 kg/m^3 error
	    varOb.setWeight(1., 5);
	    varOb.setOb((rhoa-rhoBar)*100);
	    varOb.setError(std::stof(configHash["aeri_rhoa_error"]));
	    obVector.push_back(varOb);
	    varOb.setWeight(0., 5);
	  }
	  break;
	}

	case (MetObs::crsim):
		{
		varOb.setType(MetObs::crsim);

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
		if (configHash["horizontal_radar_appx"] == "true")
			wWgt = 0;

		// Fall speed
		real Z = metOb.getReflectivity();
		real w_term = 0.0;
		real ZZ = -999.0;
		if (Z > -999.0) {
			real H = metOb.getAltitude();
			ZZ=pow(10.0,(Z*0.1));
			real melting_zone = 1000 * std::stof(configHash["melting_zone_width"]);
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
			real mixed_dbz = std::stof(configHash["mixed_phase_dbz"]);
			real rain_dbz = std::stof(configHash["rain_dbz"]);
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
			real DopplerError = fabs(wWgt)*std::stof(configHash["radar_fallspeed_error"]);
			if (DopplerError < std::stof(configHash["radar_min_error"]))
				DopplerError = std::stof(configHash["radar_min_error"]);
			varOb.setError(DopplerError);
			varOb.setOb(Vdopp);

			real maxel;
			if (configHash.exists("max_radar_elevation") == false) {
				 maxel = 90.0;
			} else {
				 maxel = std::stof(configHash["max_radar_elevation"]);
			 }
			if (fabs(metOb.getElevation() <= maxel))
				obVector.push_back(varOb);

			varOb.setWeight(0., 0);
			varOb.setWeight(0., 1);
			varOb.setWeight(0., 2);
		}

		// Reflectivity observations
		std::string gridref = configHash["qr_variable"];
		real qr = 0.;
		if (ZZ > 0) {
			if (gridref == "qr") {
				// Do the gridding as part of the variational synthesis using Z-M relationships
				// Z-M relationships from Gamache et al (1993) JAS
				real H = metOb.getAltitude();
				real melting_zone = 1000 * std::stof(configHash["melting_zone_width"]);
				real hlow= zeroClevel;
				real hhi= hlow + melting_zone;
				real rainmass = pow(ZZ/14630.,(real)0.6905);
				real icemass = pow(ZZ/670.,(real)0.5587);
				real mixed_dbz = std::stof(configHash["mixed_phase_dbz"]);
				real rain_dbz = std::stof(configHash["rain_dbz"]);
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
			real iROI = std::stof(configHash["i_reflectivity_roi"]) / iincr;
			real jROI = std::stof(configHash["j_reflectivity_roi"]) / jincr;
			real kROI = std::stof(configHash["k_reflectivity_roi"]) / kincr;
			real Rsquare = (iincr*iROI)*(iincr*iROI) + (jincr*jROI)*(jincr*jROI) + (kincr*kROI)*(kincr*kROI);
		#pragma omp parallel for
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
			int64_t bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
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

	case(MetObs::model):
			varOb.setType(MetObs::model);
			u = metOb.getZonalVelocity();
			v = metOb.getMeridionalVelocity();
			w = metOb.getVerticalVelocity();
			rho = metOb.getModelMoistDensity();
			rhoa = metOb.getModelAirDensity();
			qv = metOb.getModelQv();
			tempk = metOb.getTemperature();
			real Z = metOb.getReflectivity();
			real ZZ=pow(10.0,(Z*0.1));

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
		    varOb.setError(1);
		    obVector.push_back(varOb);
		    varOb.setWeight(0., 0);

		    varOb.setWeight(1., 1);
		    if (runMode == XYZ) {
		      rhov = rho*(v - Vm);
		    } else if (runMode == RTZ) {
		      rhov = rho*(-(u - Um)*obY + (v - Vm)*obX)/obRadius;
		    }
		    varOb.setOb(rhov);
		    varOb.setError(1);
		    obVector.push_back(varOb);
		    varOb.setWeight(0., 1);

		  }
		  if ((w != -999) and (rho != -999)) {
		    // rho w 1.5 m/s error
		    varOb.setWeight(1., 2);
		    rhow = rho*w;
		    varOb.setOb(rhow);
		    varOb.setError(1);
		    obVector.push_back(varOb);
		    varOb.setWeight(0., 2);
		  }

			if (tempk != -999) {
				// temperature 1 K error
				varOb.setWeight(1., 3);
				varOb.setOb(tempk - tBar);
				varOb.setError(std::stof(configHash["aeri_tempk_error"]));
				obVector.push_back(varOb);
				varOb.setWeight(0., 3);
			}
			if (qv != -999) {
				// Qv 0.5 g/kg error
				varOb.setWeight(1., 4);
				qv = refstate->bhypTransform(qv);
				varOb.setOb(qv-qBar);
				varOb.setError(std::stof(configHash["aeri_qv_error"]));
				obVector.push_back(varOb);
				varOb.setWeight(0., 4);
			}
			if (rhoa != -999) {
				// Rho prime .1 kg/m^3 error
				varOb.setWeight(1., 5);
				varOb.setOb((rhoa-rhoBar)*100);
				varOb.setError(std::stof(configHash["aeri_rhoa_error"]));
				obVector.push_back(varOb);
				varOb.setWeight(0., 5);
			}

			// Reflectivity observations
		  std::string gridref = configHash["qr_variable"];
		  real qr = 0.;
		  if (ZZ > 0) {
		    if (gridref == "qr") {
		      // Do the gridding as part of the variational synthesis using Z-M relationships
		      // Z-M relationships from Gamache et al (1993) JAS
		      real H = metOb.getAltitude();
		      real melting_zone = 1000 * std::stof(configHash["melting_zone_width"]);
		      real hlow= zeroClevel;
		      real hhi= hlow + melting_zone;
		      real rainmass = pow(ZZ/14630.,(real)0.6905);
		      real icemass = pow(ZZ/670.,(real)0.5587);
		      real mixed_dbz = std::stof(configHash["mixed_phase_dbz"]);
		      real rain_dbz = std::stof(configHash["rain_dbz"]);
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
		    real iROI = std::stof(configHash["i_reflectivity_roi"]) / iincr;
		    real jROI = std::stof(configHash["j_reflectivity_roi"]) / jincr;
		    real kROI = std::stof(configHash["k_reflectivity_roi"]) / kincr;
		    real Rsquare = (iincr*iROI)*(iincr*iROI) + (jincr*jROI)*(jincr*jROI) + (kincr*kROI)*(kincr*kROI);
	#pragma omp parallel for
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
				int64_t bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
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

    } // for everything in metData

    // Show a summary of what got tossed out

    cout << "Observation problem: " << obsProblem << ", Time problem: " << timeProblem
	 << ", Coordinate problem: " <<  coordProblem << ", Domain problem: " << domainProblem
	 << ", Radius problem: " << radiusProblem << endl;
    std::cout << "obVector size: " << obVector.size() << std::endl;

    int newobs = obVector.size() - prevobs;
    // if (metData->size() > 0) {
    if (newobs > 0) {
      cout << "Processed " << newobs << " observations from " << metData->size() << " entries ("
	   << 100.0*(float)newobs/(6*(float)metData->size()) << "%) file: " << file << std::endl;
    } else {
      cout << "No valid observations in file\n";
    }
    cout << obVector.size() << " total observations." << " ( " << attemptedFiles << " of " << totalFiles << " files processed ) " << endl;
  }

  delete metData;

  // Finish reflectivity interpolation

  Observation varOb;
  varOb.setTime(std::stoi(configHash["ref_time"]));
  real pseudow_weight = std::stof(configHash["dbz_pseudow_weight"]);
  real mc_weight = std::stof(configHash["mc_weight"]);

	// Initialize a MetObs for terrain
// 	std::vector<MetObs>* terrainData = new std::vector<MetObs>;
	std::string fullpath = dataPath + "/" + "terrain.hgt";
    std::ifstream terrainFile(fullpath);

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
		    if (ihalf and jhalf and khalf and (imu != 0) and (jmu != 0) and (kmu != 0)){
		      int bgI = (iIndex+1)*2 + (imu+1)/2;
		      int bgJ = (jIndex+1)*2 + (jmu+1)/2;
		      int bgK = (kIndex+1)*2 + (kmu+1)/2;
		      			int64_t bIndex = numVars*(idim+1)*2*(jdim+1)*2*bgK + numVars*(idim+1)*2*bgJ +numVars*bgI;
		      if (bgWeights[bIndex] != 0) {
			bgU[bIndex +6] /= bgWeights[bIndex];
		      }
		      if (configHash["qr_variable"] == "dbz") {
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
    
		// // Set a lower boundary condition for W
		// // Ideally use a terrain map here, but just use Z=0 for now
        
		if ((pseudow_weight > 0.0) and (!terrainFile.is_open())) {
			std::cout << "No input terrain file ... setting the lower boundary z = 0" <<std::endl;
			varOb.setAltitude(0);
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

#if IO_WRITEOBS
  GPTLstart("VarDriver3D::preprocessMetObs->writeobs");
  // Write the Obs to a summary text file
  std::string obFilename = dataPath + "/samurai_Observations.in";
  ofstream obstream(obFilename);
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
  GPTLstop("VarDriver3D::preprocessMetObs->writeobs");
#endif

  // Load the observations into a vector
  int64_t vector_size = (obVector.size() * (7 + numVars * numDerivatives));
  obs = new real[vector_size];
  for (int64_t m=0; m < obVector.size(); m++) {
    int64_t n = m * (7 + numVars * numDerivatives);
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
	      int64_t wgt_index = n + ( 7 * (d + 1)) + var;
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

  GPTLstop("VarDriver3D::preprocessMetObs");
  return true;
}

/* Load the meteorological observations from a file into a vector */

bool VarDriver3D::loadPreProcessMetObs()
{
    GPTLstart("VarDriver3D::loadPreprocessMetObs");
    Observation varOb;
    real wgt[numVars][4];
    real iPos, jPos, kPos, ob, error;
    int type;
    int64_t time;
    cout << "Loading preprocessed observations from samurai_Observations.in" << endl;

    // Open and read the file
    std::string obFilename = dataPath + "/samurai_Observations.in";
    ifstream obstream(obFilename);
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
                varOb.setWeight(wgt[var][d], var, d);
            }
        }
        obVector.push_back(varOb);
    }

    // Load the observations into the vector
    obs = new real[obVector.size() * (7 + numVars * numDerivatives)];
    for (int m=0; m < obVector.size(); m++) {
        int n = m * (7 + numVars * numDerivatives);
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
                int wgt_index = n + (7 * (d + 1)) + var;
                obs[wgt_index] = ob.getWeight(var, d);
            }
        }
    }

    GPTLstop("VarDriver3D::loadPreprocessMetObs");
    return true;
}


// Background Observations can come from
// - a samurai_Background.in file,
// - a FRACTL generated netcdf file
// - a WRF netcdf flie
// - from passed arguments to the run() function.

int VarDriver3D::loadBackgroundObs()
{
  if(bkgdAdapter == NULL) {
    std::cerr << "Error: VarDriver3D::loadBackgroundObs() called without adapter initialization."
	      << std::endl;
    exit(1);
  }

  // Get a background obs loader to load the background observations
  // Default is Spline.

  BkgdObsLoader::bg_loader_t loaderType =  BkgdObsLoader::BG_LOADER_SPLINE;
  if (configHash["bkgd_obs_interpolation"] == "kd_tree")
    loaderType = BkgdObsLoader::BG_LOADER_KD;
  else if (configHash["bkgd_obs_interpolation"] == "fractl")
    loaderType = BkgdObsLoader::BG_LOADER_FRACTL;
  BkgdObsLoader *bkgdObsLoader = BkgdObsLoaderFactory::createBkgdObsLoader(loaderType);
  if (bkgdObsLoader == NULL)
    return -1;

  bkgdObsLoader->initialize(&configHash, frameVector, bkgdAdapter, projection,
			    refstate, uStateSize, bgU,
			    numVars,
			    idim, jdim, kdim,
			    imin, jmin, kmin,
			    imax, jmax, kmax,
			    iincr, jincr, kincr);


  if (! bkgdObsLoader->loadBkgdObs(bgIn)) {
    std::cerr << "Failed to load background observations" << std::endl;
    return -1;
  }

  return  bgIn.size() * 7 / 11;
}

bool VarDriver3D::adjustBackground()
{
  // Set the minimum filter length to the background resolution, not the analysis resolution
  // to avoid artifacts when running interpolating to small mesoscale grids

  START_TIMER(timeab);

  // Load the observations into a vector
  int numbgObs = bgIn.size() * 7 / 11;
  if (std::stof(configHash["mc_weight"]) > 0) {
    numbgObs += idim * jdim  * kdim;
  }
  bgObs = new real[numbgObs * (7 + numVars * numDerivatives)];
  for (unsigned int m = 0; m < numbgObs * (7 + numVars * numDerivatives); m++)
    bgObs[m] = 0.;

  int p = 0;
  real obX, obY, obRadius, obTheta;
  obX = obY = obRadius = obTheta = -32768.;
  for (int m = 0; m < bgIn.size(); m += 11) {
    if (runMode == XYZ) {
      obX = bgIn[m];
      obY = bgIn[m + 1];
    } else if (runMode == RTZ) {
      obRadius = bgIn[m];
      obTheta = bgIn[m + 1];
    }
    real obZ = exp(bgIn[m + 2]);
    real obTime = bgIn[m + 3];
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
      bgObs[p] = bgIn[m + 4 + n];
      if ((n == 6) and (configHash["qr_variable"] == "dbz")) {
	// Convert to dBZ control variable
	real dbzavg = 10 * log10(bgIn[m + 4 + n]);
	bgObs[p] = (dbzavg + 35.) * 0.1;
      }
      // Default error of background = 0.1
      if ((configHash.exists("bg_obs_error") or (std::stof(configHash["bg_obs_error"]) <= 0.0))) {
	bgObs[p + 1] = 100.;
      } else {
	bgObs[p + 1] = 1.0 / std::stof(configHash["bg_obs_error"]);
      }
      if (runMode == XYZ) {
	bgObs[p + 2] = obX;
	bgObs[p + 3] = obY;
      } else if (runMode == RTZ) {
	bgObs[p + 2] = obRadius;
	bgObs[p + 3] = obTheta;
      }
      bgObs[p + 4] = obZ;
      // Null type
      bgObs[p + 5] = -1;
      bgObs[p + 6] = obTime;
      bgObs[p + 7 + n] = 1.;
      p += (7 + numVars * numDerivatives);
    }
  }

  // Add mass continuity constraint
  if (std::stof(configHash["mc_weight"]) > 0) {
    for (int iIndex = 0; iIndex < idim; iIndex++) {
      real i = imin + iincr * iIndex;
      if (i > ((idim - 1) * iincr + imin)) continue;
      for (int jIndex = 0; jIndex < jdim; jIndex++) {
	real j = jmin + jincr * jIndex;
	if (j > ((jdim - 1) * jincr + jmin)) continue;
	for (int kIndex = 0; kIndex < kdim; kIndex++) {
	  real k = kmin + kincr * kIndex;
	  if (k > ((kdim - 1)*kincr + kmin)) continue;
	  bgObs[p] = 0.0;
	  bgObs[p+1] = std::stof(configHash["mc_weight"]);
	  bgObs[p+2] = i;
	  bgObs[p+3] = j;
	  bgObs[p+4] = k;
	  // Null type
	  bgObs[p + 5] = -1;
	  bgObs[p + 6] = std::stoi(configHash["ref_time"]);
	  if (runMode == XYZ) {
	    bgObs[p + (7 * (1 + 1))] = 1.0;
	    bgObs[p + (7 * (2 + 1)) + 1] = 1.0;
	    bgObs[p + (7 * (3 + 1)) + 2] = 1.0;
	  } else if (runMode == RTZ) {
	    if (i > 0) {
	      real rInverse = 180.0/(i * Pi);
	      bgObs[p + 7] = 1.0/i;
	      bgObs[p + (7 * (1 + 1))] = 1.0;
	      bgObs[p + (7 * (2 + 1)) + 1] = rInverse;
	      bgObs[p + (7 * (3 + 1)) + 2] = 1.0;
	    }
	  }
	  p  += (7 + numVars * numDerivatives);
	}
      }
    }
  }

  // Store and set the background errors
  std::string bgError[7];
  bgError[0] = configHash["bg_rhou_error"];
  bgError[1] = configHash["bg_rhov_error"];
  bgError[2] = configHash["bg_rhow_error"];
  bgError[3] = configHash["bg_tempk_error"];
  bgError[4] = configHash["bg_qv_error"];
  bgError[5] = configHash["bg_rhoa_error"];
  bgError[6] = configHash["bg_qr_error"];

  std::string bg_interpolation_error = "1.0";
  if (configHash.exists("bg_interpolation_error")) {
    bg_interpolation_error = configHash["bg_interpolation_error"];
    cout << "Setting background interpolation error to " << bg_interpolation_error << "\n";
  } else {
    cout << "Using default background interpolation error of 1.0\n";
  }
  configHash.update("bg_rhou_error", bg_interpolation_error);
  configHash.update("bg_rhov_error", bg_interpolation_error);
  configHash.update("bg_rhow_error", bg_interpolation_error);
  configHash.update("bg_tempk_error", bg_interpolation_error);
  configHash.update("bg_qv_error", bg_interpolation_error);
  configHash.update("bg_rhoa_error", bg_interpolation_error);
  configHash.update("bg_qr_error", bg_interpolation_error);
  configHash.update("save_mish", "true");

  // Adjust the background field to the spline mish

  if (runMode == XYZ) {
    if (std::stof(configHash["output_pressure_increment"]) > 0) {
      bgCost3D = new CostFunctionXYP(projection, numbgObs, bStateSize);
    } else {
      bgCost3D = new CostFunctionXYZ(projection, numbgObs, bStateSize);
    }
  } else if (runMode == RTZ) {
    bgCost3D = new CostFunctionRTZ(projection, numbgObs, bStateSize);
  }
  bgCost3D->initialize(&configHash, bgU, bgObs, refstate);

  // Set the iteration to zero --
  // this will prevent writing the background file until after the adjustment
  // which is presumably what you want most of the time. Otherwise, you would not be here

  int bgIter = 1;
  bgCost3D->initState(bgIter);
  bgCost3D->minimize();

  // Increment the variables
  bgCost3D->updateBG();
  bgCost3D->finalize();

  delete bgCost3D;
  delete[] bgObs;

  // Reset the background errors

  configHash.update("bg_rhou_error", bgError[0]);
  configHash.update("bg_rhov_error", bgError[1]);
  configHash.update("bg_rhow_error", bgError[2]);
  configHash.update("bg_tempk_error", bgError[3]);
  configHash.update("bg_qv_error", bgError[4]);
  configHash.update("bg_rhoa_error", bgError[5]);
  configHash.update("bg_qr_error", bgError[6]);
  configHash.update("save_mish", "false");

  // Convert the dBZ back to Z for further processing

  if (configHash["qr_variable"] == "dbz") {
    for (int ki = -1; ki < (kdim); ki++) {
      for (int kmu = -1; kmu <= 1; kmu += 2) {
	for (int ii = -1; ii < (idim); ii++) {
	  for (int imu = -1; imu <= 1; imu += 2) {
	    for (int ji = -1; ji < (jdim); ji++) {
	      for (int jmu = -1; jmu <= 1; jmu += 2) {
		int bgI = (ii + 1) * 2 + (imu + 1) / 2;
		int bgJ = (ji + 1) * 2 + (jmu + 1) / 2;
		int bgK = (ki + 1) * 2 + (kmu + 1) / 2;
		int bIndex = numVars * (idim + 1) * 2 * (jdim + 1) * 2 * bgK
		  + numVars * (idim + 1) * 2 * bgJ + numVars * bgI;
		real dbZ = bgU[bIndex + 6] * 10. - 35.;
		real ZZ = pow(10.0, (dbZ * 0.1));
		bgU[bIndex + 6] = ZZ;
	      }
	    }
	  }
	}
      }
    }
  }
  PRINT_TIMER("adjustBackground", timeab);

  return true;
}

/* Any updates needed for additional analysis iterations go here */

void VarDriver3D::updateAnalysisParams(const int& iteration)
{
    std::string iter = std::to_string(iteration);

    std::string key = "bg_rhou_error_" + iter;
    std::string val = configHash[key];
    configHash.update("bg_rhou_error", val);

    key = "bg_rhov_error_" + iter;
    val = configHash[key];
    configHash.update("bg_rhov_error", val);

    key = "bg_rhow_error_" + iter;
    val = configHash[key];
    configHash.update("bg_rhow_error", val);

    key = "bg_tempk_error_" + iter;
    val = configHash[key];
    configHash.update("bg_tempk_error", val);

    key = "bg_qv_error_" + iter;
    val = configHash[key];
    configHash.update("bg_qv_error", val);

    key = "bg_rhoa_error_" + iter;
    val = configHash[key];
    configHash.update("bg_rhoa_error", val);

    key = "bg_qr_error_" + iter;
    val = configHash[key];
    configHash.update("bg_qr_error", val);

    key = "mc_weight_" + iter;
    val = configHash[key];
    configHash.update("mc_weight", val);

    key = "i_filter_length_" + iter;
    val = configHash[key];
    configHash.update("i_filter_length", val);

    key = "j_filter_length_" + iter;
    val = configHash[key];
    configHash.update("j_filter_length", val);

    key = "k_filter_length_" + iter;
    val = configHash[key];
    configHash.update("k_filter_length", val);

    key = "i_spline_cutoff_" + iter;
    val = configHash[key];
    configHash.update("i_spline_cutoff", val);

    key = "j_spline_cutoff_" + iter;
    val = configHash[key];
    configHash.update("j_spline_cutoff", val);

    key = "k_spline_cutoff_" + iter;
    val = configHash[key];
    configHash.update("k_spline_cutoff", val);
}

/* This routine validates that all required parameters are present
 It currently does not check the validity of a particular parameter, just that it exists */

bool VarDriver3D::validateConfig()
{
    // Validate the hash -- multiple passes are not validated currently
    std::set<std::string> configKeys;
    configKeys.insert("mc_weight");
    configKeys.insert("i_filter_length");
    configKeys.insert("j_filter_length");
    configKeys.insert("k_filter_length");
    configKeys.insert("i_spline_cutoff");
    configKeys.insert("j_spline_cutoff");
    configKeys.insert("k_spline_cutoff");
    configKeys.insert("i_rhou_bcL");
    configKeys.insert("i_rhou_bcR");
    configKeys.insert("j_rhou_bcL");
    configKeys.insert("j_rhou_bcR");
    configKeys.insert("k_rhou_bcL");
    configKeys.insert("k_rhou_bcR");
    configKeys.insert("i_rhov_bcL");
    configKeys.insert("i_rhov_bcR");
    configKeys.insert("j_rhov_bcL");
    configKeys.insert("j_rhov_bcR");
    configKeys.insert("k_rhov_bcL");
    configKeys.insert("k_rhov_bcR");
    configKeys.insert("i_rhow_bcL");
    configKeys.insert("i_rhow_bcR");
    configKeys.insert("j_rhow_bcL");
    configKeys.insert("j_rhow_bcR");
    configKeys.insert("k_rhow_bcL");
    configKeys.insert("k_rhow_bcR");
    configKeys.insert("i_tempk_bcL");
    configKeys.insert("i_tempk_bcR");
    configKeys.insert("j_tempk_bcL");
    configKeys.insert("j_tempk_bcR");
    configKeys.insert("k_tempk_bcL");
    configKeys.insert("k_tempk_bcR");
    configKeys.insert("i_qv_bcL");
    configKeys.insert("i_qv_bcR");
    configKeys.insert("j_qv_bcL");
    configKeys.insert("j_qv_bcR");
    configKeys.insert("k_qv_bcL");
    configKeys.insert("k_qv_bcR");
    configKeys.insert("i_rhoa_bcL");
    configKeys.insert("i_rhoa_bcR");
    configKeys.insert("j_rhoa_bcL");
    configKeys.insert("j_rhoa_bcR");
    configKeys.insert("k_rhoa_bcL");
    configKeys.insert("k_rhoa_bcR");
    configKeys.insert("i_qr_bcL");
    configKeys.insert("i_qr_bcR");
    configKeys.insert("j_qr_bcL");
    configKeys.insert("j_qr_bcR");
    configKeys.insert("k_qr_bcL");
    configKeys.insert("k_qr_bcR");
    configKeys.insert("data_directory");
    configKeys.insert("output_directory");

    if (fixedGrid) {
    	configKeys.insert("i_min");
    	configKeys.insert("i_max");
    	configKeys.insert("i_incr");
    	configKeys.insert("j_min");
    	configKeys.insert("j_max");
    	configKeys.insert("j_incr");
    	configKeys.insert("k_min");
    	configKeys.insert("k_max");
    	configKeys.insert("k_incr");
		}

    for (auto &key : configKeys) {
        if ( configHash.exists(key) == false ) {
            std::cout <<	"No configuration found for <" << key << "> aborting..." << std::endl;
            return false;
        }
    }

    // Add default values here

    if ( configHash.exists("bkgd_obs_interpolation") == false)
      configHash.insert("bkgd_obs_interpolation", "spline");

    if ( configHash.exists("bkgd_kd_num_neighbors") == false)
      configHash.insert("bkgd_kd_num_neighbors", "6");

    if ( configHash.exists("bkgd_kd_max_distance") == false) {
      // TODO What should the default max distance for a nearest neighbor be?
      configHash.insert("bkgd_kd_max_distance", "100");
    }

    // All done

    return true;
}

bool VarDriver3D::loadBackgroundCoeffs()
{
	// Load a set of coefficients directly for the same grid
	cout << "Loading previous coefficients from samurai_Coefficients.in" << endl;

	std::ifstream bgFile(dataPath + "/samurai_Coefficients.in");
	if (!bgFile.is_open())
		return false;

	int numbgCoeffs = 0;
	std::string line;

	while (std::getline(bgFile, line)) {
    std::istringstream iss(line);
		if (line.find("Variable") != std::string::npos) {
			continue;
		} else {
			int var;
			iss >> var;
			int iIndex;
			iss >> iIndex;
			int jIndex;
			iss >> jIndex;
			int kIndex;
			iss >> kIndex;
			int bIndex = numVars*(idim+2)*(jdim+2)*kIndex + numVars*(idim+2)*jIndex +numVars*iIndex + var;
			float value;
			iss >> value;
			bgU[bIndex] = value;
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
  // Grid specs

  imin = i_min;
  imax = i_max;
  iincr = i_incr;

  jmin = j_min;
  jmax = j_max;
  jincr = j_incr;

  // k grid dimension. This is pretty arbitrary
  kmin = 0.0;
  // kmax = sigmas[0] / 1000;	// do we ned to round this up?
  kmax = *std::max_element(sigmas, sigmas + (int) n_z) / 1000;
  kincr = 0.5;

  // Array dimensions

  idim = n_x;
  jdim = n_y;
  kdim = n_z;

  sigmaTable = sigmas;

  return validateGrid();
}

bool VarDriver3D::validateFractlGrid()
{
  // Need to read grid from fractl file.
  Nc3Error err(Nc3Error::verbose_nonfatal);
  std::string fname = configHash["fractl_nc_file"];

   // Open the file.
  Nc3File dataFile(fname.c_str(), Nc3File::ReadOnly);

   // Check to see if the file was opened.
   if(!dataFile.is_valid()) {
     std::cout << "Failed to read FRACTL generated nc file " << fname << std::endl;
     return false;
   }

   Nc3Dim *timeDim = dataFile.get_dim("time");
   if (! timeDim)
     return false;

   Nc3Dim *z0Dim = dataFile.get_dim("z0");
   if (! z0Dim)
     return false;

   Nc3Dim *y0Dim = dataFile.get_dim("y0");
   if (! y0Dim)
     return false;

   Nc3Dim *x0Dim = dataFile.get_dim("x0");
   if (! x0Dim)
     return false;


   Nc3Att *sam_idim = dataFile.get_att("sam_idim");
   if (! sam_idim)
     return false;
   Nc3Att *sam_jdim = dataFile.get_att("sam_jdim");
   if (! sam_jdim)
     return false;
   Nc3Att *sam_kdim = dataFile.get_att("sam_kdim");
   if (! sam_kdim)
     return false;

   Nc3Att *sam_imin = dataFile.get_att("sam_imin");
   if (! sam_imin)
     return false;
   Nc3Att *sam_imax = dataFile.get_att("sam_imax");
   if (! sam_imax)
     return false;
   Nc3Att *sam_iincr = dataFile.get_att("sam_iincr");
   if (! sam_iincr)
     return false;

   Nc3Att *sam_jmin = dataFile.get_att("sam_jmin");
   if (! sam_jmin)
     return false;
   Nc3Att *sam_jmax = dataFile.get_att("sam_jmax");
   if (! sam_jmax)
     return false;
   Nc3Att *sam_jincr = dataFile.get_att("sam_jincr");
   if (! sam_jincr)
     return false;

   Nc3Att *sam_kmin = dataFile.get_att("sam_kmin");
   if (! sam_kmin)
     return false;
   Nc3Att *sam_kmax = dataFile.get_att("sam_kmax");
   if (! sam_kmax)
     return false;
   Nc3Att *sam_kincr = dataFile.get_att("sam_kincr");
   if (! sam_kincr)
     return false;

   // long ntime = timeDim->size();
   long nz0 = z0Dim->size();
   long ny0 = y0Dim->size();
   long nx0 = x0Dim->size();

   Nc3Var *z0 = dataFile.get_var("z0");
   if (! z0)
     return false;

   Nc3Var *y0 = dataFile.get_var("y0");
   if (! y0)
     return false;
   Nc3Var *x0 = dataFile.get_var("x0");
   if (! x0)
     return false;

   double *xs = new double[nx0];
   double *ys = new double[ny0];
   double *zs = new double[nz0];

   bool success = true;
   if ( (xs == NULL) || (ys == NULL) || (zs == NULL) ) {
     std::cout << "Failed to allocate memory for Fractl grid coordinates." << std::endl;
     success = false;
   }

   if (success && ! z0->get(zs, nz0)) {
     std::cout << "Failed to read z scale." << std::endl;
     success = false;
   }
   if (success && ! y0->get(ys, ny0)) {
     std::cout << "Failed to read y scale." << std::endl;
     success = false;
   }
   if (success && ! x0->get(xs, nx0)) {
     std::cout << "Failed to read x scale." << std::endl;
     success = false;
   }

   if (success) {
     idim = sam_idim->as_long(0);
     jdim = sam_jdim->as_long(0);
     kdim = sam_kdim->as_long(0);

     imin  = sam_imin->as_double(0);
     imax  = sam_imax->as_double(0);
     iincr = sam_iincr->as_double(0);

     jmin  = sam_jmin->as_double(0);
     jmax  = sam_jmax->as_double(0);
     jincr = sam_jincr->as_double(0);

     kmin  = sam_kmin->as_double(0);
     kmax  = sam_kmax->as_double(0);
     kincr = sam_kincr->as_double(0);

     // TODO Do we need to read in the sigmas?

     // CostFunction3D reads grid dims from the configHash

     configHash["i_min"]  = std::to_string(imin);
     configHash["i_max"]  = std::to_string(imax);
     configHash["i_incr"] = std::to_string(iincr);
     configHash["j_min"]  = std::to_string(jmin);
     configHash["j_max"]  = std::to_string(jmax);
     configHash["j_incr"] = std::to_string(jincr);
     configHash["k_min"]  = std::to_string(kmin);
     configHash["k_max"]  = std::to_string(kmax);
     configHash["k_incr"] = std::to_string(kincr);

     success = validateGrid();
   }

   delete[] xs;
   delete[] ys;
   delete[] zs;

   return success;
}

// Create centers at 1 second increments, from -2 seconds to +2 seconds of the computed date/time
// All centers are at the given lat and lon.

void VarDriver3D::fillRunCenters(char *cdtg, int delta, int iter, float lat, float lon)
{
  // get rid of previous centers
  clearCenters();

  // Compute ref_time
  // NCAR - question for CSU - do we deal with time zones other than UTC?  That was specified in the original file (Qt::UTC)
  std::string tString;
  datetime cDateTime = ParseDate(cdtg, "%Y%m%d%H") + std::chrono::seconds(delta * iter);

  configHash.insert("ref_time", PrintTime(cDateTime));
  configHash.insert("ref_lat",  std::to_string(lat));
  configHash.insert("ref_lat",  std::to_string(lon));

  // create 6 "centers" centered around the ref time, at 1 second intervals

  float zero = 0.0;	// Framecenter constructor needs a reference... why?
  for(int delta2 = -2; delta2 <= 2; delta2 += 1) {
    datetime tTime = ParseDate(cdtg, "%Y%m%d%H") + std::chrono::seconds(delta * iter) + std::chrono::seconds(delta2);
    frameVector.push_back(FrameCenter(tTime, lat, lon, zero, zero));
  }
}

// fill up variables from the config hash  and call the common validateGrid()

bool VarDriver3D::validateFixedGrid()
{
  // Define the grid dimensions
  imin = std::stof(configHash["i_min"]);
  imax = std::stof(configHash["i_max"]);
  iincr = std::stof(configHash["i_incr"]);
  idim = (int)((imax - imin) / iincr) + 1;

  jmin = std::stof(configHash["j_min"]);
  jmax = std::stof(configHash["j_max"]);
  jincr = std::stof(configHash["j_incr"]);
  jdim = (int)((jmax - jmin) / jincr) + 1;

  kmin = std::stof(configHash["k_min"]);
  kmax = std::stof(configHash["k_max"]);
  kincr = std::stof(configHash["k_incr"]);
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
  // Read in the meteorological observations, process them into weights and positions
  // Either preprocess from raw observations or load an already processed Observations.in file

  std::string preprocess = configHash["preprocess_obs"];
  if (preprocess == "true") { // it should be true, testing
    bgWeights = new real[uStateSize];
    bool success = preProcessMetObs();
    delete[] bgWeights;

    if (! success) {
      cout << "Error pre-processing observations\n";
      return false;
    }

  } else {
    if (!loadPreProcessMetObs()) {
      cout << "Error loading observations\n";
      return false;
    }
  }

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
		if (std::stof(configHash["output_pressure_increment"]) > 0) {
      obCost3D = new CostFunctionXYP(projection, obVector.size(), bStateSize);
    } else if (configHash["output_COAMPS"] == "true") {
      CostFunctionCOAMPS *cf = new CostFunctionCOAMPS(projection, obVector.size(), bStateSize);
      cf->setSigmas(sigmaTable, kdim); // TODO kdim (grid) vs. number of sigmas. Same?
      obCost3D = cf;
    } else {
      obCost3D = new CostFunctionXYZ(projection, obVector.size(), bStateSize);
    }
  } else if (runMode == RTZ) {
    obCost3D = new CostFunctionRTZ(projection, obVector.size(), bStateSize);
  }

  obCost3D->initialize(&configHash, bgU, obs, refstate);
  return true;
}

void VarDriver3D::dumpBgu()
{
  std::cout << "------------- Dump of bgU after initialization ---------------" << std::endl;

  for(int64_t i = 0; i < uStateSize; i++) {
    if( (i % 20) == 0)
      std::cout << std::endl;
    std::cout << bgU[i] << " ";
  }
  std::cout << std::endl << "------ Done with bgU dump -----";
}

void VarDriver3D::dumpBgIn()
{
 std::cout << "------------- Dump of bgIn after initialization ---------------" << std::endl;
 for(int i = 0; i < bgIn.size(); i++) {
   if( (i % 11) == 0)
     std::cout << std::endl;
   std::cout << bgIn[i] << " ";
 }
 std::cout << std::endl << "------ Done with bgIn dump -----";
}
