#include <cmath>
//#include <netcdfcpp.h>
#include <Ncxx/Nc3File.hh>

#include "BkgdObsLoaders.h"
#include "BSpline.h"
#include "timers.h"
#include "datetime.h"

// ---------------- The Background Obs Loader factory --------------------

BkgdObsLoader *BkgdObsLoaderFactory::createBkgdObsLoader(BkgdObsLoader::bg_loader_t t)
{
  switch(t) {
  case BkgdObsLoader::BG_LOADER_PRE:
    return new BkgdObsPreLoader();
  case BkgdObsLoader::BG_LOADER_SPLINE:
    return new BkgdObsSplineLoader();  //Default should be: BkgdObsSplineLoader, If I want to load WRF obs, I need to change it to BkgdObsWRFLoader
  case BkgdObsLoader::BG_LOADER_KD:
    return new BkgdObsKDLoader();
  case BkgdObsLoader::BG_LOADER_FRACTL:
    return new BkgdObsFractlLoader();
	// case BkgdObsLoader::BG_LOADER_WRF:
	//   return new BkgdObsWRFLoader();
  default:
    std::cerr << "Unsupported BkgdObsLoader type in BkgdObsLoaderFactory::createBkgdObsLoader." << std::endl;
    return NULL;
  }
}

// ---------------- The base BkgdObsLoader --------------------

BkgdObsLoader::BkgdObsLoader()
{
}

BkgdObsLoader::~BkgdObsLoader()
{
}

// Set various variables used by member functions

bool BkgdObsLoader::initialize(HashMap *config,
			       std::vector<FrameCenter> frames,
			       BkgdAdapter *adapter,
			       Projection proj,
			       ReferenceState *refSt,
			       int64_t uStSize,
			       real *bgu,
			       unsigned int varsNum,
			       int id, int jd, int kd,
			       float imn, float jmn, float kmn,
			       float imx, float jmx, float kmx,
			       real iinc, real jinc, real kinc)
{
  configHash = config;
  frameVector = frames;
  bkgdAdapter = adapter;
  projection = proj;
  refstate = refSt;
  uStateSize = uStSize;
  bgWeights =  NULL;
  bgU = bgu;
  numVars = varsNum;

  idim = id;
  jdim = jd;
  kdim = kd;

  iincr = iinc;
  jincr = jinc;
  kincr = kinc;

  imin = imn;
  jmin = jmn;
  kmin = kmn;

  imax = imx;
  jmax = jmx;
  kmax = kmx;

  // Validate the run geometry. Vardriver should have already checked if the value is valid
  // TODO: DRY in VarDriver3D

  if ((*configHash)["mode"] == "XYZ") {
    runMode = RUN_MODE_XYZ;
  } else if ((*configHash)["mode"] == "RTZ") {
    runMode = RUN_MODE_RTZ;
  }
  interp_mode = (*configHash)["bg_interpolation"];

  return true;
}

// Attempt to fill holes everywhere after interpolation

bool BkgdObsLoader::fillHoles(std::vector<int> &emptybg)
{
  if (emptybg.size() > 0) {
    // Attempt to fill holes everywhere
    std::cout << "** trying to fill in " << emptybg.size() << " holes in bgU" << std::endl;

    int neighbors[6];
    real avg[numVars];

    for (std::size_t i = 0; i < emptybg.size(); i++) {
      for (unsigned int var = 0; var < numVars; var++)
	avg[var] = 0.0;
      // Check the neighbors
      int bIndex = emptybg[i];
      neighbors[0] = bIndex + numVars * (idim + 1) * 2 * (jdim + 1) * 2;
      neighbors[1] = bIndex - numVars * (idim + 1) * 2 * (jdim + 1) * 2;
      neighbors[2] = bIndex + numVars * (idim + 1) * 2;
      neighbors[3] = bIndex - numVars * (idim + 1) * 2;
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
	for (unsigned int var = 0; var < numVars; var++)
	  bgU[bIndex + var] = avg[var] / count;
      } else {
	if ((*configHash)["allow_background_missing_values"] != "true") {
	  std::cout << "Too large a hole in the background field!\n";
	  std::cout << "Please check your background file or ROI values.\n";
	  std::cout << "If you want to allow missing (zero) values, add the following line to <options>:\n";
	  std::cout << "<allow_background_missing_values>true</allow_background_missing_values>\n";
	  exit(-1);
	}
      }
    }
  }
  return true;
}

// Check that the given time is within the given time frame

bool BkgdObsLoader::timeCheck(real time, datetime&startTime, datetime &endTime, int &tci)
{
  std::string bgTimestring, tcstart, tcend;
  datetime bgTime;

  // NCAR note: I'm still confused about times - sometimes they're reals, sometimes ints??  Treating as seconds here
  // bgTime.setTime_t(time);
  // bgTime.setTimeSpec(Qt::UTC);
  unsigned int itime = (unsigned int)(time);
  //std::cout << "iTime: " << itime << std::endl;
  bgTime = std::chrono::time_point<std::chrono::system_clock>(std::chrono::duration<unsigned int>((unsigned int)(time)));

  // bgTimestring = bgTime.toString(Qt::ISODate);
  // tcstart = startTime.toString(Qt::ISODate);
  // tcend = endTime.toString(Qt::ISODate);
  bgTimestring = PrintDate(bgTime);
	tcstart = PrintDate(startTime);
  tcend = PrintDate(endTime);

  if ((bgTime < startTime) or (bgTime > endTime)) {
    std::cout << "TIME: " << time << std::endl;
    std::cout << std::endl << "time: " << bgTimestring
	      << ", start: " << tcstart
	      << ", end: " << tcend
	      << std::endl;
    exit(1);
    return false;
  }

  //tci = startTime.secsTo(bgTime);
  tci = std::chrono::duration_cast<std::chrono::seconds>(bgTime - startTime).count();
  if ((tci < 0) or (tci > (int)frameVector.size())) {
    std::cout << "Time problem with observation " << tci << " secs more than center entries" << std::endl;
    return false;
  }

  return true;
}

// This uses the kd tree implementation in the lrose-core distribution
// The returned tree has no observation data. It is built only with observation coordinates
//

KD_tree *BkgdObsLoader::buildKDTree(std::vector<double> &bgIn)
{
  long numPoints;
  int numDims = 3;

  bool debug = isTrue("debug_kd_build");

  if (debug)
    std::cout << "--- start of kdtree matrix (bgX, bgY, bgZ from bgIn. using exp(logZ). " << std::endl;

  numPoints = bgIn.size() / 11;			// 11 entries per observations
  KD_real **matrix = new KD_real*[numPoints];	// numPoints * ndim

  long bgIndex = 0;
  for (long i = 0; i < numPoints; i++) {
    matrix[i] = new KD_real[numDims];

    matrix[i][0] = bgIn[bgIndex + 0];		// bgX
    matrix[i][1] = bgIn[bgIndex + 1];		// bgY
    matrix[i][2] = exp(bgIn[bgIndex + 2]);	// exp(logZ)

    if (debug)
	std::cout <<   "bgX: " <<  matrix[i][0]
		  << ", bgY: " <<  matrix[i][1]
		  << ", bgZ: " <<  matrix[i][2]
		  << std::endl;

    bgIndex += 11;
  }
  if (debug)
    std::cout << "--- end of kdtree matrix" << std::endl;;

  return new KD_tree( (const KD_real **) matrix,
		       numPoints,
		       numDims);
}

void BkgdObsLoader::dumpBgIn(int from, int to, std::vector<real> &bgIn)
{
  std::cout << std::endl
	    << "---------------------- bgIn[" << from << ":" << to << "] ----------------"
	    << std::endl << std::endl;
  for(int i = from; i < to; i++) {
    int index = i * 11;
    std::cout << i << ": x: " << bgIn[index] << ", y: " << bgIn[index + 1] << ", z: " << bgIn[index + 2]
	      << ", time: " << bgIn[index + 3] << ", rhou: " << bgIn[index + 4]
	      << ", rhov: " << bgIn[index + 5] << ", rhow: " << bgIn[index + 6]
	      << std::endl;
  }
}

void BkgdObsLoader::dumpBgU(int from, int to, real *bgu)
{
  std::cout << std::endl
	    << "---------------------- start of bgU[" << from << ":" << to << "] ----------------"
	    << std::endl << std::endl;
  for(int i = from; i < to; i++) {
    int index = i * 7;
    if (index >= uStateSize)
      break;
    std::cout << index << ": U: " << bgu[index] << ", V: " << bgu[index + 1] << ", W: " << bgu[index + 2]
	      << ", tprime: " << bgu[index + 3] << ", qvPrime: " << bgu[index + 4]
	      << ", rhoPrime: " << bgu[index + 5] << ", qr: " << bgu[index + 6]
	      << std::endl;
  }
  std::cout << "---------------------- end of bgU[" << from << ":" << to << "] ----------------" << std::endl;
}

void BkgdObsLoader::bgu2nc(const char *fname, real *bgu)
{

  std::cout << "--- Dumping bgU in " << fname << std::endl;

  Nc3Error err(Nc3Error::verbose_nonfatal);
  // int NC_ERR = 0;

  // Create the file.
  Nc3File dataFile(fname, Nc3File::Replace);

  // Check to see if the file was created.
  if(!dataFile.is_valid()) {
    std::cout << "Failed to open " << fname << std::endl;
    return;
  }

  // Create a 3D array to be written out

  int iDim = (idim + 1) * 2;
  int jDim = (jdim + 1) * 2;
  int kDim = (kdim + 1) * 2;

  std::cout << "idim: " << idim << ", jdim: " << jdim << ", kdim: " << kdim << std::endl;
  std::cout << "iDim: " << iDim << ", jDim: " << jDim << ", kDim: " << kDim << std::endl;


  double *data = (double *) malloc(iDim * jDim * kDim * sizeof(double));

  if (data == NULL) {
    std::cout << "BkgdObsLoader::dumpBguAsNCDF: Failed to allocate data array"
	      << std::endl;
    return;
  }

  Nc3Dim *dim1 = dataFile.add_dim("IDIM", iDim);
  Nc3Dim *dim2 = dataFile.add_dim("JDIM", jDim);
  Nc3Dim *dim3 = dataFile.add_dim("KDIM", kDim);

  if ( (! dim1) || (! dim2) || (! dim3) ) {
    delete data;
    std::cout << "Failed to add dim(s)" << std::endl;
    return;
  }

  string var_names[] = { "rhou", "rhov", "rhow", "tprime", "qvprime", "rhoprime", "logZ" };

  // For each variable we want to dump

  for (unsigned int var = 0; var < numVars; var++) {

    for (int ii = -1; ii < (idim); ii++) {
      for (int imu = -1; imu <= 1; imu += 2) {

	for (int ji = -1; ji < (jdim); ji++) {
	  for (int jmu = -1; jmu <= 1; jmu += 2) {

	    for (int ki = -1; ki < (kdim); ki++) {
	      for (int kmu = -1; kmu <= 1; kmu += 2) {

		int bgI = (ii + 1) * 2 + (imu + 1) / 2;
		int bgJ = (ji + 1) * 2 + (jmu + 1) / 2;
		int bgK = (ki + 1) * 2 + (kmu + 1) / 2;

		// index into bgU (flat array)
		int bIndex = numVars * (idim + 1) * 2 * (jdim + 1) * 2 * bgK
		  + numVars * (idim + 1) * 2 * bgJ + numVars * bgI;

		// Index into data array (also a flat array)
		// int dIndex = bgK + kDim * (bgJ + jDim * bgI);

		int dIndex = bgI + iDim * (bgJ + jDim * bgK);   // lets reverse it for Ncdf

		data[dIndex] = bgU[bIndex + var];

		// std::cout << "bIndex: " << bIndex << ", dIndex: " << dIndex
		// 	  << ", val: " << bgU[bIndex]
		// 	  << ", bgI: " << bgI << ", bgJ: " << bgJ << ", bgK: " << bgK
		// 	  << std::endl;

	      }
	    }
	  }
	}
      }
    }

    Nc3Var *the_var = dataFile.add_var(var_names[var].c_str(), nc3Double, dim3, dim2, dim1);
    if ( ! the_var) {
      std::cout << "Failed to add_var" << std::endl;
      delete data;
      return;
    }

    for(long rec = 0; rec < kDim; rec++)
      if ( ! the_var->put_rec(&data[iDim * jDim * rec], rec) )
	std::cout << "Failed to write data at " << rec << " elevation" << std::endl;
  }
  delete data;
}


// ---------------- The Preprocessed Background Obs Loader --------------------

BkgdObsPreLoader::BkgdObsPreLoader()
{
}

BkgdObsPreLoader::~BkgdObsPreLoader()
{
}

bool BkgdObsPreLoader::loadBkgdObs(std::vector<real> &bgIn) {
  return false;
}

// ---------------- The Spline Background Obs Loader --------------------

BkgdObsSplineLoader::BkgdObsSplineLoader()
{

}

BkgdObsSplineLoader::~BkgdObsSplineLoader()
{
  delete[] bgWeights;
}

bool BkgdObsSplineLoader::initialize(HashMap *config,
				     std::vector<FrameCenter> frames,
				     BkgdAdapter *adapter,
				     Projection proj,
				     ReferenceState *refSt,
				     int64_t uStSize,
				     real *bgu,
				     unsigned int varsNum,
				     int id, int jd, int kd,
				     float imn, float jmn, float kmn,
				     float imx, float jmx, float kmx,
				     real iinc, real jinc, real kinc)
{
  // Call initialize from the base class to set all the variables

  BkgdObsLoader::initialize(config,
			    frames,
			    adapter,
			    proj,
			    refSt,
			    uStSize,
			    bgu,
			    varsNum,
			    id,  jd, kd,
			    imn, jmn, kmn,
			    imx, jmx, kmx,
			    iinc, jinc,kinc);

  // Now we can set Spline specific stuff

  if(uStateSize <= 0) {
    std::cerr << "BkgdObsSplineLoader::initialize: uStateSize is 0" << std::endl;
    return false;
  }
  bgWeights = new real[uStateSize];
  if (bgWeights == NULL) {
    std::cerr << "BkgdObsSplineLoader::initialize: failed to allocate bgWeights" << std::endl;
    return false;
  }

  for (int i = 0; i < uStateSize; i++) {
    bgU[i] = 0.0;
    bgWeights[i] = 0.0;
  }

  return true;
}

// Fill bgIn with observations within time and domain
// bgIn has 11 fields (see bgIn <<)
//
// Fill bgU with gridded values extrapolated from bgIn

bool BkgdObsSplineLoader::loadBkgdObs(std::vector<real> &bgIn)
{
  // Turn Debug on if there are problems with the vertical spline interpolation,
  // Eventually this should be replaced with the internal spline code
  // SplineD::Debug(1);

  // Geographic functions
  // GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();
  real referenceLon = std::stof((*configHash)["ref_lon"]);

  std::string time;
  real lat, lon, alt, u, v, w, t, qv, rhoa, qr, terrain_hgt;

  real Pi = acos(-1);

  bool debugDomain = isTrue("debug_domain");

  // bgZ = -32768.;  TODO. What was that about?

  // background is in km, ROI is gridpoints

  iROI = std::stof((*configHash)["i_background_roi"]) / iincr;
  jROI = std::stof((*configHash)["j_background_roi"]) / jincr;
  maxGridDist = 3.0;

  std::string interp_mode = (*configHash)["bg_interpolation"];
  if (interp_mode == "Cressman")
    maxGridDist = 1.0;

  if (frameVector.size() == 0) {
    std::cout << "Frame Vector is not initialized."  << std::endl;
    return false;
  }

  std::cout << "Loading background onto Gaussian mish with " << iROI
	    << " grid length radius of influence in i direction" << std::endl;
  std::cout << "and " << jROI << " grid length radius of influence in j direction" << std::endl;

  datetime startTime = frameVector.front().getTime();
  datetime endTime = frameVector.back().getTime();

  std::cout << "Start time: " << Date(startTime) << ",  " << PrintDate(startTime)  << std::endl;
  std::cout << "End  time:  " << Date(endTime) << ", " << PrintDate(endTime)  << std::endl;

  int timeProblem = 0;
  int domainProblem = 0;
  int radiusProblem = 0;
  int levelProblem = 0;
  int splineProblem = 0;

  std::cout << "iROI: " << iROI << ", jROI: " << jROI << ", maxGridDist: "
	    << maxGridDist << std::endl;
  std::cout << "imin: " << imin << ", iincr: " << iincr << std::endl;
  std::cout << "jmin: " << jmin << ", jincr: " << jincr << std::endl;
  std::cout << "kmin: " << kmin << ", kincr: " << kincr << std::endl;

  while( bkgdAdapter->next(time, lat, lon, alt, u, v, w, t, qv, rhoa, qr, terrain_hgt) ) {
    int tci;
		datetime dtime = ParseDate(time, "%Y-%m-%d_%H:%M:%S");
		int it = Date(dtime);
    if (! timeCheck(it, startTime, endTime, tci)) {
      timeProblem++;
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
		// std::cout << "bgX = " << bgX << std::endl;
    bgY = (metY - tcY) / 1000.;

    real heightm = (alt > 10.0) ? alt : 10;	// don't want to take log of zero below...
		// std::cout << "heightm = " << heightm << std::endl;
    bgZ = heightm / 1000.;
    bgRadius = sqrt(bgX * bgX + bgY * bgY);
    bgTheta = 180.0 * atan2(bgY, bgX) / Pi;
    if ((*configHash)["allow_negative_angles"] != "true")
      if (bgTheta < 0)
				bgTheta += 360.0;


    // Make sure the ob is in the Interpolation domain

    if (runMode == RUN_MODE_XYZ) {
      if ((bgX < (imin - iincr - (iROI * iincr * maxGridDist))) or
	  (bgX > (imax + iincr + (iROI * iincr * maxGridDist))) or
	  (bgY < (jmin - jincr - (jROI * jincr * maxGridDist))) or
	  (bgY > (jmax + jincr + (jROI * jincr * maxGridDist))) or
	  (bgZ < kmin)) { //Allow for higher values for interpolation purposes
	domainProblem++;

	if (debugDomain)
	  std::cout << "bgX: " << bgX << " bgY: " << bgY << " bgZ: " << bgZ << std::endl;
	continue;
      }
    } else if (runMode == RUN_MODE_RTZ) {
      if ((bgRadius < (imin - iincr - (iROI * iincr * maxGridDist))) or
	  (bgRadius > (imax + iincr + (iROI * iincr * maxGridDist))) or
	  (bgTheta < jmin - jincr - (jROI * jincr * maxGridDist)) or
	  (bgTheta > jmax + jincr + (jROI * jincr * maxGridDist)) or
	  (bgZ < kmin)) { //Exceeding the Theta domain only makes sense for sectors
	domainProblem++;
	continue;
      }
      real cylUm = (Um * bgX + Vm * bgY) / bgRadius;
      real cylVm = (-Um * bgY + Vm * bgX) / bgRadius;
      Um = cylUm;
      Vm = cylVm;
    }

    // Reference states
    real rhoBar = refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm);
    real qBar = refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm);
    real tBar = refstate->getReferenceVariable(ReferenceVariable::tempref, heightm);

    real rho = rhoa + rhoa * qv / 1000.;
    real rhou = rho * (u - Um);
    real rhov = rho * (v - Vm);
    real rhow = rho * w;
    real tprime = t - tBar;
    qv = refstate->bhypTransform(qv);
    real qvprime = qv - qBar;
    real rhoprime = (rhoa - rhoBar) * 100;
    real logZ = log(bgZ);

    if ((*configHash)["qr_variable"] == "qr")
      qr = refstate->bhypTransform(qr);

    // We assume here that the background precipitation field is always zero
    // real qr = 0.;

    if (runMode == RUN_MODE_XYZ) {
			bgIn.push_back(bgX);
			bgIn.push_back(bgY);
			bgIn.push_back(logZ);
			bgIn.push_back(it);
			bgIn.push_back(rhou);
			bgIn.push_back(rhov);
			bgIn.push_back(rhow);
			bgIn.push_back(tprime);
			bgIn.push_back(qvprime);
			bgIn.push_back(rhoprime);
			bgIn.push_back(qr);
    } else if (runMode == RUN_MODE_RTZ) {
			bgIn.push_back(bgRadius);
			bgIn.push_back(bgTheta);
			bgIn.push_back(logZ);
			bgIn.push_back(it);
			bgIn.push_back(rhou);
			bgIn.push_back(rhov);
			bgIn.push_back(rhow);
			bgIn.push_back(tprime);
			bgIn.push_back(qvprime);
			bgIn.push_back(rhoprime);
			bgIn.push_back(qr);
		}

    // Push values for an entire column, then solve for the spline
    if ( (logheights.size() != 0) && (logZ <= logheights.back()) ) {
      // Not first level, not same column, Solve for the spline
      if (logheights.size() == 1) {
	std::cerr << "Error at " << lat << ", " << lon << ", " << bgZ << std::endl
		  << "Only one level found in background spline setup. " << std::endl
		  << "Please check Background.in to ensure sorting by descending Z coordinate and re-run."
		  << std::endl;
	levelProblem++;
	return false;
      }
      splineSolver(0);	// This will also clear logheights, uBG, vBG, ...
    }

    logheights.push_back(logZ);
    uBG.push_back(rhou);
    vBG.push_back(rhov);
    wBG.push_back(rhow);
    tBG.push_back(tprime);
    qBG.push_back(qvprime);
    rBG.push_back(rhoprime);
    zBG.push_back(qr);

  } // each obs
	// std::cout << "height size = " << logheights.size() << std::endl;
  std::cout << "timeProblem: " << timeProblem << ", domainProblem: " << domainProblem
       << ", radiusProblem: " << radiusProblem << ", levelProblem: " << levelProblem
       << ", splineProblem: " << splineProblem
       << std::endl;

  if (!logheights.size()) {
    // Error reading in the background field
    std::cout << "No background estimates read in. "
	 << "Please check the time and location of your background field.\n";
#if 0
    std::cout << "Observation window: " << tcstart.toStdString() << " to " << tcend.toStdString() << "\n";
    std::cout << "Background time: " << bgTimestring.toStdString() << "\n";
#endif
    return false;
  }

  // std::cout << "------------------------2\n";

  // Solve for the last spline
  if (logheights.size() == 1) {
    std::cerr << "Only one level found in background spline setup. "
	 << "Please check Background.in to ensure sorting by Z coordinate and re-run." << std::endl;
    return false;
  }

  // Exponential interpolation in horizontal, b-Spline interpolation on log height in vertical

  splineSolver(2);  // This will also clear logheights, uBG, vBG, ...

  int numbgObs = bgIn.size() / 11; // 11 entries per ob
  if (isTrue("debug_bgIn"))
    dumpBgIn(0, numbgObs, bgIn);

  std::vector<int> emptybg;

  // TODO Do we need to check interpolation and fill hole for KD too?
  // bgWeight is a byproduct of calling the spline solver, so that's not available.
  // what about fill holes? It depends on emptyBg which is also a byproduct of the spline solver.

  // Check interpolation

  if (numbgObs > 0) {

    START_TIMER(timeci);

    for (int ki = -1; ki < (kdim); ki++) {
      for (int kmu = -1; kmu <= 1; kmu += 2) {
	real kPos = kmin + kincr * (ki + (0.5 * sqrt(1. / 3.) * kmu + 0.5));

	for (int ii = -1; ii < (idim); ii++) {
	  for (int imu = -1; imu <= 1; imu += 2) {
	    real iPos = imin + iincr * (ii + (0.5 * sqrt(1. / 3.) * imu + 0.5));

	    for (int ji = -1; ji < (jdim); ji++) {
	      for (int jmu = -1; jmu <= 1; jmu += 2) {
		real jPos = jmin + jincr * (ji + (0.5 * sqrt(1. / 3.) * jmu + 0.5));

		int bgI = (ii + 1) * 2 + (imu + 1) / 2;
		int bgJ = (ji + 1) * 2 + (jmu + 1) / 2;
		int bgK = (ki + 1) * 2 + (kmu + 1) / 2;

		int bIndex = numVars * (idim + 1) * 2 * (jdim + 1) * 2 * bgK
		  + numVars * (idim + 1) * 2 * bgJ + numVars * bgI;

		for (unsigned int var = 0; var < numVars; var++) {
		  if (bgWeights[bIndex] != 0) {
		    bgU[bIndex + var] /= bgWeights[bIndex];
		  } else {
		    emptybg.push_back(bIndex);
		    if (emptybg.size() < 15) {
		      std::cout << "Empty background mish for variable "
				<< var << " at " << iPos << ", " << jPos << ", " << kPos << std::endl;
		    } else if (emptybg.size() == 15) {
		      std::cout << "Too many empty mish points, will no longer report.\n";
		    }
		  }
		}
		if ((*configHash)["qr_variable"] == "dbz") {
		  if (bgU[bIndex + 6] > 0) {
		    real dbzavg = 10 * log10(bgU[bIndex + 6]);
		    bgU[bIndex + 6] = (dbzavg + 35.) * 0.1;
		  } else {
		    bgU[bIndex + 6] = 0.0;
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    fillHoles(emptybg);

    PRINT_TIMER("Interpolation check", timeci);
  } else {
    std::cout << "No background observations loaded" << std::endl;
    return false;
  }

  std::cout << numbgObs << " background observations loaded" << std::endl;

  if (isTrue("debug_bgu")) {
    dumpBgU(0, uStateSize, bgU);
  }

  if (configHash->exists("debug_bgu_nc"))
    bgu2nc((*configHash)["debug_bgu_nc"].c_str(), bgU);

  return true;
}

// Exponential interpolation in horizontal,
// b-Spline interpolation on log height in vertical

// Pre:
//   logheights, uBG, vBG, ..., zBG (contain all the background fields for a column)
// Post:
//   some bgU values get updated

bool BkgdObsSplineLoader::splineSolver(int waveLen)
{
  real Rsquare = (iincr * iROI) * (iincr * iROI) + (jincr * jROI) * (jincr * jROI);

  int radiusProblem = 0;
  int splineProblem = 0;

#pragma omp parallel for

  for (int ki = -1; ki < kdim; ki++) {

    SplineD* bgSpline = new SplineD(&logheights.front(), logheights.size(),
				    uBG.data(), waveLen, SplineBase::BC_ZERO_SECOND);
    if (!bgSpline->ok()) {
      std::cerr << "bgSpline setup failed." << std::endl;
      splineProblem++;
      continue; //return -1;
    }

    for (int kmu = -1; kmu <= 1; kmu += 2) {
      real kPos = kmin + kincr * (ki + (0.5 * sqrt(1. / 3.) * kmu + 0.5));
      if (kPos < 0) kPos = 0.001;
      real logzPos = log(kPos);
      if (logzPos < logheights[0]) logzPos = logheights[0];
      //if (fabs(kPos-obZ) > kincr*ROI*2.) continue;

      for (int ii = -1; ii < idim; ii++) {

	for (int imu = -1; imu <= 1; imu += 2) {
	  real iPos = imin + iincr * (ii + (0.5 * sqrt(1. / 3.) * imu + 0.5));
	  if (runMode == RUN_MODE_XYZ) {
	    if (fabs(iPos - bgX) > iincr * iROI  * maxGridDist) {
	      ++radiusProblem;
	      continue;
	    }
	  } else if (runMode == RUN_MODE_RTZ) {
	    if (fabs(iPos - bgRadius) > iincr * iROI * maxGridDist) {
	      ++radiusProblem;
	      continue;
	    }
	  }

	  for (int ji = -1; ji < jdim; ji++) {

	    for (int jmu = -1; jmu <= 1; jmu += 2) {
	      real jPos = jmin + jincr * (ji + (0.5 * sqrt(1. / 3.) * jmu + 0.5));
	      real rSquare = 0.0;
	      if (runMode == RUN_MODE_XYZ) {
		if (fabs(jPos - bgY) > jincr * jROI * maxGridDist) {
		  ++radiusProblem;
		  continue;
		}
		rSquare = (bgX - iPos) * (bgX - iPos) + (bgY - jPos) * (bgY - jPos);
	      } else if (runMode == RUN_MODE_RTZ) {
		real dTheta = fabs(jPos - bgTheta);
		if (dTheta > 360.) dTheta -= 360.;
		if (dTheta > jincr * jROI*2.) {
		  ++radiusProblem;
		  continue;
		}
		rSquare = (bgRadius - iPos) * (bgRadius - iPos) + (dTheta) * (dTheta);
	      }

	      // Add one extra index to account for buffer zone in analysis
	      int bgI = (ii + 1) * 2 + (imu + 1) / 2;
	      int bgJ = (ji + 1) * 2 + (jmu +1 ) / 2;
	      int bgK = (ki + 1) * 2 + (kmu + 1) / 2;

	      int bIndex = numVars * (idim + 1) * 2 * (jdim + 1) * 2 * bgK
		+ numVars * (idim + 1) * 2 * bgJ + numVars * bgI;

	      if (rSquare < Rsquare * maxGridDist) {
		real weight = exp(-2.302585092994045 * rSquare / Rsquare);
		if (interp_mode == "Cressman") {
		  weight = (Rsquare - rSquare) / (Rsquare + rSquare);
		}
		if (logzPos > logheights.front()) {
		  bgSpline->solve(uBG.data());
		  bgU[bIndex] += weight * (bgSpline->evaluate(logzPos));
		  bgSpline->solve(vBG.data());
		  bgU[bIndex + 1] += weight * (bgSpline->evaluate(logzPos));
		  bgSpline->solve(wBG.data());
		  bgU[bIndex + 2] += weight * (bgSpline->evaluate(logzPos));
		  bgSpline->solve(tBG.data());
		  bgU[bIndex + 3] += weight * (bgSpline->evaluate(logzPos));
		  bgSpline->solve(qBG.data());
		  bgU[bIndex + 4] += weight * (bgSpline->evaluate(logzPos));
		  bgSpline->solve(rBG.data());
		  bgU[bIndex + 5] += weight * (bgSpline->evaluate(logzPos));
		  bgSpline->solve(zBG.data());
		  bgU[bIndex + 6] += weight * (bgSpline->evaluate(logzPos));
		  bgWeights[bIndex] += weight;
		} else {
		  // Below the spline interpolation
		  bgU[bIndex]     += weight * uBG.front();
		  bgU[bIndex + 1] += weight * vBG.front();
		  bgU[bIndex + 2] += weight * wBG.front();
		  bgU[bIndex + 3] += weight * tBG.front();
		  bgU[bIndex + 4] += weight * qBG.front();
		  bgU[bIndex + 5] += weight * rBG.front();
		  bgU[bIndex + 6] += weight * zBG.front();
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

  return true;
}

// ---------------- The K-D Tree Background Obs Loader ------------------

BkgdObsKDLoader::BkgdObsKDLoader() {}
BkgdObsKDLoader::~BkgdObsKDLoader() {}

// Debug support
// Overwrite some U values in the bgU array

bool BkgdObsKDLoader::overwriteBgu(const char *fname)
{
  std::ifstream ifs(fname);
  if (! ifs.good()) {
    std::cout << "Error opening " << fname << " for reading." << std::endl;
    return false;
  }
  std::cout << "Overwriting bgU with values from " << fname << std::endl;

  int idx;
  float u, v, w, tp, qvp, rp, qr;
  std::string dummy;

  int count = 0;

  while (ifs >> idx
	>> dummy >> dummy >> u
	>> dummy >> dummy >> v
	>> dummy >> dummy >> w
	>> dummy >> dummy >> tp
	>> dummy >> dummy >> qvp
	>> dummy >> dummy >> rp
	>> dummy >> dummy >> qr	) {
    ifs.ignore(100, '\n');

    bgU[idx + 0] = u;
    bgU[idx + 1] = v;
    bgU[idx + 2] = w;
    bgU[idx + 3] = tp;
    bgU[idx + 4] = qvp;
    bgU[idx + 5] = rp;
    bgU[idx + 6] = qr;
    count++;
  }
  std::cout << "Overwrote " << count << " values in bgU" << std::endl;
  return true;
}

bool BkgdObsKDLoader::loadBkgdObs(std::vector<real> &bgIn)
{
  std::string time;
  real lat, lon, alt, u, v, w, t, qv, rhoa, qr, terrain_hgt;
  real Pi = acos(-1);

  KD_tree *kdTree;

  datetime startTime = frameVector.front().getTime();
  datetime endTime = frameVector.back().getTime();

  std::cout << "Start time: " << PrintTime(startTime) << std::endl;
  std::cout << "End  time:  " << PrintTime(endTime)  << std::endl;

  // background is in km, ROI is gridpoints
  // TODO that comment doesn't seem correct. These are set to 45.0 in the config file...

  iROI = std::stof((*configHash)["i_background_roi"]) / iincr;
  jROI = std::stof((*configHash)["j_background_roi"]) / jincr;

  maxGridDist = 3.0;

  int timeProblem = 0;
  int domainProblem = 0;

  int debugKd = 0;
  if (configHash->exists("debug_kd")) {
    debugKd = std::stoi((*configHash)["debug_kd"]);
    std::cout << "*** dbugKD: " << debugKd << std::endl;
  }

  int debugKdStep = 0;
  if (configHash->exists("debug_kd_step")) {
    debugKdStep = std::stoi((*configHash)["debug_kd_step"]);
    std::cout << "*** Debug KD Step: " << debugKdStep << std::endl;
  }

  while( bkgdAdapter->next(time, lat, lon, alt, u, v, w, t, qv, rhoa, qr, terrain_hgt) ) {
    // Process the bkgdObs into Observations
    int tci;
		datetime dtime = ParseDate(time, "%Y-%m-%d_%H:%M:%S");
		int it = Date(dtime);
    if (! timeCheck(it, startTime, endTime, tci)) {
      timeProblem++;
      continue;
    }

    real Um = frameVector[tci].getUmean();
    real Vm = frameVector[tci].getVmean();

    // Get the X, Y & Z
    real tcX, tcY, metX, metY;
    real referenceLon = std::stof((*configHash)["ref_lon"]);

    projection.Forward(referenceLon, frameVector[tci].getLat() , frameVector[tci].getLon(),
		       tcX, tcY);
    projection.Forward(referenceLon, lat, lon , metX, metY);
    bgX = (metX - tcX) / 1000.;
    bgY = (metY - tcY) / 1000.;

    real heightm = (alt > 10.0) ? alt : 10;	// don't want to take log of zero below...

    bgZ = heightm / 1000.;
    bgRadius = sqrt(bgX * bgX + bgY * bgY);
    bgTheta = 180.0 * atan2(bgY, bgX) / Pi;
    if ((*configHash)["allow_negative_angles"] != "true")
      if (bgTheta < 0)
	bgTheta += 360.0;

    // Make sure the ob is in the Interpolation domain

    if (runMode == RUN_MODE_XYZ) {
      if ((bgX < (imin - iincr - (iROI * iincr * maxGridDist))) or
	  (bgX > (imax + iincr + (iROI * iincr * maxGridDist))) or
	  (bgY < (jmin - jincr - (jROI * jincr * maxGridDist))) or
	  (bgY > (jmax + jincr + (jROI * jincr * maxGridDist))) or
	  (bgZ < kmin)) { //Allow for higher values for interpolation purposes
	domainProblem++;
	continue;
      }
    } else if (runMode == RUN_MODE_RTZ) {
      if ((bgRadius < (imin - iincr - (iROI * iincr * maxGridDist))) or
	  (bgRadius > (imax + iincr + (iROI * iincr * maxGridDist))) or
	  (bgTheta < jmin - jincr - (jROI * jincr * maxGridDist)) or
	  (bgTheta > jmax + jincr + (jROI * jincr * maxGridDist)) or
	  (bgZ < kmin)) { //Exceeding the Theta domain only makes sense for sectors
	domainProblem++;
	continue;
      }
      real cylUm = (Um * bgX + Vm * bgY) / bgRadius;
      real cylVm = (-Um * bgY + Vm * bgX) / bgRadius;
      Um = cylUm;
      Vm = cylVm;
    }

    // Reference states
    real rhoBar = refstate->getReferenceVariable(ReferenceVariable::rhoaref, heightm);
    real qBar = refstate->getReferenceVariable(ReferenceVariable::qvbhypref, heightm);
    real tBar = refstate->getReferenceVariable(ReferenceVariable::tempref, heightm);


    real rho = rhoa + rhoa * qv / 1000.;
    real rhou = rho * (u - Um);
    real rhov = rho * (v - Vm);
    real rhow = rho * w;
    real tprime = t - tBar;
    qv = refstate->bhypTransform(qv);
    real qvprime = qv - qBar;
    real rhoprime = (rhoa - rhoBar) * 100;
    real logZ = log(bgZ);

    if ((*configHash)["qr_variable"] == "qr")
      qr = refstate->bhypTransform(qr);

    // We assume here that the background precipitation field is always zero
    // real qr = 0.;

    if (runMode == RUN_MODE_XYZ) {
			bgIn.push_back(bgX);
			bgIn.push_back(bgY);
			bgIn.push_back(logZ);
			bgIn.push_back(it);
			bgIn.push_back(rhou);
			bgIn.push_back(rhov);
			bgIn.push_back(rhow);
			bgIn.push_back(tprime);
			bgIn.push_back(qvprime);
			bgIn.push_back(rhoprime);
			bgIn.push_back(qr);
    } else if (runMode == RUN_MODE_RTZ) {
			bgIn.push_back(bgRadius);
			bgIn.push_back(bgTheta);
			bgIn.push_back(logZ);
			bgIn.push_back(it);
			bgIn.push_back(rhou);
			bgIn.push_back(rhov);
			bgIn.push_back(rhow);
			bgIn.push_back(tprime);
			bgIn.push_back(qvprime);
			bgIn.push_back(rhoprime);
			bgIn.push_back(qr);
		}
  } // each obs

  std::cout << "timeProblem: " << timeProblem << ", domainProblem: " << domainProblem
            << std::endl;

  int numbgObs = bgIn.size() / 11; // 11 entries per ob

  if (numbgObs <= 0) {
    std::cout << "No background observations loaded" << std::endl;
    return false;
  }

  std::cout << numbgObs << " background observations loaded" << std::endl;

  if (isTrue("debug_bgIn"))
    dumpBgIn(0, numbgObs, bgIn);

  // All the obs are now loaded in bgIn. Build a KD tree

  kdTree = buildKDTree(bgIn);
  std::vector<int> emptybg;

  KD_real *centerLoc = new KD_real[3];	// 3 dim
  int nbrMax  = std::stoi((*configHash)["bkgd_kd_num_neighbors"]);
  if (nbrMax <= 0) {
    std::cerr << "Number of neighbors must be greater than 0." << std::endl;
    return false;
  }
  float maxDist = std::stof((*configHash)["bkgd_kd_max_distance"]);
  // KD tree returns the square of the distance...
  maxDist *= maxDist;

  // For each point on the mish

  std::cout << "*** idim: " << idim << ", jdim: " << jdim << ", kdim: " << kdim << std::endl;

  if (debugKd)
    std::cout << "--- start of debug kd" << std::endl;

  for (int ii = -1; ii < (idim); ii++) {
    for (int imu = -1; imu <= 1; imu += 2) {
      real iPos = imin + iincr * (ii + (0.5 * sqrt(1. / 3.) * imu + 0.5));

      for (int ji = -1; ji < (jdim); ji++) {
	for (int jmu = -1; jmu <= 1; jmu += 2) {
	  real jPos = jmin + jincr * (ji + (0.5 * sqrt(1. / 3.) * jmu + 0.5));

	  for (int ki = -1; ki < (kdim); ki++) {
	    for (int kmu = -1; kmu <= 1; kmu += 2) {
	      real kPos = kmin + kincr * (ki + (0.5 * sqrt(1. / 3.) * kmu + 0.5));

	      int bgI = (ii + 1) * 2 + (imu + 1) / 2;
	      int bgJ = (ji + 1) * 2 + (jmu + 1) / 2;
	      int bgK = (ki + 1) * 2 + (kmu + 1) / 2;

	      // std::cout << "iPos: " << iPos << ", jPos: " << jPos << ", kPos: " << kPos << std::endl;
	      // std::cout << "bgI: " << bgI << ", bgJ: " << bgJ << ", bgK " << bgK << std::endl;
	      // continue;

	      // index into bgU (flat array)
	      int bIndex = numVars * (idim + 1) * 2 * (jdim + 1) * 2 * bgK
		+ numVars * (idim + 1) * 2 * bgJ + numVars * bgI;

	      // Do the KD tree here.
	      // give me the n nearest neighbors and average.

	      centerLoc[0] = iPos;
	      centerLoc[1] = jPos;
	      centerLoc[2] = kPos;

	      fillBguEntry(centerLoc, nbrMax, maxDist, bgIn, kdTree, bgU, bIndex, emptybg, debugKd);
	      if (debugKd > 0)
		debugKd -= debugKdStep;
	    }
	  }
	}
      }
    }
  }

  if (emptybg.size() > 0)
    std::cout << "** KD left " << emptybg.size() << " holes in bgU" << std::endl;

  if (debugKd)
    std::cout << "--- end of debug kd" << std::endl;

  if (configHash->exists("debug_bgu_overwrite"))
    overwriteBgu((*configHash)["debug_bgu_overwrite"].c_str());

  if (isTrue("debug_bgu")) {
    dumpBgU(0, uStateSize, bgU);
  }

  if (configHash->exists("debug_bgu_nc"))
    bgu2nc((*configHash)["debug_bgu_nc"].c_str(), bgU);
  return true;
}

bool BkgdObsKDLoader::fillBguEntry(KD_real *centerLoc, int nbrMax, float maxDistance,
			     std::vector<real> &bgIn, KD_tree *kdTree, real *bgU, int bIndex,
			     std::vector<int> &emptybg,
			     int debug)
{
  int nbrIxs[nbrMax];
  std::fill(nbrIxs, nbrIxs + nbrMax, -1);
  KD_real nbrDistSqs[nbrMax];

  if (debug > 0)
    std::cout << "Point (" << centerLoc[0] << ", " << centerLoc[1] << ", " << centerLoc[2]
	      << ")" << std::endl;

  // Fill nbrIxs with indices into bgIn of the nearest neighbors
  kdTree->nnquery(centerLoc, nbrMax, KD_EUCLIDEAN, 1, nbrIxs, nbrDistSqs);

  int	numNbrActual = 0;
  float sumU = 0;
  float sumV = 0;
  float sumW = 0;
  float sumTprime = 0;
  float sumQVprime = 0;
  float sumRhoPrime = 0;
  float sumQr = 0;

  for(int inbr = 0; inbr < nbrMax; inbr++) {
    // do we have to worry about aircraft distance? (Do we even have that info?

    if (nbrIxs[inbr] < 0) {
      std::cerr << "nbrIxs < 0" << std::endl;
      continue;
    }
    if (nbrIxs[inbr] > bgIn.size() / 11) {
      std::cerr << "nbrIxs > size of bgIn" << std::endl;
      continue;
    }

    // nbrIxs[inbr] is the index of the observation. We need to convert that to
    // an index into the bgIn array

    long bgInIdx = nbrIxs[inbr] * 11;

    float bgX      = bgIn[bgInIdx + 0];
    float bgY      = bgIn[bgInIdx + 1];
    float logZ     = bgIn[bgInIdx + 2];
    float rhou     = bgIn[bgInIdx + 4];
    float rhov     = bgIn[bgInIdx + 5];
    float rhow     = bgIn[bgInIdx + 6];
    float tprime   = bgIn[bgInIdx + 7];
    float qvprime  = bgIn[bgInIdx + 8];
    float rhoprime = bgIn[bgInIdx + 9];
    float qr       = bgIn[bgInIdx + 10];

    float localDist = nbrDistSqs[inbr];

    if (debug > 0) {
      std::cout << "     bgIn Index: " << bgInIdx / 11
		<< "  neighbor " << inbr << " ("
		<< bgX << ", " << bgY << ", " << exp(logZ) << ") dist: " << sqrt(nbrDistSqs[inbr])
		<< ", U: " << rhou << ", V: " << rhov << ", W: " << rhow
		<< std::endl;
    }

    if(localDist < maxDistance) {
      // TODO: add values of fields to corresponding Stat objects?
      numNbrActual++;
      sumU += rhou;
      sumV += rhov;
      sumW += rhow;
      sumTprime += tprime;
      sumQVprime += qvprime;
      sumRhoPrime += rhoprime;
      sumQr += qr;
    } else {
      // omitted point
    }
  }

  if (numNbrActual > 0) {	// we have enough neighbors
    // set the ob mean fields to the corresponding stat / corresponding start number of good
    bgU[bIndex]		= sumU / numNbrActual;
    bgU[bIndex + 1]	= sumV / numNbrActual;
    bgU[bIndex + 2]	= sumW / numNbrActual;
    bgU[bIndex + 3]	= sumTprime / numNbrActual;
    bgU[bIndex + 4]	= sumQVprime / numNbrActual;
    bgU[bIndex + 5]	= sumRhoPrime / numNbrActual;
    bgU[bIndex + 6]	= sumQr / numNbrActual;
  } else {
    emptybg.push_back(bIndex);
  }

  return true;
}


// ---------------- The FRACTL Mish Background Obs Loader ------------------

// This loader assumes that Fractl was run with the GRID_MISH option.

BkgdObsFractlLoader::BkgdObsFractlLoader() {}
BkgdObsFractlLoader::~BkgdObsFractlLoader() {}

bool BkgdObsFractlLoader::loadBkgdObs(std::vector<real> &bgIn)
{
  Nc3Error err(Nc3Error::verbose_nonfatal);
  std::string fname = (*configHash)["fractl_nc_file"];
  // Open the file.
  Nc3File dataFile(fname.c_str(), Nc3File::ReadOnly);

  // Check to see if the file was opened.
  if(!dataFile.is_valid()) {
    std::cout << "Failed to read FRACTL generated nc file " << fname
	      << std::endl;
    return false;
  }

  // TODO The same size variables are also read in VarDriver3D::validateFractlGrid()
  //      Combine into common code?

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

  long ntime = timeDim->size();
  long nz0 = z0Dim->size();
  long ny0 = y0Dim->size();
  long nx0 = x0Dim->size();

  Nc3Var *time0 = dataFile.get_var("time");
  if (! time0)
    return false;

  Nc3Var *z0 = dataFile.get_var("z0");
  if (! z0)
    return false;

  Nc3Var *y0 = dataFile.get_var("y0");
  if (! y0)
    return false;
  Nc3Var *x0 = dataFile.get_var("x0");
  if (! x0)
    return false;
  Nc3Var *lat0 = dataFile.get_var("lat0");
  if (! lat0)
    return false;
  Nc3Var *lon0 = dataFile.get_var("lon0");
  if (! lon0)
    return false;

  Nc3Var *upward_air_velocity = dataFile.get_var("W");
  if (! upward_air_velocity)
    return false;
  Nc3Var *northward_wind = dataFile.get_var("V");
  if (! northward_wind)
    return false;
  Nc3Var *eastward_wind = dataFile.get_var("U");
  if (! eastward_wind)
    return false;
  Nc3Var *meanNbrDbz = dataFile.get_var("DBZ");
  if (! meanNbrDbz)
    return false;

  Nc3Var *meanNeighborNcp = dataFile.get_var("NCP");
  if (! meanNeighborNcp)
    return false;

  Nc3Var *W_std = dataFile.get_var("W_std");
  if (! W_std)
    return false;
  Nc3Var *V_std = dataFile.get_var("V_std");
  if (! V_std)
    return false;
  Nc3Var *U_std = dataFile.get_var("U_std");
  if (! U_std)
    return false;

  long numBgObs = ntime * nz0 * ny0 * nx0;

  if (numBgObs <= 0) {
    std::cout << "No background observations loaded" << std::endl;
    return false;
  }
  std::cout << "Loading " << numBgObs << " background observations" << std::endl;

  // Allocate data arrays

  double *times = new double[ntime];

  double *lats = new double[nx0 * ny0]; // should be ny0*nx0?
  double *lons = new double[nx0 * ny0]; // should be ny0*nx0?
  double *alts = new double[nz0];

  double *upWind = new double[numBgObs];
  double *northWind = new double[numBgObs];
  double *eastWind = new double[numBgObs];
  double *meanDbz =  new double[numBgObs];

  double *wStd = new double[numBgObs];
  double *vStd = new double[numBgObs];
  double *uStd = new double[numBgObs];

  bool success = true;

  // Make sure table allocations worked.
  // alts was allocated last so if any are going to fail that would be the one.

  if (alts == NULL) {
    std::cout << "Failed to allocate altitude table." << std::endl;
    success = false;
  }

  // Get the time array

  if ( success && ! time0->get(times, 1, ntime) ) {
    std::cout << "Failed to read in the times." << std::endl;
    success = false;
  }

  // Get the heights in meter

  if ( success && ! z0->get(alts, nz0) ) {
    std::cout << "Failed to read in altitudes." << std::endl;
    success = false;
  }

  for (int rec = 0; rec < ntime; rec++) {
    if ( ! upward_air_velocity->set_cur(rec, 0, 0, 0)) {
      std::cout << "Failed to read in upward air velocity" << std::endl;
      success = false;
    }
    if ( success && ! upward_air_velocity->get(upWind, 1, nz0, ny0, nx0)) {
      std::cout << "Failed to read in upward air velocity" << std::endl;
      success = false;
    }

    if ( success && ! northward_wind->set_cur(rec, 0, 0, 0)) {
      std::cout << "Failed to read in northward air velocity" << std::endl;
      success = false;
    }
    if ( success && ! northward_wind->get(northWind, 1, nz0, ny0, nx0)) {
      std::cout << "Failed to read in northward air velocity" << std::endl;
      success = false;
    }

    if ( success && ! eastward_wind->set_cur(rec, 0, 0, 0)) {
      std::cout << "Failed to read in eastward air velocity" << std::endl;
      success = false;
    }
    if ( success && ! eastward_wind->get(eastWind, 1, nz0, ny0, nx0)) {
      std::cout << "Failed to read in eastward air velocity" << std::endl;
      success = false;
    }

    if ( success && ! meanNbrDbz->set_cur(rec, 0, 0, 0)) {
      std::cout << "Failed to read in mean dbZ" << std::endl;
      success = false;
    }
    if ( success && ! meanNbrDbz->get(meanDbz, 1, nz0, ny0, nx0)) {
      std::cout << "Failed to read in mean dbZ" << std::endl;
      success = false;
    }

    if ( success && ! W_std->set_cur(rec, 0, 0, 0)) {
      std::cout << "Failed to read in W_std" << std::endl;
      success = false;
    }

    if ( success && ! W_std->get(wStd, 1, nz0, ny0, nx0)) {
      std::cout << "Failed to read in W_std" << std::endl;
      success = false;
    }

    if ( success && ! V_std->set_cur(rec, 0, 0, 0)) {
      std::cout << "Failed to read in V_std" << std::endl;
      success = false;
    }

    if ( success && ! V_std->get(vStd, 1, nz0, ny0, nx0)) {
      std::cout << "Failed to read in V_std" << std::endl;
      success = false;
    }

    if ( success && ! U_std->set_cur(rec, 0, 0, 0)) {
      std::cout << "Failed to read in U_std" << std::endl;
      success = false;
    }

    if ( success && ! W_std->get(uStd, 1, nz0, ny0, nx0)) {
      std::cout << "Failed to read in U_std" << std::endl;
      success = false;
    }

  }

  double fillVal = 0.0;
  // double fillVal = std::nan(""); // testing. Quick 2 iteration

  for (int t = 0; t < ntime; t++) {
    // double currentTime = times[t];

    for (int z = 0; z < nz0; z++) {

      float heightm = alts[z];	// TODO check unit
      float rho = refstate->getReferenceVariable(ReferenceVariable::rhoref, heightm);

      for (int y = 0; y < ny0; y++) {
	for (int x = 0; x < nx0; x++) {
	  int vIndex = x + nx0 * (y + ny0 * z);		// variables index
	  int bIndex = vIndex * numVars;		// index into bgU

	  double rhou = eastWind[vIndex] * rho;
	  if (std::isnan(rhou))
	    rhou = fillVal;
	  double rhov = northWind[vIndex] * rho;
	  if (std::isnan(rhov))
	    rhov = fillVal;
	  double rhow = upWind[vIndex] * rho;
	  if (std::isnan(rhow))
	    rhow = fillVal;
	  double dbz  = meanDbz[vIndex];
	  if (std::isnan(dbz))
	    dbz = -35.0;		// clear air.

	  bgU[bIndex + 0] = rhou;
	  bgU[bIndex + 1] = rhov;
	  bgU[bIndex + 2] = rhow;
	  bgU[bIndex + 3] = 0.0;			// t'
	  bgU[bIndex + 4] = 0.0;			// qv'
	  bgU[bIndex + 5] = 0.0;			// rho'
	  bgU[bIndex + 6] = (dbz + 35.) * 0.1;		// transform to scaled value [0..7]

	  // W_std, V_std, U_std are handled separatively in CostFunction3D
	}
      }
    }
  } // ntime

  if (isTrue("debug_bgu"))
    dumpBgU(0, uStateSize, bgU);

  if (configHash->exists("debug_bgu_nc"))
    bgu2nc((*configHash)["debug_bgu_nc"].c_str(), bgU);

   delete[] lats;
   delete[] lons;
   delete[] alts;
   delete[] upWind;
   delete[] northWind;
   delete[] eastWind;
   delete[] meanDbz;

   delete[] wStd;
   delete[] vStd;
   delete[] uStd;

   return success;
}

// ---------------- The WRF output Background Obs Loader ------------------

BkgdObsWRFLoader::BkgdObsWRFLoader() {}
BkgdObsWRFLoader::~BkgdObsWRFLoader() {}

bool BkgdObsWRFLoader::loadBkgdObs(std::vector<real> &bgIn)
{
		Nc3Error err(Nc3Error::verbose_nonfatal);
		double Rd = 287.058; //J/kg/K
		double Rv = 461.495; // J/kg/K
		std::string dataPath = "/Users/tingyu/devel/samurai_terrain/testDEM";
    std::string wrffile = "/Users/tingyu/Dropbox/CSU/Research/PhD/samurai/hill2d/wrfout_d01_0001-01-01_00_00_00";
		Nc3File dataFile(wrffile.c_str(), Nc3File::ReadOnly);
		if(!dataFile.is_valid()) {
			std::cout << "Failed to read WRF output netCDF file " << wrffile
					<< std::endl;
			return false;}

		Nc3Var *timeVar = dataFile.get_var("Times");
	  if (! timeVar)
	    return false;

	  Nc3Var *phVar = dataFile.get_var("PH");
	  if (! phVar)
	    return false;

		Nc3Var *phbVar = dataFile.get_var("PHB");
		if (! phbVar)
		  return false;

		Nc3Var *pVar = dataFile.get_var("P");
		if (! pVar)
			return false;

		Nc3Var *pbVar = dataFile.get_var("PB");
		if (! phbVar)
			return false;

	  Nc3Var *latVar = dataFile.get_var("XLAT");
	  if (! latVar)
	    return false;

	  Nc3Var *lonVar = dataFile.get_var("XLONG");
	  if (! lonVar)
	    return false;

	  Nc3Var *wwindVar = dataFile.get_var("W");
	  if (! wwindVar)
	    return false;

	  Nc3Var *vwindVar = dataFile.get_var("V");
	  if (! vwindVar)
	    return false;

	  Nc3Var *uwindVar = dataFile.get_var("U");
	  if (! uwindVar)
	    return false;

		Nc3Var *tempVar = dataFile.get_var("T");
		if (! tempVar)
		  return false;
		// Nc3Bool test2 = dataFile.get_var("QRAIN")->is_valid();
		Nc3Var *qrVar = dataFile.get_var("P_HYD");
		if (!(dataFile.get_var("QRAIN"))){
			Nc3Var *qrVar = dataFile.get_var("P_HYD");
			std::cout << "QRAIN Variable is not found, so QRAIN will be set to 0." << std::endl;}
		else{
			Nc3Var *qrVar = dataFile.get_var("QRAIN");}

		Nc3Var *terrainHgtVar = dataFile.get_var("HGT");
		if (! terrainHgtVar)
			return false;

	  Nc3Var *qvVar = dataFile.get_var("QVAPOR");
		if (! qvVar)
		  return false;

			//Set up dimensions
		int NREC = 5;//timeVar->get_dim(0)->size();
		long NLAT = latVar->get_dim(1)->size();
		long NLON = lonVar->get_dim(2)->size();
	  long NLVL = phVar->get_dim(1)->size()-1;

		 // Allocate data arrays
		 char *times = new char[1];
		 double *lats =  new double[NLAT*NLON];
		 double *lons = new double[NLAT*NLON];
		 double *u = new double[NLON*NLAT*NLVL];
		 double *v = new double[NLON*NLAT*NLVL];
		 double *w = new double[NLON*NLAT*NLVL];
		 double *z = new double[NLON*NLAT*NLVL];
		 double *p = new double[NLON*NLAT*NLVL];
		 double *t = new double[NLON*NLAT*NLVL];
		 double *qv = new double[NLON*NLAT*NLVL];
		 double *rhoa = new double[NLON*NLAT*NLVL];
		 double *qr = new double[NLON*NLAT*NLVL];
		 double *terrain_hgt = new double[NLON*NLAT];


		 double *uwind_in = new double[NLVL*NLAT*(NLON+1)];
		 double *vwind_in = new double[NLVL*(NLAT+1)*NLON];
		 double *wwind_in = new double[(NLVL+1)*NLAT*NLON];
		 double *ph_in = new double[(NLVL+1)*NLAT*NLON];
		 double *phb_in = new double[(NLVL+1)*NLAT*NLON];
		 double *p_in = new double[NLVL*NLAT*NLON];
		 double *pb_in = new double[NLVL*NLAT*NLON];
		 double *t_in = new double[NLVL*NLAT*NLON];

		 std::cout << "*** SUCCESS allocate data arrays" << std::endl;

		// Read the data. Since we know the contents of the file we know
	  // that the data arrays in this program are the correct size to
	  // hold one timestep.
		// I used the 5th timestep record
					if (!timeVar->set_cur(NREC)){
					std::cout << "Failed to read in timeVar" << std::endl;
				return false;}
					if (!timeVar->get(times, 1, 19, NREC)){
				return false;}
					if (!latVar->get(lats, 1, NLAT, NLON)){
				return false;}
					if (!lonVar->get(lons, 1, NLAT, NLON)){
				return false;}
					if (!terrainHgtVar->get(terrain_hgt, 1, NLAT, NLON)){
				return false;}
	  for (int rec = NREC; rec <= NREC;  rec++)
	  {
			 if (!terrainHgtVar->set_cur(rec, 0, 0)){
			 std::cout << "Failed to read in terrainHgtVar" << std::endl;
	  return false;}
	     if (!uwindVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in uwindVar" << std::endl;
	 	return false;}
	     if (!vwindVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in vwindVar" << std::endl;
	  return false;}
	     if (!wwindVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in wwindVar" << std::endl;
	  return false;}
	     if (!phVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in phVar" << std::endl;
	 	return false;}
	     if (!phbVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in phbVar" << std::endl;
		return false;}
			 if (!pVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in pVar" << std::endl;
		return false;}
			 if (!pbVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in pbVar" << std::endl;
	 	return false;}
			 if (!tempVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in tempVar" << std::endl;
	  return false;}
			 if (!qvVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in qvVar" << std::endl;
	  return false;}
			 if (!qrVar->set_cur(rec, 0, 0, 0)){
			 std::cout << "Failed to read in qrVar" << std::endl;
 		return 0;}


		// Get 1 record of NLVL by NLAT by NLON values for each variable.
		   if (!uwindVar->get(uwind_in, 1, NLVL, NLAT, NLON+1)){
			 std::cout << "Failed to get uwind_in" << std::endl;
 		return false;}
			 if (!vwindVar->get(vwind_in, 1, NLVL, NLAT+1, NLON)){
			 std::cout << "Failed to get vwind_in" << std::endl;
 		return false;}
			 if (!wwindVar->get(wwind_in, 1, NLVL+1, NLAT, NLON)){
			 std::cout << "Failed to get wwind_in" << std::endl;
 		return false;}
			 if (!phVar->get(ph_in, 1, NLVL, NLAT, NLON)){
			 std::cout << "Failed to get ph_in" << std::endl;
 		return false;}
			 if (!phbVar->get(phb_in, 1, NLVL, NLAT, NLON)){
			 std::cout << "Failed to get phb_in" << std::endl;
 		return false;}
			 if (!pVar->get(p_in, 1, NLVL, NLAT, NLON)){
			 std::cout << "Failed to get p_in" << std::endl;
	  return false;}
			 if (!pbVar->get(pb_in, 1, NLVL, NLAT, NLON)){
			 std::cout << "Failed to get pb_in" << std::endl;
	  return false;}
			 if (!tempVar->get(t_in, 1, NLVL, NLAT, NLON)){
			 std::cout << "Failed to get temperature" << std::endl;
	  return false;}
			 if (!qvVar->get(qv, 1, NLVL, NLAT, NLON)){
			 std::cout << "Failed to get qvapor" << std::endl;
	  return false;}
			 if (!qrVar->get(qr, 1, NLVL, NLAT, NLON)){
			 std::cout << "No qrain" << std::endl;
 		return false;}
	} // next record


		std::cout << "*** SUCCESS record Data" << std::endl;


		for (int i = 0; i < NLON; i++) {
 	   for (int j = 0; j < NLAT; j++) {
 	      for (int k = 0; k < NLVL - 1; k++) {
 			  u[k*NLON*NLAT + j*NLON + i] = (uwind_in[k*(NLON*1)*NLAT + j*(NLON+1)+ i] + uwind_in[k*(NLON+1)*NLAT + j*(NLON+1) + i+1])/2.0;
 			  v[k*NLON*NLAT + j*NLON + i] = (vwind_in[k*NLON*(NLAT+1) + j*NLON + i] + vwind_in[k*NLON*(NLAT+1) + (j+1)*NLON + i])/2.0;
 			  w[k*NLON*NLAT + j*NLON + i] = (wwind_in[k*NLON*NLAT + j*NLON + i] + wwind_in[(k+1)*NLON*NLAT + j*NLON + i])/2.0;
 			  z[k*NLON*NLAT + j*NLON + i] = (ph_in[k*NLON*NLAT + j*NLON + i] + ph_in[(k+1)*NLON*NLAT + j*NLON + i]
 				  + phb_in[k*NLON*NLAT + j*NLON + i] + phb_in[(k+1)*NLON*NLAT + j*NLON + i])/(2.0*9.81);
				p[k*NLON*NLAT + j*NLON + i] = p_in[k*NLON*NLAT + j*NLON + i] + pb_in[k*NLON*NLAT + j*NLON + i];
				t[k*NLON*NLAT + j*NLON + i] = t_in[k*NLON*NLAT + j*NLON + i] + 300.0;
				double PRESV = (p[k*NLON*NLAT + j*NLON + i] * qv[k*NLON*NLAT + j*NLON + i])/(0.622+qv[k*NLON*NLAT + j*NLON + i]);
				double PRESD = p[k*NLON*NLAT + j*NLON + i] - PRESV;
				rhoa[k*NLON*NLAT + j*NLON + i] = PRESD/(Rd*t[k*NLON*NLAT + j*NLON + i]) + PRESV/(Rv*t[k*NLON*NLAT + j*NLON + i]);
				qr[k*NLON*NLAT + j*NLON + i] = 0.001;
 		  }
 	  }
   }

 	  std::cout << "*** SUCCESS get Data" << std::endl;
		real spec = 2000;
		real x, y;
		real gateLat, gateLon;
		// std::vector<real> wrf_lat, wrf_lon;
		// for (int i = 0; i < NLON; i++){
		// 	for(int j = 0; j < 5; j++){
		// 		x = i*spec;
		// 		y = j*spec;
		// 		projection.Reverse(148, x, y, gateLat, gateLon);
		// 		wrf_lat.push_back(gateLat);
		// 		wrf_lon.push_back(gateLon);
		// 	}
		//
		// }

		// std::string timet = "978325200";
		std::string timet = "2001-01-01_05:00:00";
		datetime reftime = ParseTime("2001-01-01_05:00:00", "%Y-%m-%d_%H:%M:%S");
		std::cout << "Found matching reference time " << PrintDate(reftime) << ", " << PrintTime(reftime) << std::endl;

		std::string obFilename = dataPath + "/samurai_terrain.in";
		std::cout << "ObFile = " << obFilename << std::endl;
		std::cout << "Dimensions = (" << NLON << ", " << NLAT << ", " << NLVL << ") " << std::endl;
	 	ofstream obstream(obFilename);
		ostream_iterator<double> od(obstream, "\t ");
 		ostream_iterator<string> oi(obstream, "\t ");
		ostream_iterator<double> on(obstream, "\n");
		real refX, refY;
		projection.Forward(0, 0, 0, refX, refY);
		if (NLAT < 5){
		for (int i = 0; i < NLON; i++) {
			for (int j = 0; j < 6; j++) {
				x = i*spec;
				y = j*spec;
				projection.Reverse(0.0, refX + x, refY + y, gateLat, gateLon);
				// for (int k = 0; k < NLVL; k++) {
				//	*oi++ = timet;
					*od++ = gateLat;
					*od++ = gateLon;
					*od++ = x;
					*od++ = y;
					// *od++ = z[k*NLON*NLAT + 0*NLON + i]; // in m
					// *od++ = u[k*NLON*NLAT + 0*NLON + i];
					// *od++ = v[k*NLON*NLAT + 0*NLON + i];
					// *od++ = w[k*NLON*NLAT + 0*NLON + i];
					// *od++ = t[k*NLON*NLAT + 0*NLON + i];
					// *od++ = qv[k*NLON*NLAT + 0*NLON + i];
					// *od++ = rhoa[k*NLON*NLAT + 0*NLON + i];
					// *on++ = qr[k*NLON*NLAT + 0*NLON + i];
					*on++ = terrain_hgt[0*NLON + i]; // in m
				// }
			}
			}
		}
		else{
 		for (int i = 0; i < NLON; i++) {
 		 for (int j = 0; j < NLAT; j++) {
			  x = i*spec;
			  y = j*spec;
				// real x0, y0;
				//
				projection.Reverse(0.0, refX + x, refY + y, gateLat, gateLon);
 				for (int k = 0; k < NLVL; k++) {
 					*oi++ = timet;
					*od++ = gateLat;
					*od++ = gateLon;
					*od++ = z[k*NLON*NLAT + j*NLON + i];
 					*od++ = u[k*NLON*NLAT + j*NLON + i];
					*od++ = v[k*NLON*NLAT + j*NLON + i];
					*od++ = w[k*NLON*NLAT + j*NLON + i];
					*od++ = t[k*NLON*NLAT + j*NLON + i];
					*od++ = qv[k*NLON*NLAT + j*NLON + i];
					*od++ = rhoa[k*NLON*NLAT + j*NLON + i];
					*on++ = qr[k*NLON*NLAT + j*NLON + i];
					*on++ = terrain_hgt[j*NLON + i];
 			}
 		}
	}}
	 obstream << endl;
 	 std::cout << "*** SUCCESS print Data to samurai_Background.in" << std::endl;


	 return true;
  }
	// std::cout

  // TODO The same size variables are also read in VarDriver3D::validateFractlGrid()
  //      Combine into common code?


	// TODO Specify the timestep of the WRF file that you want to look at
	// Add it to the Args file
  // Nc3Dim *timeDim = dataFile.get_dim("TIME");
  // if (! timeDim)
  //   return false;

//   Nc3Dim *z0Dim = dataFile.get_dim("z0");
//   if (! z0Dim)
//     return false;
//
//   Nc3Dim *y0Dim = dataFile.get_dim("XLAT");
//   if (! y0Dim)
//     return false;
//
//   Nc3Dim *x0Dim = dataFile.get_dim("XLONG");
//   if (! x0Dim)
//     return false;
//
//   long ntime = timeDim->size();
//   long nz0 = z0Dim->size();
//   long ny0 = y0Dim->size();
//   long nx0 = x0Dim->size();
//
//   Nc3Var *time0 = dataFile.get_var("time");
//   if (! time0)
//     return false;
//
//   Nc3Var *z0 = dataFile.get_var("z0");
//   if (! z0)
//     return false;
//
//   Nc3Var *y0 = dataFile.get_var("y0");
//   if (! y0)
//     return false;
//   Nc3Var *x0 = dataFile.get_var("x0");
//   if (! x0)
//     return false;
//   Nc3Var *lat0 = dataFile.get_var("lat0");
//   if (! lat0)
//     return false;
//   Nc3Var *lon0 = dataFile.get_var("lon0");
//   if (! lon0)
//     return false;
//
//   Nc3Var *upward_air_velocity = dataFile.get_var("W");
//   if (! upward_air_velocity)
//     return false;
//   Nc3Var *northward_wind = dataFile.get_var("V");
//   if (! northward_wind)
//     return false;
//   Nc3Var *eastward_wind = dataFile.get_var("U");
//   if (! eastward_wind)
//     return false;
//   Nc3Var *meanNbrDbz = dataFile.get_var("DBZ");
//   if (! meanNbrDbz)
//     return false;
//
//   Nc3Var *meanNeighborNcp = dataFile.get_var("NCP");
//   if (! meanNeighborNcp)
//     return false;
//
//   Nc3Var *W_std = dataFile.get_var("W_std");
//   if (! W_std)
//     return false;
//   Nc3Var *V_std = dataFile.get_var("V_std");
//   if (! V_std)
//     return false;
//   Nc3Var *U_std = dataFile.get_var("U_std");
//   if (! U_std)
//     return false;
//
//   long numBgObs = ntime * nz0 * ny0 * nx0;
//
//   if (numBgObs <= 0) {
//     std::cout << "No background observations loaded" << std::endl;
//     return false;
//   }
//   std::cout << "Loading " << numBgObs << " background observations" << std::endl;
//
//   // Allocate data arrays
//
//   double *times = new double[ntime];
//
//   double *lats = new double[nx0 * ny0];
//   double *lons = new double[nx0 * ny0];
//   double *alts = new double[nz0];
//
//   double *upWind = new double[numBgObs];
//   double *northWind = new double[numBgObs];
//   double *eastWind = new double[numBgObs];
//   double *meanDbz =  new double[numBgObs];
//
//   double *wStd = new double[numBgObs];
//   double *vStd = new double[numBgObs];
//   double *uStd = new double[numBgObs];
//
//   bool success = true;
//
//   // Make sure table allocations worked.
//   // alts was allocated last so if any are going to fail that would be the one.
//
//   if (alts == NULL) {
//     std::cout << "Failed to allocate altitude table." << std::endl;
//     success = false;
//   }
//
//   // Get the time array
//
//   if ( success && ! time0->get(times, 1, ntime) ) {
//     std::cout << "Failed to read in the times." << std::endl;
//     success = false;
//   }
//
//   // Get the heights in meter
//
//   if ( success && ! z0->get(alts, nz0) ) {
//     std::cout << "Failed to read in altitudes." << std::endl;
//     success = false;
//   }
//
//   for (int rec = 0; rec < ntime; rec++) {
//     if ( ! upward_air_velocity->set_cur(rec, 0, 0, 0)) {
//       std::cout << "Failed to read in upward air velocity" << std::endl;
//       success = false;
//     }
//     if ( success && ! upward_air_velocity->get(upWind, 1, nz0, ny0, nx0)) {
//       std::cout << "Failed to read in upward air velocity" << std::endl;
//       success = false;
//     }
//
//     if ( success && ! northward_wind->set_cur(rec, 0, 0, 0)) {
//       std::cout << "Failed to read in northward air velocity" << std::endl;
//       success = false;
//     }
//     if ( success && ! northward_wind->get(northWind, 1, nz0, ny0, nx0)) {
//       std::cout << "Failed to read in northward air velocity" << std::endl;
//       success = false;
//     }
//
//     if ( success && ! eastward_wind->set_cur(rec, 0, 0, 0)) {
//       std::cout << "Failed to read in eastward air velocity" << std::endl;
//       success = false;
//     }
//     if ( success && ! eastward_wind->get(eastWind, 1, nz0, ny0, nx0)) {
//       std::cout << "Failed to read in eastward air velocity" << std::endl;
//       success = false;
//     }
//
//     if ( success && ! meanNbrDbz->set_cur(rec, 0, 0, 0)) {
//       std::cout << "Failed to read in mean dbZ" << std::endl;
//       success = false;
//     }
//     if ( success && ! meanNbrDbz->get(meanDbz, 1, nz0, ny0, nx0)) {
//       std::cout << "Failed to read in mean dbZ" << std::endl;
//       success = false;
//     }
//
//     if ( success && ! W_std->set_cur(rec, 0, 0, 0)) {
//       std::cout << "Failed to read in W_std" << std::endl;
//       success = false;
//     }
//
//     if ( success && ! W_std->get(wStd, 1, nz0, ny0, nx0)) {
//       std::cout << "Failed to read in W_std" << std::endl;
//       success = false;
//     }
//
//     if ( success && ! V_std->set_cur(rec, 0, 0, 0)) {
//       std::cout << "Failed to read in V_std" << std::endl;
//       success = false;
//     }
//
//     if ( success && ! V_std->get(vStd, 1, nz0, ny0, nx0)) {
//       std::cout << "Failed to read in V_std" << std::endl;
//       success = false;
//     }
//
//     if ( success && ! U_std->set_cur(rec, 0, 0, 0)) {
//       std::cout << "Failed to read in U_std" << std::endl;
//       success = false;
//     }
//
//     if ( success && ! W_std->get(uStd, 1, nz0, ny0, nx0)) {
//       std::cout << "Failed to read in U_std" << std::endl;
//       success = false;
//     }
//
//   }
//
//   double fillVal = 0.0;
//   // double fillVal = std::nan(""); // testing. Quick 2 iteration
//
//   for (int t = 0; t < ntime; t++) {
//     // double currentTime = times[t];
//
//     for (int z = 0; z < nz0; z++) {
//
//       float heightm = alts[z];	// TODO check unit
//       float rho = refstate->getReferenceVariable(ReferenceVariable::rhoref, heightm);
//
//       for (int y = 0; y < ny0; y++) {
// 	for (int x = 0; x < nx0; x++) {
// 	  int vIndex = x + nx0 * (y + ny0 * z);		// variables index
// 	  int bIndex = vIndex * numVars;		// index into bgU
//
// 	  double rhou = eastWind[vIndex] * rho;
// 	  if (std::isnan(rhou))
// 	    rhou = fillVal;
// 	  double rhov = northWind[vIndex] * rho;
// 	  if (std::isnan(rhov))
// 	    rhov = fillVal;
// 	  double rhow = upWind[vIndex] * rho;
// 	  if (std::isnan(rhow))
// 	    rhow = fillVal;
// 	  double dbz  = meanDbz[vIndex];
// 	  if (std::isnan(dbz))
// 	    dbz = -35.0;		// clear air.
//
// 	  bgU[bIndex + 0] = rhou;
// 	  bgU[bIndex + 1] = rhov;
// 	  bgU[bIndex + 2] = rhow;
// 	  bgU[bIndex + 3] = 0.0;			// t'
// 	  bgU[bIndex + 4] = 0.0;			// qv'
// 	  bgU[bIndex + 5] = 0.0;			// rho'
// 	  bgU[bIndex + 6] = (dbz + 35.) * 0.1;		// transform to scaled value [0..7]
//
// 	  // W_std, V_std, U_std are handled separatively in CostFunction3D
// 	}
//       }
//     }
//   } // ntime
//
//   if (isTrue("debug_bgu"))
//     dumpBgU(0, uStateSize, bgU);
//
//   if (configHash->exists("debug_bgu_nc"))
//     bgu2nc((*configHash)["debug_bgu_nc"].c_str(), bgU);
//
//    delete[] lats;
//    delete[] lons;
//    delete[] alts;
//    delete[] upWind;
//    delete[] northWind;
//    delete[] eastWind;
//    delete[] meanDbz;
//
//    delete[] wStd;
//    delete[] vStd;
//    delete[] uStd;
//
//    return success;
// }
