#ifndef BKGD_OBS_LOADER_H
#define BKGD_OBS_LOADER_H

#include <iostream>
#include <vector>

#include <string>
#include <unordered_map>

#include <kd/kd.hh>	// KD Tree implementation from Lrose

#include "FrameCenter.h"
#include "BkgdAdapter.h"
#include "Projection.h"
#include "ReferenceState.h"

#include "precision.h"

// Various loader for background observations
// This allows having different interpolation methods

class BkgdObsLoader {

 public:

  BkgdObsLoader();
  virtual ~BkgdObsLoader();

  virtual bool loadBkgdObs(std::vector<real> &bgIn) = 0;

  bool timeCheck(real time, datetime& startTime, datetime& endTime, int &tci);

  virtual bool initialize(std::unordered_map<std::string, std::string> *config,
			  std::vector<FrameCenter> frames,
			  BkgdAdapter *adapter,
			  Projection proj,
			  ReferenceState *refSt,
			  int uStSize,
			  real *bgu,			// TODO s this right? vector<real> maybe
			  unsigned int varsNum,
			  int id, int jd, int kd,
			  float imn, float jmn, float kmn,
			  float imx, float jmx, float kmx,
			  real iinc, real jinc, real kinc);
  
  KD_tree *buildKDTree(std::vector<double> &bgIn);
  
  bool fillHoles(std::vector<int> &emptybg);
  bool isTrue(const char *flag_in) { 
		std::string flag = flag_in;
		if (configHash->find(flag) != configHash->end()) {
			if ((*configHash)[flag] == "true") {
				return true;
			} 
		}
		return false;
	}
  void dumpBgIn(int from, int to, std::vector<real> &bgIn);
  void dumpBgU(int from, int to, real *bgu);
  void bgu2nc(const char *fname, real *bgu);
  
  typedef BSpline<real> SplineD;
  typedef BSplineBase<real> SplineBase;
  
  typedef enum {
    BG_LOADER_PRE,
    BG_LOADER_SPLINE,
    BG_LOADER_KD,
    BG_LOADER_FRACTL
  } bg_loader_t;

  typedef enum {
    RUN_MODE_XYZ,
    RUN_MODE_RTZ
  } run_mode_t;
  
 protected:

  BkgdAdapter *bkgdAdapter;

  real iincr, jincr, kincr;
  int idim, jdim, kdim;
  float imin, jmin, kmin;
  float imax, jmax, kmax;
  
  std::unordered_map<std::string, std::string> *configHash;
  std::vector<FrameCenter> frameVector;
  Projection projection;

  run_mode_t runMode;

  ReferenceState *refstate;
  unsigned int numVars;
  int uStateSize;

  real *bgWeights;
  real *bgU;

  std::vector<real> logheights, uBG, vBG, wBG, tBG, qBG, rBG, zBG;

  std::string interp_mode;

  real bgX, bgY, bgRadius, bgTheta;
  real bgZ;

  real maxGridDist;
  real iROI, jROI;
  
};

// ----------------- Load preprocessed background obs --------------------

class BkgdObsPreLoader : public BkgdObsLoader {

 public:

  BkgdObsPreLoader();
  ~BkgdObsPreLoader();
  
  bool loadBkgdObs(std::vector<real> &bgIn);  
};


// --------------- Load using a spline -----------------------

// This is the traditional Samurai method.
// Uses
//   exponential interpolation in horizontal dimensions
//   b-spline interpolation on log height in vertical

class BkgdObsSplineLoader : public BkgdObsLoader {

 public:

  BkgdObsSplineLoader();
  ~BkgdObsSplineLoader();

  bool loadBkgdObs(std::vector<real> &bgIn);
  
  bool initialize(std::unordered_map<std::string, std::string> *config,
		  std::vector<FrameCenter> frames,
		  BkgdAdapter *adapter,
		  Projection proj,
		  ReferenceState *refSt,
		  int uStSize,
		  real *bgu,
		  unsigned int varsNum,
		  int id, int jd, int kd,
		  float imn, float jmn, float kmn,
		  float imx, float jmx, float kmx,
		  real iinc, real jinc, real kinc);

 protected:

  bool splineSolver(int waveLen);

};

// --------------- Load using KD Tree -----------------------
// This is the method used by FRACL (previously RadarWind)
// Uses a KD tree to find nearest neighbors

class BkgdObsKDLoader : public BkgdObsLoader {

 public:

  BkgdObsKDLoader();
  ~BkgdObsKDLoader();
  
  // bool initialize();
  bool loadBkgdObs(std::vector<real> &bgIn);
  
 protected:

  bool fillBguEntry(KD_real *centerLoc, int nbrMax, float maxDistance,
		    std::vector<real> &bgIn, KD_tree *kdTree, real *bgU, int bIndex,
		    std::vector<int> &emptybg, int debug);
  bool overwriteBgu(const char *fname);
};

// --------------- Load FRACTL generated Netcdf file -----------------------

class BkgdObsFractlLoader : public BkgdObsLoader {

 public:

  BkgdObsFractlLoader();
  ~BkgdObsFractlLoader();
  
  //  bool initialize();
  bool loadBkgdObs(std::vector<real> &bgIn);
  bool fillBguEntry(std::vector<real> &bgIn, real *bgU, int iIndex, int bIndex);
};

// ------------- Object Factory --------------
// Return a subclass of BkgdObsLoader depending on the value of t

class BkgdObsLoaderFactory {

 public:

  static  BkgdObsLoader *createBkgdObsLoader(BkgdObsLoader::bg_loader_t t);
};

#endif
