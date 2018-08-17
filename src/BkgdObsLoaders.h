#ifndef BKGD_OBS_LOADER_H
#define BKGD_OBS_LOADER_H

#include <iostream>
#include <vector>

#include <QString>
#include <QList>
#include <QVector>

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

  virtual bool loadBkgdObs(QList<real> &bgIn) = 0;

  bool timeCheck(real time, QDateTime &startTime, QDateTime &endTime, int &tci);

  virtual bool initialize(QHash<QString, QString> *config,
			  std::vector<FrameCenter> frames,
			  BkgdAdapter *adapter,
			  Projection proj,
			  ReferenceState *refSt,
			  int uStSize,
			  real *bgu,			// TODO s this right? QList<real> maybe
			  unsigned int varsNum,
			  int id, int jd, int kd,
			  float imn, float jmn, float kmn,
			  float imx, float jmx, float kmx,
			  real iinc, real jinc, real kinc);

  bool fillHoles(QVector<int> &emptybg);
  bool isTrue(const char *flag) { return configHash->contains(flag) && configHash->value(flag) == "true"; }  
  void dumpBgIn(int from, int to, QList<real> &bgIn);
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
  
  QHash<QString, QString> *configHash;
  std::vector<FrameCenter> frameVector;
  Projection projection;

  run_mode_t runMode;

  ReferenceState *refstate;
  unsigned int numVars;
  int uStateSize;

  real *bgWeights;
  real *bgU;

  QVector<real> logheights, uBG, vBG, wBG, tBG, qBG, rBG, zBG;

  QString interp_mode;

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
  
  bool loadBkgdObs(QList<real> &bgIn);  
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

  bool loadBkgdObs(QList<real> &bgIn);
  
  bool initialize(QHash<QString, QString> *config,
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
  bool loadBkgdObs(QList<real> &bgIn);
  
 protected:

  KD_tree *buildKDTree(QList<double> &bgIn);

  bool fillIt(KD_real *centerLoc, int nbrMax, float maxDistance,
	      QList<real> &bgIn, KD_tree *kdTree, real *bgU, int bIndex,
	      QVector<int> &emptybg, int debug);
  bool overwriteBgu(const char *fname);
};

// --------------- Load FRACTL generated Netcdf file -----------------------

class BkgdObsFractlLoader : public BkgdObsLoader {

 public:

  BkgdObsFractlLoader();
  ~BkgdObsFractlLoader();
  
  //  bool initialize();
  bool loadBkgdObs(QList<real> &bgIn);
  
};

// ------------- Object Factory --------------
// Return a subclass of BkgdObsLoader depending on the value of t

class BkgdObsLoaderFactory {

 public:

  static  BkgdObsLoader *createBkgdObsLoader(BkgdObsLoader::bg_loader_t t);
};

#endif
