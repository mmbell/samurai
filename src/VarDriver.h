/*
 *  VarDriver.h
 *  samurai
 *
 *  Created by Michael Bell on 4/12/08.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVER_H
#define VARDRIVER_H

#include "BSpline.h"
#include "Observation.h"
#include "MetObs.h"
#include "FrameCenter.h"
#include "ReferenceState.h"
#include "Xml.h"
#include "HashMap.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>


#include <string>

#include "precision.h"
#include "Projection.h"
#include "samurai.h"

using namespace std;

class VarDriver
{

 public:

  VarDriver();
  virtual ~VarDriver();
  virtual bool initialize(const XMLNode& configuration) = 0;
  virtual bool initialize(const samurai_config &configSam) = 0;
  virtual bool initialize() = 0;

  virtual bool run() = 0;
  virtual bool run(int nx, int ny, int nsigma,

		   // ----- new -----
		   char cdtg[10],	// "12Z oct 4 2015 -> "2015100412"
		   int delta,	// delta * iter1 past cdtg
		   int iter1,
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
		   float *psam) = 0;

  virtual bool finalize() = 0;

  void clearCenters();
  void appendCenter(std::string date, std::string time, float lat, float lon, float Vm, float Um);
  void popCenter();

	// bool readTerrain(std::string &filename, std::vector<MetObs>* metData);

  HashMap *getConfigHash() { return &configHash; }

 protected:

  bool fixedGrid;	// Indicates if the grid dims come from the config file or run call

  real CoriolisF;
  real Pi;
  unsigned int numVars;
  unsigned int numHeights;
  unsigned int maxHeights;
  unsigned int numDerivatives;
  unsigned int maxJdim;
  unsigned int maxKdim;
  vector<FrameCenter> frameVector;
  HashMap configHash;
  ReferenceState* refstate;
  // Data Processing
  std::unordered_map<std::string, int> dataSuffix;
  std::string dataPath;
  enum dataFormats {
    unknown,
    cen,
    frd,
    cls,
    sec,
    ten,
    swp,
    sfmr,
    wwind,
    qscat,
    ascat,
    nopp,
    eol,
    cimss,
    dwl,
    insitu,
    mtp,
    mesonet,
    classnc,
    qcf,
    aeri,
    rad,
    cfrad,
    terrain,
    model,
    crsim,
    hrdradial,
    hdob
  };
  Projection projection;

  bool read_met_obs_file(int suffix, std::string &filename, std::vector<MetObs>* metObVector);
  bool read_frd(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_cls(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_wwind(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_eol(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_sec(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_ten(std::string& filename, std::vector<MetObs>*metObVector);
  bool read_dorade(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_sfmr(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_qscat(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_ascat(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_nopp(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_cimss(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_dwl(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_insitu(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_mtp(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_mesonet(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_classnc(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_qcf(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_aeri(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_rad(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_cfrad(std::string &fileName, std::vector<MetObs>* metObVector);
  bool read_terrain(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_model(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_crsim(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_hrdradial(std::string& filename, std::vector<MetObs>* metObVector);
  bool read_hdobs(std::string& filename, std::vector<MetObs>* metObVector); 
  bool readFrameCenters();
  bool parseXMLconfig(const XMLNode& config);
  bool parseSamuraiConfig(const samurai_config &config);
  Projection::ProjectionType projectionFromConfig();
};

#endif
