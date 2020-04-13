#include <iostream>

#include "samurai.h"
#include "VarDriver3D.h"
#include "Xml.h"

//
// This is the Samurai API definition.
// These functions are in the the samurai library. So you can access them
// from any language that let you call external functions.
// I have tested the interface with Fortran, Python, and Julia
// See https://github.com/mjskier/samurai_scripts
//

using namespace std;

extern "C" {

  // Constructors
  
  // Call this one if you want to pass a config structure.
  // Very error prone. The Fortran and the C++ data structure must match exactly...
  
  VarDriver3D *create_vardriver3D(const samurai_config *config, bool fixedGrid ) {
    cout << "C API, create_vardriver3D" << endl;
    VarDriver3D *driver = new VarDriver3D();
    if (driver == NULL) {
      std:: cout << "Failed to create a VarDriver3D" << std::endl;
      return NULL;
    }
    
    driver->setGridFlag(fixedGrid);
    if(!driver->initialize(*config)) {
      std:: cout << "Failed to initialize the VarDriver3D" << std::endl;
      delete driver;
      return NULL;
    }
    return driver;
  }

  // Pass the path to a config file.
  // I recommend this one. A lots simpler, and no need to match
  // data structure from 2 different languages

	XMLDocument xml;
  VarDriver3D *create_vardriver3D_From_File(const char *config_path, bool fixedGrid) {
    cout << "C API, create_vardriver3D_From_File" << endl;
    XMLError ec = xml.LoadFile(config_path);
    if (ec != XML_SUCCESS) {
      std::cout << "Error opening XML file: " << config_path << std::endl;
      return NULL;
    }

    XMLElement* root = xml.FirstChildElement("samurai");


    VarDriver3D *driver = new VarDriver3D();
    if (driver == NULL) {
      std:: cout << "Failed to create a VarDriver3D" << std::endl;
      return NULL;
    }
    driver->setGridFlag(fixedGrid);
    
    if ( ! driver->initialize(*root)) {
      std:: cout << "Failed to initialize VarDriver3D with content of '" << config_path << "'"  << std::endl;
      delete driver;
      return NULL;
    }      
    return driver;
  }

  // Destructor
  void delete_vardriver3D(VarDriver3D *d) {
    cout << "C API, delete_vardriver3D" << endl;
    if (d != NULL) 
      delete d;
  }
  
  // Run the analysis
  int run_vardriver3D(VarDriver3D *d,
		      // These are input values
		      int nx, int ny, int nsigma,

		      // ----- new -----
		      char cdtg[10],	// "12Z oct 4 2015 -> "2015100412"
		      int delta,	// delta * iter1 past cdtg
		      int iter1,
		      float imin, float imax, float iincr, // used to come from config
		      float jmin, float jmax, float jincr,
		      // ----- new -----
		      
		      float *sigmas,	// 1D array (nsigma)
		      float *latitude,	// 2D array (nx, ny)
		      float *longitude,	// 2D array
		      float *u1,	// 3D array (nx, ny, nsigma)  
		      float *v1,	// 3D array
		      float *w1,	// 3D array
		      float *th1,	// 3D array
		      float *p1,	// 3D array

		      // These are output values
		      float *usam,	// 3D array
		      float *vsam,	// 3D array
		      float *wsam,	// 3D array
		      float *thsam,	// 3D array
		      float *psam	// 3D array
		      ) {
    if (d == NULL) {
      cout << "run_vardriver3D was called with a NULL driver handle" << endl;
      return 0;
    }
    
    if (d->isFixedGrid()) {
      cout << "This version of run can only be called with non-fixed grids" << endl;
      return 0;
    }
    
    return d->run(nx, ny, nsigma,
		  cdtg, delta, iter1,
		  imin, imax, iincr,
		  jmin, jmax, jincr,
		  sigmas, latitude, longitude,
		  u1, v1, w1, th1, p1,
		  usam, vsam, wsam, thsam, psam);
  }

  // Call this to clear the center vector

  void clear_centers(VarDriver3D *d) {
    if (d == NULL)
      return;
    d->clearCenters();
  }
  
  void pop_center(VarDriver3D *d) {
    if (d == NULL)
      return;
    d->popCenter();
  }

  void append_center(VarDriver3D *d, char *date, char *time,
		     float lat, float lon,
		     float vm, float um) {
    if (d == NULL)
      return;
    d->appendCenter(date, time, lat, lon, vm, um);
  }
  
  // Call this to debug
  void dump_hash(std::unordered_map<std::string, std::string> &hash) {
		for (auto &entry : hash) {
			std::cout << entry.first << " : " << entry.second << std::endl;
		}
  }
}
