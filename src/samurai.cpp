#include <iostream>
#include <QtCore>
#include <QtXml>

#include "samurai.h"
#include "VarDriver3D.h"
#include "Xml.h"

using namespace std;

extern "C" {

  // Constructors
  
  // Call this one if you want to pass a config structure.
  // Very error prone. The Fortran and the C++ data structure must match exactly...
  
  VarDriver3D *create_vardriver3D(const samurai_config *config) {
    cout << "C API, create_vardriver3D" << endl;
    VarDriver3D *driver = new VarDriver3D();
    if (driver == NULL) {
      std:: cout << "Failed to create a VarDriver3D" << std::endl;
      return NULL;
    }
    
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

  VarDriver3D *create_vardriver3D_From_File(const char *config_path) {
    cout << "C API, create_vardriver3D_From_File" << endl;
    QDomDocument *domDoc = readXmlConfig(config_path);
    
    if (domDoc == NULL)
      return NULL;

    VarDriver3D *driver = new VarDriver3D();
    if (driver == NULL) {
      std:: cout << "Failed to create a VarDriver3D" << std::endl;
      delete domDoc;
      return NULL;
    }

    if ( ! driver->initialize(domDoc->documentElement() )) {
      std:: cout << "Failed to initialize VarDriver3D with content of '" << config_path << "'"  << std::endl;
      delete driver;
      delete domDoc;
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
		      float dx, float dy,
		      float *sigmas,	// 2D array (nsigma)
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
    if (d == NULL)
      return 0;
    return d->run(nx, ny, nsigma, dx, dy,
			   sigmas, latitude, longitude,
			   u1, v1, w1, th1, p1,
			   usam, vsam, wsam, thsam, psam);
  }

  // Call this to debug
  void dump_hash(QHash<QString, QString> &hash) {
    QHash<QString, QString>::iterator it;
    for(it = hash.begin(); it != hash.end(); ++it)
      std::cout << it.key().toLatin1().data() << ": "
		<< it.value().toLatin1().data() << std::endl;
  }
}
