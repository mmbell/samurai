
#include "precision.h"
#include "VarDriver3D.h"
#include "Xml.h"
#include <iostream>
#include "Args.h"
#include "Strings.h" // Helper functions, eg: 'endsWith' (which will be in C++20)


void usage() {
  std::cout << "Usage: samurai <samurai_configuration.xml>" << std::endl;
}

int main (int argc, char *argv[]) {
	// Set up our XML document:
  XMLDocument xml;

  // Read the command line argument to get the XML configuration file

  // Generic driver which will be instanced by the configuration specification
  VarDriver *driver = new VarDriver3D();
  if (driver == NULL)
    return EXIT_FAILURE;

  switch (argc) {
    
  case 1:  // no argument given
    
    usage();
    return EXIT_FAILURE;
    break;
    
  default: // assume old way of specifying an .xml file

    std::string fname = argv[1];
    if (endsWith(fname, ".xml") ) {
      XMLError ec = xml.LoadFile(fname.c_str());
			if (ec != XML_SUCCESS) {
				std::cout << "Error opening XML file: " << fname << std::endl;
				return EXIT_FAILURE;
			}
    
      XMLElement* root = xml.FirstChildElement("samurai");
    
      if(!driver->initialize(*root)) {
				delete driver;
				return EXIT_FAILURE;
      }
    } else {	// New way, goes through the TDRP system
      Args args;
      if (args.parseArgs(argc, argv)) {
	args.paramsToHash(driver->getConfigHash());
	if (! driver->initialize()) {
	  delete driver;
	  return EXIT_FAILURE;
	}
      } else {
	return EXIT_FAILURE;
			}
    }
    // Parse was fine. Fall through
  }
  
  // Do the analysis
	
  if(!driver->run()) {
    delete driver;
    return EXIT_FAILURE;
  }
  if(!driver->finalize()) {
    delete driver;
    return EXIT_FAILURE;
  }
    
  delete driver;
  std::cout << "Analysis successful!\n";
  return EXIT_SUCCESS;
}
