
#include "precision.h"
#include "VarDriver3D.h"
#include "Xml.h"
#include <iostream>
#include "timing/gptl.h"
#include "Args.h"
#include "Strings.h" // Helper functions, eg: 'endsWith' (which will be in C++20)


void usage() {
  std::cout << "Usage: samurai <samurai_configuration.xml>" << std::endl;
}

int main (int argc, char *argv[]) {
  // Initialize our timers:
  GPTLinitialize();
  GPTLstart("Total");
	// Set up our XML document:
  XMLDocument xml;

    GPTLstart("Main::Init");

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
  GPTLstop("Main::Init");
  
  // Do the analysis

#if !IO_BENCHMARK  
  GPTLstart("Main::Run");
  if(!driver->run()) {
    delete driver;
    return EXIT_FAILURE;
  }
  GPTLstop("Main::Run");
#endif
  GPTLstart("Main::Finalize");
  if(!driver->finalize()) {
    delete driver;
    return EXIT_FAILURE;
  }
  GPTLstop("Main::Finalize");

  //delete driver; // NCAR - uncommenting this results in a crash on GPU.. must still be something allocated in GPU space?
  std::cout << "Analysis successful!\n";
  GPTLstop("Total");
  GPTLpr(0);
  GPTLfinalize();
  return EXIT_SUCCESS;
}
