
#include "precision.h"
#include "VarDriver1d.h"
#include "VarDriver2d.h"
#include "VarDriverRZ.h"
#include "VarDriverXY.h"
#include "VarDriverVAR.h"
#include "VarDriver3D.h"
#include <iostream>
#include <QApplication>

int main (int argc, char *argv[]) {
	
	//QApplication app(argc, argv);
	if (argc >=2) {
		QString arg(argv[1]);
		if (arg == "PX") {
			VarDriverXY driver;
			driver.initialize();
			driver.run();
		} else if (arg == "XY") {
			VarDriverVAR driver;
			driver.initialize();
			driver.run();
		} else if (arg == "RZ") {
			VarDriverRZ driver;
			driver.initialize();
			driver.run();
		} else if (arg == "XYZ") {
			VarDriver3D driver;
			driver.initialize();
			driver.run();
		}
		std::cout << "Analysis complete!\n";
	} else {
		std::cout << "Usage: samurai <mode>\n Available modes: XYZ, XY, PX, RZ\n";
	}
		
	return 0;
	
}

