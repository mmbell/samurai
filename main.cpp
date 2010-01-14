
#include "precision.h"
#include "VarDriver1d.h"
#include "VarDriver2d.h"
#include "VarDriverRZ.h"
#include "VarDriverXY.h"
#include <iostream>
#include <QApplication>

int main (int argc, char *argv[]) {
	
	//QApplication app(argc, argv);
	if (argc >=2) {
		QString arg(argv[1]);
		if (arg == "XY") {
			VarDriverXY driver;
			driver.initialize();
			driver.run();
		} else {
			VarDriverRZ driver;
			driver.initialize();
			driver.run();
		}
	} else {
		VarDriverRZ driver;
		driver.initialize();
		driver.run();
	}
	std::cout << "Analysis complete!\n";
		
	return 0;
	
}

