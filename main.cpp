
#include "VarDriver1d.h"
#include "VarDriver2d.h"
#include <iostream>
#include <QApplication>

int main (int argc, char *argv[]) {
	
	//QApplication app(argc, argv);
	//VarDriver1d driver;
	VarDriver2d driver;
	if (!driver.run()) {
		std::cout << "Failed to run driver\n";
	}
	
	return 0;
	
}

