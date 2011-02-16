
#include "precision.h"
#include "VarDriverRZ.h"
#include "VarDriverXY.h"
#include "VarDriverVAR.h"
#include "VarDriverXYZ.h"
#include <iostream>
#include <QtCore>
#include <QtXml>
//#include <QApplication>

int main (int argc, char *argv[]) {
	
	//QApplication app(argc, argv); It would be nice to wrap a GUI around this at some point
	
	// Declare a generic driver which will be instanced by the configuration specification
	VarDriver *driver = NULL;
	
	// Read the command line argument to get the XML configuration file
	if (argc >=2) {
		
		// Check to make sure the argument has the right suffix
		QString xmlfile(argv[1]);
		if (xmlfile.right(3) != "xml") {
			std::cout << xmlfile.toStdString() << " does not look like an XML file\n";
			return EXIT_FAILURE;
		}
		
		// Open the file
		QFile file(xmlfile);
		if (!file.open(QIODevice::ReadOnly)) {
			std::cout << "Error Opening Configuration File, Check Permissions on " << xmlfile.toStdString() << "\n";
			return EXIT_FAILURE;
		}
		
		// Create a DOM document with contents from the configuration file
		QDomDocument domDoc;
		QString errorStr;
		int errorLine;
		int errorColumn;
		if (!domDoc.setContent(&file, true, &errorStr, &errorLine, &errorColumn)) {
			// Exit on malformed XML
			QString errorReport = QString("XML Parse Error in "+xmlfile+" at Line %1, Column %2:\n%3")
			.arg(errorLine)
			.arg(errorColumn)
			.arg(errorStr);
			std::cout << errorReport.toStdString() << "\n";
			file.close();
			return EXIT_FAILURE;
		}
		
		// Successful file read
		file.close();
		
		// Check the root node to make sure this is really a SAMURAI configuration file
		QDomElement root = domDoc.documentElement();
		if (root.tagName() != "samurai") {
			std::cout << "The XML file " << xmlfile.toStdString() << " is not an SAMURAI configuration file\n.";
			return EXIT_FAILURE;
		}
		
		// Get the run 'mode' child node to instance the proper driver
		QDomNodeList nodeList = root.childNodes();
		for (int i = 0; i < nodeList.count(); i++) {
			QDomNode currNode = nodeList.item(i);
			QDomElement element = currNode.toElement();
			if (element.tagName() == "mode") {
				QString mode = element.text();
				if (mode == "XYZ") {
					driver = new VarDriverXYZ();
				} else {
					std:: cout << "Unsupported run mode " << mode.toStdString() << std::endl;
					return EXIT_FAILURE;
				}
				break;
			}
		}
		
		// Make sure we were able to create a driver, then drive
		if (driver != NULL) {
			// Do the analysis
			if(!driver->initialize(root)) {
				delete driver;
				return EXIT_FAILURE;
			}
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
		} else {
			return EXIT_FAILURE;
		}
	
	} else {
		std::cout << "Usage: samurai <samurai_configuration.xml>\n";
		return EXIT_SUCCESS;
	}
	
}

