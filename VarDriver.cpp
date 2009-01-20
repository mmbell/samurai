/*
 *  VarDriver.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriver.h"
#include "Dorade.h"
#include <fstream>
#include <cmath>
#include <QTextStream>

VarDriver::VarDriver()
{
	CoriolisF = 6e-5;
	Pi = acos(-1);
	
	// Set up the datatype hash
	dataSuffix["frd"] = frd;
	dataSuffix["cls"] = cls;
	dataSuffix["sec"] = sec;
	dataSuffix["ten"] = ten;
	dataSuffix["swp"] = swp;
	
}

VarDriver::~VarDriver()
{

}

bool VarDriver::readTCcenters()
{
	QFile centerFile("./centerfile.txt");
	if (!centerFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;
	
	QTextStream in(&centerFile);
	while (!in.atEnd()) {
		// Hard code for now
		QDate date(2003, 9, 12);
		
		QString line = in.readLine();
		QStringList lineparts = line.split(QRegExp("\\s+"));
		QString timestr = lineparts[0];
		if (timestr.left(2).toInt() > 23) {
			date = date.addDays(1);
		}
		QTime time = QTime::fromString(timestr, "HHmmss");
		QDateTime datetime = QDateTime(date, time, Qt::UTC);
		float lat = lineparts[1].toFloat();
		float lon = lineparts[2].toFloat();
		float Vm = lineparts[3].toFloat();
		float Um = lineparts[4].toFloat();
		TCcenter center(datetime, lat, lon, Um, Vm);
		tcVector.push_back(center);
	}
	
	return true;
}

bool VarDriver::read_frd(QFile& metFile, QList<MetObs>* metObVector)
{
	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	QString datestr, timestr, aircraft;
	QDateTime datetime;
	bool start = false;
	while (!in.atEnd()) {
		QString line = in.readLine();
		if (line.startsWith(" Aircraft:")) {
			aircraft = line.mid(11,11);
		} else if (line.startsWith(" Date:")) {
			datestr = "20" + line.mid(7,6);
		} else if (line.startsWith(" Time:")) {
			timestr = line.mid(7,6);
		} else if (line.startsWith("IX")) {
			// Start reading data
			start = true;
			QDate date = QDate::fromString(datestr, "yyyyMMdd");
			QTime time = QTime::fromString(timestr, "HHmmss");
			datetime = QDateTime(date, time, Qt::UTC);
		} else if (start) {
			MetObs ob;
			ob.setStationName(aircraft);
			QStringList lineparts = line.split(QRegExp("\\s+"));
			int msec = (int)lineparts[1].toFloat()*1000;
			ob.setTime(datetime.addMSecs(msec));
			ob.setLat(lineparts[17].toFloat());
			ob.setLon(lineparts[18].toFloat());
			float altT = lineparts[5].toFloat();
			float altW = lineparts[12].toFloat();
			ob.setAltitude((altT + altW)/2);
			ob.setPressure(lineparts[2].toFloat());
			if (lineparts[3] != "-999.00") {
				ob.setTemperature(lineparts[3].toFloat() + 273.15);
			} else {
				ob.setTemperature(-999.);
			}
			ob.setRH(lineparts[4].toFloat());
			ob.setWindDirection(lineparts[6].toFloat());
			ob.setWindSpeed(lineparts[7].toFloat());
			ob.setVerticalVelocity(lineparts[11].toFloat());
			ob.setObType(MetObs::dropsonde);
			metObVector->push_back(ob);
		}
	}
	metFile.close();
	return true;
}

bool VarDriver::read_cls(QFile& metFile, QList<MetObs>* metObVector)
{
	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;
	
	QTextStream in(&metFile);
	QString datestr, timestr, aircraft;
	QDateTime datetime;
	bool start = false;
	while (!in.atEnd()) {
		QString line = in.readLine();
		if (line.startsWith("Launch Site Type")) {
			QStringList lineparts = line.split(":");
			aircraft = lineparts[1].trimmed();
		} else if (line.startsWith("GMT")) {
			datestr = line.mid(35,12);
			timestr = line.mid(49,8);
			QDate date = QDate::fromString(datestr, "yyyy, MM, dd");
			QTime time = QTime::fromString(timestr, "HH:mm:ss");
			datetime = QDateTime(date, time, Qt::UTC);
		} else if (line.startsWith("------")) {
			// Start reading data
			start = true;
		} else if (start) {
			MetObs ob;
			ob.setStationName(aircraft);
			QStringList lineparts = line.split(QRegExp("\\s+"));
			int msec = (int)lineparts[1].toFloat()*1000;
			ob.setTime(datetime.addMSecs(msec));
			if (lineparts[11].toFloat() != 999.) { 
				ob.setLon(lineparts[11].toFloat());
			} else {
				ob.setLon(-999.);
			}
			if (lineparts[12].toFloat() != 999.) { 
				ob.setLat(lineparts[12].toFloat());
			} else {
				ob.setLat(-999.);
			}
			if (lineparts[15].toFloat() != 99999.0) { 
				ob.setAltitude(lineparts[15].toFloat());
			} else {
				ob.setAltitude(99999.0);
			}
			if (lineparts[2].toFloat() != 9999.0) { 
				ob.setPressure(lineparts[2].toFloat());
			} else {
				ob.setPressure(9999.0);
			}
			if (lineparts[3].toFloat() != 999.0) {
				ob.setTemperature(lineparts[3].toFloat() + 273.15);
			} else {
				ob.setTemperature(-999.);
			}
			if (lineparts[5].toFloat() != 999.0) {
				ob.setRH(lineparts[5].toFloat());
			} else {
				ob.setRH(-999.);
			}
			if (lineparts[9].toFloat() != 999.0) {
				ob.setWindDirection(lineparts[9].toFloat());
			} else {
				ob.setWindDirection(-999.);
			}
			if (lineparts[8].toFloat() != 999.0) {
				ob.setWindSpeed(lineparts[8].toFloat());
			} else {
				ob.setWindSpeed(-999.);
			}
			if ((lineparts[10].toFloat() != 99.0) and (lineparts[2].toFloat() != 9999.0)) {
				float w = lineparts[10].toFloat()+(-0.01*lineparts[2].toFloat()+22.);
				ob.setVerticalVelocity(w);
			} else {
				ob.setVerticalVelocity(-999.);
			}
			ob.setObType(MetObs::dropsonde);
			metObVector->push_back(ob);
		}
	}
	metFile.close();
	return true;
}

bool VarDriver::read_sec(QFile& metFile, QList<MetObs>* metObVector)
{
	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;
	
	QTextStream in(&metFile);
	QString datestr, timestr, aircraft;
	aircraft = metFile.fileName().indexOf(QRegExp("HINU"));
	//datestr = metFile.fileName().left(8);	
	//QDate date = QDate::fromString(datestr, "yyyyMMdd");
	// Hard code date for now
	QDate date(2003, 9, 12);
	MetObs ob;
	while (!in.atEnd()) {
		QString line = in.readLine();
		line.remove(0, 1);
		QDateTime datetime;
		ob.setStationName(aircraft);
		timestr = line.left(6);
		QTime time = QTime::fromString(timestr, "HHmmss");
		datetime = QDateTime(date, time, Qt::UTC);
		ob.setTime(datetime);

		QStringList lineparts = line.split(QRegExp("\\s+"));
		ob.setLat(lineparts[1].toFloat());
		ob.setLon(lineparts[2].toFloat());
		ob.setAltitude(lineparts[7].toFloat());
		ob.setPressure(lineparts[8].toFloat());
		ob.setTemperature(lineparts[11].toFloat() + 273.15);
		ob.setDewpoint(lineparts[12].toFloat() + 273.15);
		ob.setWindDirection(lineparts[9].toFloat());
		ob.setWindSpeed(lineparts[10].toFloat());
		ob.setVerticalVelocity(lineparts[16].toFloat());
		ob.setObType(MetObs::flightlevel);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;	
	
}

bool VarDriver::read_ten(QFile& metFile, QList<MetObs>* metObVector)
{
	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;
	
	QTextStream in(&metFile);
	QString datestr, timestr, aircraft;
	aircraft = metFile.fileName().indexOf(QRegExp("HINU"));
	//datestr = metFile.fileName().left(8);	
	//QDate date = QDate::fromString(datestr, "yyyyMMdd");
	// Hard code date for now
	QDate date(2003, 9, 12);
	MetObs ob;
	while (!in.atEnd()) {
		QString line = in.readLine();
		QDateTime datetime;
		ob.setStationName(aircraft);
		timestr = line.left(8);
		QTime time = QTime::fromString(timestr, "HH:mm:ss");
		datetime = QDateTime(date, time, Qt::UTC);
		ob.setTime(datetime);
		
		QStringList lineparts = line.split(QRegExp("\\s+"));
		lineparts[2].chop(1);
		double latmin = lineparts[2].toFloat();
		double lat = lineparts[1].toFloat() + latmin/60;
		ob.setLat(lat);
		lineparts[4].chop(1);
		double lonmin = lineparts[4].toFloat();
		double lon = -(lineparts[3].toFloat() + lonmin/60);
		ob.setLon(lon);
		ob.setAltitude(lineparts[11].toFloat());
		// Calculate pressure from palt using HRD formula
		double palt = lineparts[10].toFloat();
		double press =  1013.25*pow((1-palt/44331.),(1/.190263));
		ob.setPressure(press);
		ob.setTemperature(lineparts[5].toFloat() + 273.15);
		ob.setDewpoint(lineparts[6].toFloat() + 273.15);
		ob.setWindDirection(lineparts[7].toFloat());
		ob.setWindSpeed(lineparts[8].toFloat()*0.514444);
		ob.setVerticalVelocity(-999);
		ob.setObType(MetObs::flightlevel);
		metObVector->push_back(ob);
	}
	
	metFile.close();
	return true;	
}

bool VarDriver::read_dorade(QFile& metFile, QList<MetObs>* metObVector)
{

	Dorade swpfile(metFile.fileName());
	if(!swpfile.readSwpfile())
		return false;

	float radarLat = swpfile.getRadarLat();
	float radarLon = swpfile.getRadarLon();
	float radarAlt = swpfile.getRadarAlt();
	
	for (int i=0; i < swpfile.getNumRays(); i++) {
		float az = swpfile.getAzimuth(i);
		float el = swpfile.getElevation(i);
		float* refdata = swpfile.getReflectivity(i);
		float* veldata = swpfile.getRadialVelocity(i);	
		float* swdata = swpfile.getSpectrumWidth(i);
		QDateTime rayTime = swpfile.getRayTime(i);

		float* gatesp = swpfile.getGateSpacing();
		for (int n=0; n < swpfile.getNumGates(); n+=4) {
		//for (int n=0; n < swpfile.getNumGates(); n++) {
			MetObs ob;
			float dz = refdata[n];
			float vr = veldata[n];
			if (vr == -32768) continue;
			float sw = swdata[n];
			float range = gatesp[n];
			if (range <= 0) continue;
			
			double relX = -range*sin(az*Pi/180.)*cos(el*Pi/180.);
			double relY = -range*cos(az*Pi/180.)*cos(el*Pi/180.);
			double relZ = range*sin(el*Pi/180.);
			double latrad = radarLat * Pi/180.0;
			double fac_lat = 111.13209 - 0.56605 * cos(2.0 * latrad)
			+ 0.00012 * cos(4.0 * latrad) - 0.000002 * cos(6.0 * latrad);
			double fac_lon = 111.41513 * cos(latrad)
			- 0.09455 * cos(3.0 * latrad) + 0.00012 * cos(5.0 * latrad);
			double gateLon = radarLon - (relX/1000)/fac_lon;
			double gateLat = radarLat - (relY/1000)/fac_lat;
			double gateAlt = relZ + radarAlt*1000;
			ob.setObType(MetObs::radar);
			ob.setLat(gateLat);
			ob.setLon(gateLon);
			ob.setAltitude(gateAlt);
			ob.setAzimuth(az);
			ob.setElevation(el);
			ob.setRadialVelocity(vr);
			ob.setReflectivity(dz);
			ob.setSpectrumWidth(sw);
			ob.setTime(rayTime);
			metObVector->push_back(ob);
		}
	}
	
	return true;
	
}

