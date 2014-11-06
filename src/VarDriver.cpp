/*
 *  VarDriver.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriver.h"
#include "Dorade.h"
#include "ReferenceState.h"
#include <fstream>
#include <cmath>
#include <QTextStream>
#include <QDomDocument>
#include <QDomNodeList>
#include <GeographicLib/TransverseMercatorExact.hpp>
#include <netcdfcpp.h>

// Constructor
VarDriver::VarDriver()
{
	// Constant for all drivers
	Pi = acos(-1);

	// Set up the datatype hash
	dataSuffix["cen"] = cen;
	dataSuffix["frd"] = frd;
	dataSuffix["cls"] = cls;
	dataSuffix["sec"] = sec;
	dataSuffix["ten"] = ten;
	dataSuffix["swp"] = swp;
	dataSuffix["sfmr"] = sfmr;
	dataSuffix["Wwind"] = wwind;
	dataSuffix["eol"] = eol;
	dataSuffix["qscat"] = qscat;
	dataSuffix["ascat"] = ascat;
	dataSuffix["nopp"] = nopp;
	dataSuffix["cimss"] = cimss;
	dataSuffix["dwl"] = dwl;
	dataSuffix["insitu"] = insitu;
    dataSuffix["mtp"] = mtp;
    dataSuffix["mesonet"] = mesonet;
    dataSuffix["classnc"] = classnc;
    dataSuffix["qcf"] = qcf;
    dataSuffix["aeri"] = aeri;
}

// Destructor
VarDriver::~VarDriver()
{

}

/* This routine reads a text file containing a list of 1 second centers for the inertial reference frame
 File must be named yyyyMMdd.cen and has the format:
 HHmmss lat lon Vm Um
 where Vm and Um are the frame motion in m/s */

bool VarDriver::readFrameCenters()
{
	// Check the data directory for a centerfile
	QStringList filenames = dataPath.entryList();
	QString centerFilename;
	for (int i = 0; i < filenames.size(); ++i) {
		QString file = filenames.at(i);
		QStringList fileparts = file.split(".");
		if (fileparts.isEmpty()) {
			continue;
		}
		QString suffix = fileparts.last();
		if (suffix == "cen") {
			// Match to centerfile
			centerFilename = file;
			break;
		}
	}

	// Open the file
	QFile centerFile(dataPath.filePath(centerFilename));
	if (!centerFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
		std::cout << "Unable to open centerfile " << centerFilename.toAscii().data() << std::endl;
		return false;
	}

	// Get the date from the filename
	QString datestr = centerFilename.left(8);
	QDate startDate = QDate::fromString(datestr, "yyyyMMdd");

	// Read the centers
	QTextStream in(&centerFile);
	while (!in.atEnd()) {
		QDate date;
		QString line = in.readLine();
		QStringList lineparts = line.split(QRegExp("\\s+"));
		QString timestr = lineparts[0];
		int hour = timestr.left(2).toInt();
		if (hour > 23) {
			date = startDate.addDays(1);
			hour -= 24;
			QString newhr;
			newhr.setNum(hour);
			if (hour < 10) {
				timestr.replace(0,1,"0");
				timestr.replace(1,1,newhr);
			} else {
				timestr.replace(0,2,newhr);
			}
		} else {
			date = startDate;
		}
		QTime time = QTime::fromString(timestr, "HHmmss");
		QDateTime datetime = QDateTime(date, time, Qt::UTC);
		float lat = lineparts[1].toFloat();
		float lon = lineparts[2].toFloat();
		float Vm = lineparts[3].toFloat();
		float Um = lineparts[4].toFloat();
		FrameCenter center(datetime, lat, lon, Um, Vm);
		frameVector.push_back(center);
	}

	return true;
}

/* This routine reads the FRD insitu format from NOAA/HRD */

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
			real altT = lineparts[5].toFloat();
			real altW = lineparts[12].toFloat();
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

/* This routine reads the old CLASS dropsonde format from NCAR */

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
			int msec = lineparts[1].toInt()*1000;
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
				real w = lineparts[10].toFloat()+(-0.01*lineparts[2].toFloat()+22.);
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

/* This routine reads the modified CLASS dropsonde format from NCAR for the TPARC/TCS08 Dataset (and maybe TREX?) */

bool VarDriver::read_wwind(QFile& metFile, QList<MetObs>* metObVector)
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
		} else if (line.startsWith("UTC")) {
			datestr = line.mid(43,12);
			timestr = line.mid(57,8);
			QDate date = QDate::fromString(datestr, "yyyy, MM, dd");
			QTime time = QTime::fromString(timestr, "HH:mm:ss");
			datetime = QDateTime(date, time, Qt::UTC);
		} else if (line.startsWith("------")) {
			// Start reading data
			start = true;
		} else if (start) {
			MetObs ob;
			ob.setStationName(aircraft);
			line = QString(" ") + line;
			QStringList lineparts = line.split(QRegExp("\\s+"));
			int sec = (int)lineparts[1].toFloat();
			ob.setTime(datetime.addSecs(sec));
			if (lineparts[15].toFloat() != -999.) {
				ob.setLon(lineparts[15].toFloat());
			} else {
				ob.setLon(-999.);
			}
			if (lineparts[16].toFloat() != -999.) {
				ob.setLat(lineparts[16].toFloat());
			} else {
				ob.setLat(-999.);
			}
			if (lineparts[14].toFloat() != -999.0) {
				ob.setAltitude(lineparts[14].toFloat());
			} else {
				ob.setAltitude(-999.);
			}
			if (lineparts[5].toFloat() != -999.0) {
				ob.setPressure(lineparts[5].toFloat());
			} else {
				ob.setPressure(-999.0);
			}
			if (lineparts[6].toFloat() != -999.0) {
				ob.setTemperature(lineparts[6].toFloat() + 273.15);
			} else {
				ob.setTemperature(-999.);
			}
			if (lineparts[8].toFloat() != -999.0) {
				ob.setRH(lineparts[8].toFloat());
			} else {
				ob.setRH(-999.);
			}
			if (lineparts[12].toFloat() != -999.0) {
				ob.setWindDirection(lineparts[12].toFloat());
			} else {
				ob.setWindDirection(-999.);
			}
			if (lineparts[11].toFloat() != -999.0) {
				ob.setWindSpeed(lineparts[11].toFloat());
			} else {
				ob.setWindSpeed(-999.);
			}
			if (lineparts[18].toFloat() != 999.0) {
				ob.setVerticalVelocity(lineparts[18].toFloat());
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

/* This routine reads the modified CLASS dropsonde format from NCAR for the PREDICT field project */

bool VarDriver::read_eol(QFile& metFile, QList<MetObs>* metObVector)
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
		} else if (line.startsWith("UTC")) {
			datestr = line.mid(43,12);
			timestr = line.mid(57,8);
			QDate date = QDate::fromString(datestr, "yyyy, MM, dd");
			QTime time = QTime::fromString(timestr, "HH:mm:ss");
			datetime = QDateTime(date, time, Qt::UTC);
		} else if (line.startsWith("------")) {
			// Start reading data
			start = true;
		} else if (start) {
			MetObs ob;
			ob.setStationName(aircraft);
			line = QString(" ") + line;
			QStringList lineparts = line.split(QRegExp("\\s+"));
			int sec = (int)lineparts[1].toFloat();
			ob.setTime(datetime.addSecs(sec));
			if (lineparts[15].toFloat() != -999.) {
				ob.setLon(lineparts[15].toFloat());
			} else {
				ob.setLon(-999.);
			}
			if (lineparts[16].toFloat() != -999.) {
				ob.setLat(lineparts[16].toFloat());
			} else {
				ob.setLat(-999.);
			}
			if (lineparts[14].toFloat() != -999.0) {
				ob.setAltitude(lineparts[14].toFloat());
			} else {
				ob.setAltitude(-999.);
			}
			if (lineparts[5].toFloat() != -999.0) {
				ob.setPressure(lineparts[5].toFloat());
			} else {
				ob.setPressure(-999.0);
			}
			if (lineparts[6].toFloat() != -999.0) {
				ob.setTemperature(lineparts[6].toFloat() + 273.15);
			} else {
				ob.setTemperature(-999.);
			}
			if (lineparts[8].toFloat() != -999.0) {
				ob.setRH(lineparts[8].toFloat());
			} else {
				ob.setRH(-999.);
			}
			if (lineparts[12].toFloat() != -999.0) {
				ob.setWindDirection(lineparts[12].toFloat());
			} else {
				ob.setWindDirection(-999.);
			}
			if (lineparts[11].toFloat() != -999.0) {
				ob.setWindSpeed(lineparts[11].toFloat());
			} else {
				ob.setWindSpeed(-999.);
			}
			ob.setObType(MetObs::dropsonde);
			metObVector->push_back(ob);
		}
	}
	metFile.close();
	return true;
}

/* This routine reads a 1 sec flight level data file similar to the format provided by NOAA/HRD
 Columns are:
  TIME       Lat       Lon     Head   Track      GSpd     TAS    GAlt    Press     WndDr     WndSp   Tempr    Dewpt   DVal     PAlt     SurfP    VtWnd
 Pitch     Roll   Drift    Theta    Theta-e SFMRDown  SFMRSide */

bool VarDriver::read_sec(QFile& metFile, QList<MetObs>* metObVector)
{
	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	QString datestr, timestr, aircraft;
	QFileInfo info(metFile);
	QString fileName = info.fileName();
	aircraft = fileName.indexOf(QRegExp("HINUE"));
	datestr = fileName.left(8);
	QDate startDate = QDate::fromString(datestr, "yyyyMMdd");
	MetObs ob;
	while (!in.atEnd()) {
		QString line = in.readLine();
		line.remove(0, 1);
		QDateTime datetime;
		QDate date;
		ob.setStationName(aircraft);
		timestr = line.left(6);
		int hour = timestr.left(2).toInt();
		if (hour > 23) {
			date = startDate.addDays(1);
			hour -= 24;
			QString newhr;
			newhr.setNum(hour);
			if (hour < 10) {
				timestr.replace(0,1,"0");
				timestr.replace(1,1,newhr);
			} else {
				timestr.replace(0,2,newhr);
			}
		} else {
			date = startDate;
		}
		QTime time = QTime::fromString(timestr, "HHmmss");
		datetime = QDateTime(date, time, Qt::UTC);
		ob.setTime(datetime);

		QStringList lineparts = line.split(QRegExp("\\s+"));
		ob.setLat(lineparts[1].toFloat());

		// Note that the longitude can sometimes be given in degrees West, which jacks up the typically negative decimal representation for US locations;
		ob.setLon(lineparts[2].toFloat());
		ob.setAltitude(lineparts[7].toFloat());
		ob.setPressure(lineparts[8].toFloat());
		if ((lineparts[11].toFloat() != -999) or (lineparts[11].toFloat() != -32767)) {
			ob.setTemperature(lineparts[11].toFloat() + 273.15);
		} else {
			ob.setTemperature(-999.);
		}
		if ((lineparts[12].toFloat() != -999) or (lineparts[12].toFloat() != -32767)) {
			ob.setDewpoint(lineparts[12].toFloat() + 273.15);
		} else {
			ob.setDewpoint(-999.);
		}
		if (lineparts[10].toFloat() >= 0) {
			ob.setWindDirection(lineparts[9].toFloat());
			ob.setWindSpeed(lineparts[10].toFloat());
		}
		if ((lineparts[16].toFloat() != -999) or (lineparts[16].toFloat() != -32767)) {
			ob.setVerticalVelocity(lineparts[16].toFloat());
		}
		ob.setObType(MetObs::flightlevel);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;

}

/* This routine reads a 10 sec flight level data file from the data provided by the USAF Hurricane Hunters */
/* GMT Time   AOA     BSP   CAS    CC     CSP    DPR     DVAL      GA     GPSA   GS      HSS     LAT      LON       PA PITCH       RA ROLL     SLP   SS     TA   TAS    TDA    TDD  THD   TRK    TT   V V     WD     WS Valid Flags Source Tags */

bool VarDriver::read_ten(QFile& metFile, QList<MetObs>* metObVector)
{
	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	QString datestr, timestr, aircraft;
	QFileInfo info(metFile);
	QString fileName = info.fileName();
	aircraft = fileName.indexOf(QRegExp("HINUE"));
	datestr = fileName.left(8);
	QDate startDate = QDate::fromString(datestr, "yyyyMMdd");
	QDate date;
	MetObs ob;
	while (!in.atEnd()) {
		QString line = in.readLine();
		QDateTime datetime;
		ob.setStationName(aircraft);
		timestr = line.left(8);
		int hour = timestr.left(2).toInt();
		if (hour > 23) {
			date = startDate.addDays(1);
			hour -= 24;
			QString newhr;
			newhr.setNum(hour);
			if (hour < 10) {
				timestr.replace(0,1,"0");
				timestr.replace(1,1,newhr);
			} else {
				timestr.replace(0,2,newhr);
			}
		} else {
			date = startDate;
		}
		QTime time = QTime::fromString(timestr, "HH:mm:ss");
		datetime = QDateTime(date, time, Qt::UTC);
		ob.setTime(datetime);

		QStringList lineparts = line.split(QRegExp("\\s+"));
		real lat = lineparts[12].toFloat();
		ob.setLat(lat);
		real lon = lineparts[13].toFloat();
		ob.setLon(lon);
		ob.setAltitude(lineparts[16].toFloat());

		// Calculate pressure from palt using HRD formula
		real press = lineparts[5].toFloat();
		//real press =  1013.25*pow((1-palt/44331.),(1/.190263));
		ob.setPressure(press);
		ob.setTemperature(lineparts[20].toFloat() + 273.15);
		ob.setDewpoint(lineparts[23].toFloat() + 273.15);
		ob.setWindDirection(lineparts[28].toFloat());
		ob.setWindSpeed(lineparts[29].toFloat()*0.514444);
		ob.setVerticalVelocity(-999);
		ob.setObType(MetObs::flightlevel);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;
}

/* This routine reads a Dorade Sweepfile
	This does not correctly read 'pure' Dorade files due to big-endian representation in that format,
	so it needs a byte swap which has not been implemented yet*/

bool VarDriver::read_dorade(QFile& metFile, QList<MetObs>* metObVector)
{

	Dorade swpfile(metFile.fileName());

	// Use a Transverse Mercator projection to map the radar gates to the grid
	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();

	QString radardbz = configHash.value("radar_dbz");
	QString radarvel = configHash.value("radar_vel");
	QString radarsw = configHash.value("radar_sw");
	if(!swpfile.readSwpfile(radardbz, radarvel, radarsw))
		return false;

	int rayskip = configHash.value("radar_skip").toInt();
	int minstride = configHash.value("radar_stride").toInt();
	bool dynamicStride = configHash.value("dynamic_stride").toInt();
	int stride = minstride;
	for (int i=0; i < swpfile.getNumRays(); i+=rayskip) {
		real radarLat = swpfile.getRadarLat(i);
		real radarLon = swpfile.getRadarLon(i);
		real radarAlt = swpfile.getRadarAlt(i);
		real az = swpfile.getAzimuth(i);
		real el = swpfile.getElevation(i);
		float* refdata = swpfile.getReflectivity(i);
		float* veldata = swpfile.getRadialVelocity(i);
		float* swdata = swpfile.getSpectrumWidth(i);
		QDateTime rayTime = swpfile.getRayTime(i);
		float* gatesp = swpfile.getGateSpacing();
		real gatelength = gatesp[1] - gatesp[0];
		real beamwidth = sin(swpfile.getBeamwidthDeg()*Pi/180.);

		for (int n=0; n < swpfile.getNumGates()-stride; n+=stride) {
			MetObs ob;
			real range = gatesp[n+stride/2];
			if (dynamicStride) {
				stride = (int)(range*beamwidth/gatelength);
				if (stride < minstride) stride = minstride;
			}
			real dz = 0;
			real vr = 0;
			real sw = 0;
			real count = 0;
			for (int g=n; g<(n+stride); g++) {
				if (veldata[g] == -32768) continue;
				if (refdata[g] == -32768) continue;
				if (swdata[g] == -32768) continue;
				if (gatesp[g] <= 0) continue;
				dz += pow(10.0,(refdata[g]*0.1));
				vr += veldata[g];
				sw += swdata[g];
				count++;
			}
			if (count > 0) {
				dz = dz/count;
				dz = 10*log10(dz);
				vr = vr/count;
				sw = sw/count;
				real relX = range*sin(az*Pi/180.)*cos(el*Pi/180.);
				real relY = range*cos(az*Pi/180.)*cos(el*Pi/180.);
				real rEarth = 6371000.0;

				// Take into account curvature of the earth for the height of the radar beam
				real relZ = sqrt(range*range + rEarth*rEarth + 2.0 * range * rEarth * sin(el*Pi/180.)) - rEarth;
				real radarX, radarY, gateLat, gateLon;
				tm.Forward(radarLon, radarLat, radarLon, radarX, radarY);
				tm.Reverse(radarLon, radarX + relX, radarY + relY, gateLat, gateLon);
				real gateAlt = relZ + radarAlt*1000;

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
				/* cout << rayTime.toString(Qt::ISODate).toStdString() << "\t"
				<< gateLat << "\t" << gateLon << "\t" << gateAlt << "\t"
				<< az << "\t" << el << "\t" << dz << "\t" << vr << "\t" << sw << endl; */
			}
		}
	}

	return true;

}


/* This routine reads an old SFMR data format from CBLAST, not sure how current this is */

bool VarDriver::read_sfmr(QFile& metFile, QList<MetObs>* metObVector)
{

	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	QString datestr, timestr;
	QDateTime datetime;
	QTime time;
	QDate date;
	MetObs ob;
	while (!in.atEnd()) {
		QString line = in.readLine();
		QStringList lineparts = line.split(QRegExp("\\s+"));
		ob.setStationName("sfmr");
		timestr = lineparts[3];
		time = QTime::fromString(timestr, "HH:mm:ss");
		if (lineparts[2].toFloat() < 10) {
			datestr = lineparts[1] + "0" + lineparts[2] + lineparts[4];
		} else {
			datestr = lineparts[1] + lineparts[2] + lineparts[4];
		}
		date = QDate::fromString(datestr, "MMMddyyyy");
		datetime = QDateTime(date, time, Qt::UTC);
		ob.setTime(datetime);
		//QString datestr = datetime.toString(Qt::ISODate);
		//cout << datestr.toAscii().data() << endl;
		ob.setLat(lineparts[5].toFloat());
		ob.setLon(lineparts[6].toFloat());
		ob.setAltitude(10.0);
		ob.setWindSpeed(lineparts[7].toFloat());
		ob.setObType(MetObs::sfmr);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;

}

/* This routine reads an ASCII dump of the Level II Quikscat data, see code for columns*/

bool VarDriver::read_qscat(QFile& metFile, QList<MetObs>* metObVector)
{

	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	QString year,jd,timestr;
	MetObs ob;
	while (!in.atEnd()) {
		QString line = in.readLine();
		QStringList lineparts = line.split(QRegExp("\\s+"));
		QDateTime datetime;
		ob.setStationName("qscat");
		year = line.mid(56,4);
		jd = line.mid(61,3);
		QDate date(year.toInt(),1,1);
		date = date.addDays(jd.toInt()-1);
		timestr = line.mid(65,12);
		QTime time = QTime::fromString(timestr, "HH:mm:ss.zzz");
		datetime = QDateTime(date, time, Qt::UTC);
		ob.setTime(datetime);
		QString lat = line.mid(12,5);
		QString lon = line.mid(18,6);
		QString speed = line.mid(28,5);
		QString dir = line.mid(34,6);
		QString rainprob = line.mid(46,5);
		real rainflag = rainprob.toFloat();
		if (rainflag > 0.50) { continue; }
		real dirdeg = dir.toFloat() + 180.;
		if (dirdeg > 360.) dirdeg -= 360.;
		ob.setLat(lat.toFloat());
		ob.setLon(lon.toFloat());
		ob.setAltitude(10.0);
		ob.setWindSpeed(speed.toFloat());
		ob.setWindDirection(dirdeg);
		ob.setObType(MetObs::qscat);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;

}

/* This routine reads an ASCII dump of the ASCAT data, see code for columns*/

bool VarDriver::read_ascat(QFile& metFile, QList<MetObs>* metObVector)
{

	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	MetObs ob;
	QString unixtime = in.readLine();
	QDateTime datetime;
	datetime = QDateTime::fromTime_t(unixtime.toUInt());
	while (!in.atEnd()) {
		QString line = in.readLine();
		QStringList lineparts = line.split(QString(","));
		ob.setStationName("ascat");
		ob.setTime(datetime);
		ob.setLat(lineparts[0].toFloat());
		ob.setLon(lineparts[1].toFloat());
		ob.setAltitude(10.0);
		real speed = lineparts[2].toFloat() / 1.94384449;
		ob.setWindSpeed(speed);
		ob.setWindDirection(lineparts[3].toFloat());
		ob.setObType(MetObs::ascat);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;

}

/* This routine reads an ASCII dump of the Quikscat data for NOPP, see code for columns*/

bool VarDriver::read_nopp(QFile& metFile, QList<MetObs>* metObVector)
{

	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	MetObs ob;
	QDate startDate;
	QDateTime datetime;
	// Skip three lines
	QString line;
	in.readLine(); in.readLine(); in.readLine(); in.readLine();
	line = in.readLine();
	QStringList lineparts = line.split(QRegExp("\\s+"));
	//ob.setStationName(lineparts.last());
	ob.setStationName("qscat");
	startDate = QDate::fromString(lineparts[1], "yyyyMMdd");
	datetime.setDate(startDate);
	while (!in.atEnd()) {
		line = in.readLine();
		QStringList lineparts = line.split(QRegExp("\\s+"));
		lineparts.removeFirst();
		ob.setLat(lineparts[0].toFloat());
		ob.setLon(lineparts[1].toFloat());
		int secs = (int)lineparts[2].toFloat();
		ob.setTime(datetime.addSecs(secs));
		ob.setAltitude(10.0);
		ob.setWindSpeed(lineparts[3].toFloat());
		ob.setWindDirection(lineparts[4].toFloat());
		ob.setObType(MetObs::qscat);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;

}

/* This routine reads the CIMSS atmospheric motion vector files.
	Pressure levels are converted to Heights via Newton's method and the background
       reference state */

bool VarDriver::read_cimss(QFile& metFile, QList<MetObs>* metObVector)
{

	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	MetObs ob;
	QDateTime datetime;
	// Skip first line
	in.readLine();

	QString line;
	while (!in.atEnd()) {
		line = in.readLine();
		QStringList lineparts = line.split(QRegExp("\\s+"));
		QString station = lineparts[0] + " " + lineparts[1];
		ob.setStationName(station);
		QDate date = QDate::fromString(lineparts[2], "yyyyMMdd");
		QTime time = QTime::fromString(lineparts[3], "HHmm");
		datetime = QDateTime(date, time, Qt::UTC);
		ob.setTime(datetime);
		ob.setLat(lineparts[4].toFloat());
		ob.setLon(-lineparts[5].toFloat());
		// Convert pressure to altitude using Newton's method
		real presslevel = 100*lineparts[6].toFloat();
		real height = 9000;
		real pmin = 1e34;
		int iter = 0;
		while ((fabs(pmin) > 1.) and (iter < 5000)) {
			real p = refstate->getReferenceVariable(ReferenceVariable::pressref, height)-presslevel;
			real pprime = (refstate->getReferenceVariable(ReferenceVariable::pressref, height+500.)
                           - refstate->getReferenceVariable(ReferenceVariable::pressref, height-500.))/1000.;
			if (pprime != 0) {
				height = height - p/pprime;
				pmin = p;
			}
			iter++;
		}
		ob.setAltitude(height);
		ob.setWindSpeed(lineparts[7].toFloat());
		ob.setWindDirection(lineparts[8].toFloat());
		ob.setObType(MetObs::AMV);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;

}

/* This routine reads the LOS Doppler Wind Lidar data from TPARC, see code for columns*/

bool VarDriver::read_dwl(QFile& metFile, QList<MetObs>* metObVector)
{

	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();
	QTextStream in(&metFile);
	// Skip two lines
	in.readLine(); in.readLine();

	QDateTime datetime;
	QDate date;
	QString line;
	bool start = false;
	while (!in.atEnd()) {
		line = in.readLine();
		if (line.contains("rawfile")) {
			// Get the datetime
			date = QDate::fromString(line.mid(10,8), "yyyyMMdd");
		} else if (line.contains("%~~")) {
			start = true;
		} else if (start) {
			QStringList lineparts = line.split(QRegExp("\\s+"));
			QTime time = QTime::fromString(lineparts[0], "HHmmss");
			datetime = QDateTime(date, time, Qt::UTC);
			real radarLat = lineparts[1].toFloat();
			real radarLon = lineparts[2].toFloat();
			real radarAlt = lineparts[3].toFloat();
			real achdg = lineparts[4].toFloat();
			real az = lineparts[6].toFloat() + achdg;
			real el = lineparts[7].toFloat();

			// Lidar data is sufficiently sparse compared to radar data to include it all for now
			// Could implement some data thinning here
			int stride = 1;

			for (int n=0; n < lineparts[10].toInt(); n+=stride) {
				MetObs ob;
				int offset = n*5;
				real range = lineparts[12+offset].toFloat();
				real dz = 0;
				real vr = 0;
				real sw = 0;
				real count = 0;
				for (int g=n; g<(n+stride); g++) {
					offset = g*5;
					real veldata = lineparts[14+offset].toFloat();
					real refdata = lineparts[15+offset].toFloat();
					real swdata = lineparts[16+offset].toFloat();
					if (veldata == -999) continue;
					dz += pow(10.0,(refdata*0.1));
					vr += veldata;
					sw += swdata;
					count++;
				}
				if (count > 0) {
					dz = dz/count;
					dz = 10*log10(dz);
					vr = vr/count;
					sw = sw/count;
					real relX = range*sin(az*Pi/180.)*cos(el*Pi/180.);
					real relY = range*cos(az*Pi/180.)*cos(el*Pi/180.);
					real rEarth = 6371000.0;
					// Take into account curvature of the earth
					real relZ = sqrt(range*range + rEarth*rEarth + 2.0 * range * rEarth * sin(el*Pi/180.)) - rEarth;
					real radarX, radarY, gateLat, gateLon;
					tm.Forward(radarLon, radarLat, radarLon, radarX, radarY);
					tm.Reverse(radarLon, radarX + relX, radarY + relY, gateLat, gateLon);
					real gateAlt = relZ + radarAlt;
					ob.setObType(MetObs::lidar);
					ob.setLat(gateLat);
					ob.setLon(gateLon);
					ob.setAltitude(gateAlt);
					ob.setAzimuth(az);
					ob.setElevation(el);
					ob.setRadialVelocity(vr);
					ob.setReflectivity(dz);
					ob.setSpectrumWidth(sw);
					ob.setTime(datetime);
					metObVector->push_back(ob);
					/*cout << datetime.toString(Qt::ISODate).toStdString() << "\t"
					 << gateLat << "\t" << gateLon << "\t" << gateAlt << "\t"
					 << az << "\t" << el << "\t" << dz << "\t" << vr << "\t" << sw << endl; */
				}
			}
		}
	}

	return true;

}

/* This routine reads a generic insitu format suitable for many different observation types
TIME       Lat       Lon    Altitude   Pressure  Temperature (C)  Dewpoint (C)  WindDir   WindSpeed  VerticalVelocity*/
bool VarDriver::read_insitu(QFile& metFile, QList<MetObs>* metObVector)
{
	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	QString datestr, timestr, platform;
	QDateTime datetime;
	QFileInfo info(metFile);
	QString fileName = info.fileName();
	datestr = fileName.left(8);
	QDate startDate = QDate::fromString(datestr, "yyyyMMdd");
	// Get the platform name
	platform = fileName;
	platform.remove(0,8);
	platform.chop(7);
	MetObs ob;
	ob.setStationName(platform);
	while (!in.atEnd()) {
		QString line = in.readLine();
		QDate date;
		timestr = line.left(6);
		int hour = timestr.left(2).toInt();
		if (hour > 23) {
			date = startDate.addDays(1);
			hour -= 24;
			QString newhr;
			newhr.setNum(hour);
			if (hour < 10) {
				timestr.replace(0,1,"0");
				timestr.replace(1,1,newhr);
			} else {
				timestr.replace(0,2,newhr);
			}
		} else {
			date = startDate;
		}
		QTime time = QTime::fromString(timestr, "HHmmss");
		datetime = QDateTime(date, time, Qt::UTC);
		ob.setTime(datetime);

		QStringList lineparts = line.split(QRegExp("\\s+"));
		if ((lineparts[1].toFloat() != -32768) or (lineparts[1].toFloat() != -32767)) {
			ob.setLat(lineparts[1].toFloat());
		} else {
			ob.setLat(-999.0);
		}
		if ((lineparts[2].toFloat() != -32768) or (lineparts[2].toFloat() != -32767)) {
			ob.setLon(lineparts[2].toFloat());
		} else {
			ob.setLon(-999.0);
		}
		if ((lineparts[3].toFloat() != -32768) or (lineparts[3].toFloat() != -32767)) {
			ob.setAltitude(lineparts[3].toFloat());
		} else {
			ob.setAltitude(-999.0);
		}
		if (lineparts[4].toFloat() > 0) {
			ob.setPressure(lineparts[4].toFloat());
		} else {
			ob.setPressure(-999.0);
		}
		if (lineparts[5].toFloat() > -273) {
			float temp = lineparts[5].toFloat();
			if (temp < 100) temp += 273.15;
			ob.setTemperature(temp);
		} else {
			ob.setTemperature(-999.0);
		}
		if (lineparts[6].toFloat() > -273) {
			float dewp = lineparts[6].toFloat();
			if (dewp < 100) dewp += 273.15;
			ob.setDewpoint(dewp);
		} else {
			ob.setDewpoint(-999.0);
		}
		if (lineparts[8].toFloat() >= 0) {
			ob.setWindDirection(lineparts[7].toFloat());
			ob.setWindSpeed(lineparts[8].toFloat());
		} else {
			ob.setWindDirection(-999.0);
			ob.setWindSpeed(-999.0);
		}
		if ((lineparts[9].toFloat() != -32768) or (lineparts[9].toFloat() != -32767)) {
			ob.setVerticalVelocity(lineparts[9].toFloat());
		}  else {
			ob.setVerticalVelocity(-999.0);
		}
		ob.setObType(MetObs::insitu);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;
}

bool VarDriver::read_mtp(QFile& metFile, QList<MetObs>* metObVector)
{
	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	QString datestr, timestr, platform;
	QDateTime datetime;
	QFileInfo info(metFile);
	QString fileName = info.fileName();
	datestr = fileName.mid(2,8);
	QDate startDate = QDate::fromString(datestr, "yyyyMMdd");
	// Get the platform name
	platform = "GV";
	MetObs ob;
	ob.setStationName(platform);
	while (!in.atEnd()) {
		QString line = in.readLine();
        QStringList lineparts = line.split(QRegExp("\\s+"));
        lineparts.removeFirst();
        if (lineparts.size() > 5) {
            QDate date;
            date = startDate;
            int UTSEC = lineparts[0].toInt();
            QTime time(0,0,0);
            time = time.addSecs(UTSEC);
            QString timestring = time.toString();
            if (UTSEC > 86400) date = startDate.addDays(1);
            datetime = QDateTime(date, time, Qt::UTC);
            ob.setTime(datetime);
            ob.setLat(lineparts[10].toFloat());
            ob.setLon(lineparts[11].toFloat());
		} else {
            if (lineparts[1].toFloat() < 9999) {
                float temp = lineparts[1].toFloat();
                ob.setTemperature(temp);
                if (lineparts[4].toFloat() != 99999.0) {
                    float pressure = temp * lineparts[4].toFloat() * 0.000138;
                    ob.setPressure(pressure);
                } else {
                    ob.setPressure(-999.0);
                }
            } else {
                ob.setTemperature(-999.0);
                ob.setPressure(-999.0);
            }
            if (lineparts[2].toFloat() < 9999) {
                float error = lineparts[2].toFloat();
                ob.setTemperatureError(error);
            } else {
                ob.setTemperatureError(1.0);
            }
            if (lineparts[3].toFloat() != 99999) {
                float altitude = lineparts[3].toFloat();
                ob.setAltitude(altitude);
            } else {
                ob.setAltitude(-999.0);
            }
            ob.setObType(MetObs::mtp);
            metObVector->push_back(ob);
        }
	}

	metFile.close();
	return true;
}

/* This routine parses the supplied XML configuration
 and validates parameters that are common to all drivers
 Mode specific configs are validated in those drivers */

bool VarDriver::parseXMLconfig(const QDomElement& config)
{

	cout << "Parsing configuration file...\n";

	// Parse the nodes to a hash
	QDomNodeList nodeList = config.childNodes();
	for (int i = 0; i < nodeList.count(); i++) {
		QDomNode currNode = nodeList.item(i);
		// Check to see if this is a set of pass parameters
		QString iter = 0;
		if (currNode.hasAttributes() and currNode.attributes().contains("iter")) {
			iter = currNode.toElement().attribute("iter");
		}
		QDomNodeList configList = currNode.childNodes();
		for (int j = 0; j < configList.count(); j++) {
			QDomNode configItem = configList.item(j);
			QDomElement group = configItem.toElement();
			QString tag = group.tagName();
			if (iter.toInt() > 1) {
				// Append the pass number to the tagName
				tag += "_" + iter;
			}
			if (!group.text().isEmpty()) {
				configHash.insert(tag, group.text());
				cout << tag.toStdString() << " => " << configHash.value(tag).toStdString() << endl;
			}
		}
	}

	// Validate the hash -- multiple passes are not validated currently
	QStringList configKeys;
	configKeys << "ref_state" << "ref_time" << // "reflat" << "reflon" are set by the VarDriver
	"qr_variable" << "i_background_roi" << "j_background_roi" <<
	"i_reflectivity_roi" << "j_reflectivity_roi" << "k_reflectivity_roi" <<
	"load_background" << "adjust_background" <<
	"radar_dbz" << "radar_vel" << "radar_sw" << "radar_skip" <<
    "radar_stride" << "dynamic_stride" << "dbz_pseudow_weight" <<
    "melting_zone_width" << "mixed_phase_dbz" << "rain_dbz" <<
	"num_iterations" << "output_mish" << "output_txt" << "output_qc" <<
    "output_netcdf" << "output_asi" << "preprocess_obs" << "mask_reflectivity" <<
    "dropsonde_rhou_error" << "dropsonde_rhov_error" << "dropsonde_rhow_error" <<
    "dropsonde_tempk_error" << "dropsonde_qv_error" << "dropsonde_rhoa_error" <<
    "flightlevel_rhou_error" << "flightlevel_rhov_error" << "flightlevel_rhow_error" <<
    "flightlevel_tempk_error" << "flightlevel_qv_error" << "flightlevel_rhoa_error" <<
    "insitu_rhou_error" << "insitu_rhov_error" << "insitu_rhow_error" <<
    "insitu_tempk_error" << "insitu_qv_error" << "insitu_rhoa_error" <<
    "sfmr_windspeed_error" << "qscat_rhou_error" << "qscat_rhov_error" <<
    "ascat_rhou_error" << "ascat_rhov_error" << "amv_rhou_error" << "amv_rhov_error" <<
    "lidar_sw_error" << "lidar_power_error" << "lidar_min_error" <<
    "radar_sw_error" << "radar_fallspeed_error" << "radar_min_error" <<
    "bg_rhou_error" << "bg_rhov_error" << "bg_rhow_error" << "bg_tempk_error" <<
    "bg_qv_error" << "bg_rhoa_error" << "bg_qr_error";
	for (int i = 0; i < configKeys.count(); i++) {
		if (!configHash.contains(configKeys.at(i))) {
			cout <<	"No configuration found for <" << configKeys.at(i).toStdString() << "> aborting..." << endl;
			return false;
		}
	}
	return true;

}

bool VarDriver::read_mesonet(QFile& metFile, QList<MetObs>* metObVector)
{
	NcError err(NcError::verbose_nonfatal);

	// Read the file.
	NcFile dataFile(metFile.fileName().toAscii(), NcFile::ReadOnly);

	// Check to see if the file was read.
	if(!dataFile.is_valid())
		return false;

    // Get the number of records
    NcDim* recnum;
    if (!(recnum = dataFile.get_dim("recNum")))
        return false;
    int NREC = recnum->size();

    NcVar *latVar, *lonVar, *altVar, *timeVar, *tempVar,
        *dewpVar, *wdirVar, *wspdVar, *pressVar;
    if (!(latVar = dataFile.get_var("latitude")))
        return false;
    if (!(lonVar = dataFile.get_var("longitude")))
        return false;
    if (!(altVar = dataFile.get_var("elevation")))
        return false;
    if (!(timeVar = dataFile.get_var("observationTime")))
        return false;
    if (!(tempVar  = dataFile.get_var("temperature")))
        return false;
    if (!(dewpVar = dataFile.get_var("dewpoint")))
        return false;
    if (!(wdirVar = dataFile.get_var("windDir")))
        return false;
    if (!(wspdVar = dataFile.get_var("windSpeed")))
        return false;
    if (!(pressVar = dataFile.get_var("stationPressure")))
        return false;

    real lat[NREC], lon[NREC], alt[NREC], obtime[NREC], temp[NREC],
        dewp[NREC], wdir[NREC], wspd[NREC], press[NREC];
    if (!latVar->get(lat, NREC))
        return false;
    if (!lonVar->get(lon, NREC))
        return false;
    if (!altVar->get(alt, NREC))
        return false;
    if (!timeVar->get(obtime, NREC))
        return false;
    if (!tempVar->get(temp, NREC))
        return false;
    if (!dewpVar->get(dewp, NREC))
        return false;
    if (!wdirVar->get(wdir, NREC))
        return false;
    if (!wspdVar->get(wspd, NREC))
        return false;
    if (!pressVar->get(press, NREC))
        return false;

    MetObs ob;
	ob.setStationName("Mesonet");
    for (int rec = 0; rec < NREC; rec++)
    {
        if ((lat[rec] != -9999.0) and (lat[rec] < 1.0e32)) {
			ob.setLat(lat[rec]);
		} else {
			ob.setLat(-999.0);
		}
        if ((lon[rec] != -9999.0) and (lon[rec] < 1.0e32)) {
			ob.setLon(lon[rec]);
		} else {
			ob.setLon(-999.0);
		}
        if ((alt[rec] != -9999.0) and (alt[rec] < 1.0e32)) {
			ob.setAltitude(alt[rec]);
		} else {
			ob.setAltitude(-999.0);
		}

        QDateTime datetime;
        datetime = QDateTime::fromTime_t(obtime[rec]);
        ob.setTime(datetime.toUTC());

        if ((temp[rec] != -9999.0) and (temp[rec] < 1.0e32)) {
			ob.setTemperature(temp[rec]);
		} else {
			ob.setTemperature(-999.0);
		}
        if ((dewp[rec] != -9999.0) and (dewp[rec] < 1.0e32)) {
			ob.setDewpoint(dewp[rec]);
		} else {
			ob.setDewpoint(-999.0);
		}
        if ((wdir[rec] != -9999.0) and (wdir[rec] < 1.0e32)) {
			ob.setWindDirection(wdir[rec]);
		} else {
			ob.setWindDirection(-999.0);
		}
        if ((wspd[rec] != -9999.0) and (wspd[rec] < 1.0e32)) {
			ob.setWindSpeed(wspd[rec]);
		} else {
			ob.setWindSpeed(-999.0);
		}
        if ((press[rec] != -9999.0) and (press[rec] < 1.0e32)) {
			ob.setPressure(press[rec]/100.0);
		} else {
			ob.setPressure(-999.0);
		}

        ob.setObType(MetObs::mesonet);
		metObVector->push_back(ob);
    }

    return true;

}

bool VarDriver::read_classnc(QFile& metFile, QList<MetObs>* metObVector)
{
	NcError err(NcError::verbose_nonfatal);

	// Read the file.
	NcFile dataFile(metFile.fileName().toAscii(), NcFile::ReadOnly);

	// Check to see if the file was read.
	if(!dataFile.is_valid())
		return false;

    // Get the number of records
    NcDim* recnum;
    if (!(recnum = dataFile.get_dim("time")))
        return false;
    int NREC = recnum->size();

    NcVar *latVar, *lonVar, *altVar, *timeVar, *tempVar,
    *rhVar, *wdirVar, *wspdVar, *pressVar;
    if (!(latVar = dataFile.get_var("lat")))
        return false;
    if (!(lonVar = dataFile.get_var("lon")))
        return false;
    if (!(altVar = dataFile.get_var("alt")))
        return false;
    if (!(timeVar = dataFile.get_var("time")))
        return false;
    if (!(tempVar  = dataFile.get_var("temp")))
        return false;
    if (!(rhVar = dataFile.get_var("rh")))
        return false;
    if (!(wdirVar = dataFile.get_var("wdir")))
        return false;
    if (!(wspdVar = dataFile.get_var("wspd")))
        return false;
    if (!(pressVar = dataFile.get_var("pres")))
        return false;

    real lat[NREC], lon[NREC], alt[NREC], obtime[NREC], temp[NREC],
        rh[NREC], wdir[NREC], wspd[NREC], press[NREC];
    if (!latVar->get(lat, NREC))
        return false;
    if (!lonVar->get(lon, NREC))
        return false;
    if (!altVar->get(alt, NREC))
        return false;
    if (!timeVar->get(obtime, NREC))
        return false;
    if (!tempVar->get(temp, NREC))
        return false;
    if (!rhVar->get(rh, NREC))
        return false;
    if (!wdirVar->get(wdir, NREC))
        return false;
    if (!wspdVar->get(wspd, NREC))
        return false;
    if (!pressVar->get(press, NREC))
        return false;

    QDateTime datetime;
    QStringList fileparts = metFile.fileName().split(".");
	QDate startDate = QDate::fromString(fileparts[1], "yyyyMMdd");
    QTime startTime = QTime::fromString(fileparts[2], "HHmmss");

	// Get the platform name
	MetObs ob;
	ob.setStationName(fileparts[0]);
    for (int rec = 0; rec < NREC; rec++)
    {
        if ((lat[rec] != -9999.0) and (lat[rec] < 1.0e32)) {
			ob.setLat(lat[rec]);
		} else {
			ob.setLat(-999.0);
		}
        if ((lon[rec] != -9999.0) and (lon[rec] < 1.0e32)) {
			ob.setLon(lon[rec]);
		} else {
			ob.setLon(-999.0);
		}
        if ((alt[rec] != -9999.0) and (alt[rec] < 1.0e32)) {
			ob.setAltitude(alt[rec]*1000.0);
		} else {
			ob.setAltitude(-999.0);
		}

        QTime time(startTime);
        time.addSecs(obtime[rec]);
        datetime = QDateTime(startDate, time, Qt::UTC);
        ob.setTime(datetime);

        if ((temp[rec] != -9999.0) and (temp[rec] < 1.0e32)) {
			ob.setTemperature(temp[rec]+273.15);
		} else {
			ob.setTemperature(-999.0);
		}
        if ((rh[rec] != -9999.0) and (rh[rec] < 1.0e32)) {
			ob.setRH(rh[rec]);
		} else {
			ob.setRH(-999.0);
		}
        if ((wdir[rec] != -9999.0) and (wdir[rec] < 1.0e32)) {
			ob.setWindDirection(wdir[rec]);
		} else {
			ob.setWindDirection(-999.0);
		}
        if ((wspd[rec] != -9999.0) and (wspd[rec] < 1.0e32)) {
			ob.setWindSpeed(wspd[rec]);
		} else {
			ob.setWindSpeed(-999.0);
		}
        if ((press[rec] != -9999.0) and (press[rec] < 1.0e32)) {
			ob.setPressure(press[rec]);
		} else {
			ob.setPressure(-999.0);
		}

        ob.setObType(MetObs::dropsonde);
		metObVector->push_back(ob);
    }

    return true;

}

/* This routine reads the qcf composite data format*/

bool VarDriver::read_qcf(QFile& metFile, QList<MetObs>* metObVector)
{

	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&metFile);
	MetObs ob;
	QDate startDate;
	QDateTime datetime;
	// Skip three lines
	QString line;
	in.readLine(); in.readLine(); in.readLine();
	while (!in.atEnd()) {
		line = in.readLine();
		QStringList lineparts = line.split(QRegExp("\\s+"));
        QDate date = QDate::fromString(lineparts[0], "yyyy/MM/dd");
        QTime time = QTime::fromString(lineparts[1], "HH:mm:ss");
		datetime = QDateTime(date, time, Qt::UTC);
        ob.setTime(datetime);
        ob.setStationName(lineparts[5]);
		ob.setLat(lineparts[6].toFloat());
		ob.setLon(lineparts[7].toFloat());
		ob.setAltitude(lineparts[9].toFloat());
        if (lineparts[11] == "G") {
            ob.setPressure(lineparts[10].toFloat());
        }
        if (lineparts[17] == "G") {
            ob.setTemperature(lineparts[16].toFloat() + 273.15);
        }
        if (lineparts[19] == "G") {
            ob.setDewpoint(lineparts[18].toFloat() + 273.15);
        }
        if (lineparts[21] == "G") {
            ob.setWindSpeed(lineparts[20].toFloat());
        }
        if (lineparts[23] == "G") {
            ob.setWindDirection(lineparts[22].toFloat());
        }
		ob.setObType(MetObs::mesonet);
		metObVector->push_back(ob);
	}

	metFile.close();
	return true;

}

bool VarDriver::read_aeri(QFile& metFile, QList<MetObs>* metObVector)
{
	NcError err(NcError::verbose_nonfatal);

	// Read the file.
	NcFile dataFile(metFile.fileName().toAscii(), NcFile::ReadOnly);

	// Check to see if the file was read.
	if(!dataFile.is_valid())
		return false;

    // Get the number of records
    NcDim* recnum;
    if (!(recnum = dataFile.get_dim("time")))
        return false;
    int NREC = recnum->size();

    // Get the number of pressure levels
    NcDim* pres_level;
    if (!(pres_level = dataFile.get_dim("pres_level")))
        return false;
    int NLVL = pres_level->size();

    NcVar *latVar, *lonVar, *altVarBase, *altVar, *timeVarBase, *timeVar,
    *tempVar, *dewpVar, *pressVar;
    if (!(latVar = dataFile.get_var("lat")))
        return false;
    if (!(lonVar = dataFile.get_var("lon")))
        return false;
    if (!(altVarBase = dataFile.get_var("alt")))
        return false;
    if (!(altVar = dataFile.get_var("height")))
        return false;
    if (!(timeVarBase = dataFile.get_var("base_time")))
        return false;
    if (!(timeVar = dataFile.get_var("time_offset")))
        return false;
    if (!(tempVar  = dataFile.get_var("ambientTemp")))
        return false;
    if (!(dewpVar = dataFile.get_var("dewpointTemp")))
        return false;
    if (!(pressVar = dataFile.get_var("pressure")))
        return false;

    real lat, lon, altBase, basetime, alt[NREC][NLVL], obtime[NREC], temp[NREC][NLVL],
    dewp[NREC][NLVL], press[NREC][NLVL];
    if (!latVar->get(&lat, 1))
        return false;
    if (!lonVar->get(&lon, 1))
        return false;
    if (!altVarBase->get(&altBase, 1))
        return false;
    if (!timeVarBase->get(&basetime, 1))
        return false;
    if (!altVar->get(&alt[0][0], NREC, NLVL))
        return false;
    if (!timeVar->get(&obtime[0], NREC, NLVL))
        return false;
    if (!tempVar->get(&temp[0][0], NREC, NLVL))
        return false;
    if (!dewpVar->get(&dewp[0][0], NREC, NLVL))
        return false;
    if (!pressVar->get(&press[0][0], NREC, NLVL))
        return false;

    QDateTime datetime;

	// Get the platform name
	MetObs ob;
    QStringList fileparts = metFile.fileName().split("_");
	ob.setStationName(fileparts[0]);
    for (int rec = 0; rec < NREC; rec++) {
        for (int p = 0; p < NLVL; p++) {
        if ((lat != -9999.0) and (lat < 1.0e32)) {
			ob.setLat(lat);
		} else {
			ob.setLat(-999.0);
		}
        if ((lon != -9999.0) and (lon < 1.0e32)) {
			ob.setLon(lon);
		} else {
			ob.setLon(-999.0);
		}
        if ((alt[rec][p] != -9999.0) and (alt[rec][p] < 1.0e32)) {
			ob.setAltitude(alt[rec][p] + altBase);
		} else {
			ob.setAltitude(-999.0);
		}

        datetime = QDateTime::fromTime_t(basetime);
        QString obstring = datetime.toString(Qt::ISODate);
        QDateTime obdatetime = datetime.addSecs(obtime[rec]).toUTC();
        obstring = obdatetime.toString(Qt::ISODate);
        ob.setTime(obdatetime);

        if ((temp[rec][p] != -9999.0) and (temp[rec][p] < 1.0e32)) {
			ob.setTemperature(temp[rec][p]);
		} else {
            //cout << temp[rec][p] << "\n";
			ob.setTemperature(-999.0);
		}
        if ((dewp[rec][p] != -9999.0) and (dewp[rec][p] < 1.0e32)) {
			ob.setDewpoint(dewp[rec][p]);
		} else {
			ob.setDewpoint(-999.0);
		}
        if ((press[rec][p] != -9999.0) and (press[rec][p] < 1.0e32)) {
			ob.setPressure(press[rec][p]);
		} else {
			ob.setPressure(-999.0);
		}

        ob.setObType(MetObs::aeri);
		metObVector->push_back(ob);
        }
    }

    return true;

}
