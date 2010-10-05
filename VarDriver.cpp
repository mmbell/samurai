/*
 *  VarDriver.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "VarDriver.h"
#include "Dorade.h"
#include <fstream>
#include <cmath>
#include <QTextStream>
#include <QDomDocument>
#include <QDomNodeList>

VarDriver::VarDriver()
{
	CoriolisF = 6e-5;
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
	
}

VarDriver::~VarDriver()
{

}

bool VarDriver::readTCcenters()
{
	// Check the data directory for a centerfile
	QDir dataPath("./vardata");
	dataPath.setFilter(QDir::Files);
	dataPath.setSorting(QDir::Name);
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
	QFile centerFile(dataPath.filePath(centerFilename));
	if (!centerFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
		std::cout << "Unable to open centerfile " << centerFilename.toAscii().data() << std::endl;
		return false;
	}
	
	QString datestr = centerFilename.left(8);	
	QDate startDate = QDate::fromString(datestr, "yyyyMMdd");
	
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


bool VarDriver::read_sec(QFile& metFile, QList<MetObs>* metObVector)
{
	/* TIME       Lat       Lon     Head   Track      GSpd     TAS    GAlt    Press     WndDr     WndSp   Tempr    Dewpt   DVal     PAlt     
	 SurfP    VtWnd     Pitch     Roll   Drift    Theta    Theta-e SFMRDown  SFMRSide */
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
		//ob.setLon(-lineparts[2].toFloat());
		ob.setLon(lineparts[2].toFloat());
		ob.setAltitude(lineparts[7].toFloat());
		ob.setPressure(lineparts[8].toFloat());
		if (lineparts[11].toFloat() != -999) {
			ob.setTemperature(lineparts[11].toFloat() + 273.15);
		} else {
			ob.setTemperature(-999.);
		}
		if (lineparts[12].toFloat() != -999) {
			ob.setDewpoint(lineparts[12].toFloat() + 273.15);
		} else {
			ob.setDewpoint(-999.);
		}
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
	QString radardbz = configHash.value("radardbz");
	QString radarvel = configHash.value("radarvel");
	QString radarsw = configHash.value("radarsw");
	if(!swpfile.readSwpfile(radardbz, radarvel, radarsw))
		return false;

	float radarLat = swpfile.getRadarLat();
	float radarLon = swpfile.getRadarLon();
	float radarAlt = swpfile.getRadarAlt();
	int rayskip = configHash.value("radarskip").toInt();
	int stride = configHash.value("radarstride").toInt();
	for (int i=0; i < swpfile.getNumRays(); i+=rayskip) {
		float az = swpfile.getAzimuth(i);
		float el = swpfile.getElevation(i);
		float* refdata = swpfile.getReflectivity(i);
		float* veldata = swpfile.getRadialVelocity(i);	
		float* swdata = swpfile.getSpectrumWidth(i);
		QDateTime rayTime = swpfile.getRayTime(i);
		float* gatesp = swpfile.getGateSpacing();
		float beamwidth = sin(Pi*0.01);
		for (int n=0; n < swpfile.getNumGates()-stride; n+=stride) {
			MetObs ob;
			float range = gatesp[n+stride/2];
			stride = (range*beamwidth)/75;
			if (stride < 5) stride = 5;
			float dz = 0;
			float vr = 0;
			float sw = 0;
			float count = 0;
			for (int g=n; g<(n+stride); g++) {
				if (veldata[g] == -32768) continue;
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
	}
	
	return true;
	
}

bool VarDriver::read_sfmr(QFile& metFile, QList<MetObs>* metObVector)
{

	if (!metFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;
	
	QTextStream in(&metFile);
	QString datestr, timestr;
	QDateTime datetime;
	QTime time;
	QDate date;
	//datestr = metFile.fileName().left(8);	
	//QDate date = QDate::fromString(datestr, "yyyyMMdd");
	// Hard code date for now
	//QDate date(2003, 9, 12);
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
		ob.setAltitude(0.0);
		ob.setWindSpeed(lineparts[7].toFloat());
		ob.setObType(MetObs::sfmr);
		metObVector->push_back(ob);
	}
	
	metFile.close();
	return true;		

}

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
		float rainflag = rainprob.toFloat();
		if (rainflag > 0.50) { continue; }
		float dirdeg = dir.toFloat() + 180.;
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
		float speed = lineparts[2].toFloat() / 1.94384449;
		ob.setWindSpeed(speed);
		ob.setWindDirection(lineparts[3].toFloat());
		ob.setObType(MetObs::ascat);
		metObVector->push_back(ob);
	}
	
	metFile.close();
	return true;		
	
}

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
		// Convert pressure to altitude
		real presslevel = 100*lineparts[6].toFloat();
		// Newton's method 
		double height = 9000;
		double pmin = 1e34;
		int iter = 0;
		while ((fabs(pmin) > 1.) and (iter < 5000)) {
			double p = getReferenceVariable(pressref, height)-presslevel;
			double pprime = (getReferenceVariable(pressref, height+500.) - getReferenceVariable(pressref, height-500.))/1000.;
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

bool VarDriver::readXMLconfig(const QString& xmlfile)
{

	//QDomDocument domDoc("Configuration");
	QFile file(xmlfile);
	if (!file.open(QIODevice::ReadOnly)) {
		cout << "Error Opening Configuration File, Check Permissions on " << xmlfile.toStdString() << "\n";
		return false;
	}
	
	QString errorStr;
	int errorLine;
	int errorColumn;
	
	// Set the DOM document contents to those of the file
	if (!domDoc.setContent(&file, true, &errorStr, &errorLine, &errorColumn)) {
		QString errorReport = QString("XML Parse Error in "+xmlfile+" at Line %1, Column %2:\n%3")
		.arg(errorLine)
		.arg(errorColumn)
		.arg(errorStr);
		cout << errorReport.toStdString() << "\n";
		file.close();
		return false;
	}
	file.close();
	
	// Basic check to see if this is really a SAMURAI configuration file
	QDomElement root = domDoc.documentElement();
	if (root.tagName() != "samurai") {
		cout << "The XML file " << xmlfile.toStdString() << " is not an SAMURAI configuration file\n.";
		return false;
	}
		
	// Parse the nodes to a hash
	QDomNodeList nodeList = root.childNodes();
	for (int i = 0; i < nodeList.count(); i++) {
		QDomNode currNode = nodeList.item(i);
		QDomElement group = currNode.toElement();
		configHash.insert(group.tagName(), group.text());
	}
	
	return true;
	
}

real VarDriver::getReferenceVariable(int refVariable, real heightm)
{
	real qvbhypcoeff[5];
	real rhoacoeff[5];
	real dpdzcoeff[5];
	
	if (referenceState == jordan) {
		qvbhypcoeff[0] = 9.5108;
		qvbhypcoeff[1] = -0.002761;
		qvbhypcoeff[2] = 3.0362e-7;
		qvbhypcoeff[3] = -1.476e-11;
		qvbhypcoeff[4] = 2.6515e-16;
		
		rhoacoeff[0] = 1.1451;
		rhoacoeff[1] = -0.00010098;
		rhoacoeff[2] = 3.1546e-09;
		rhoacoeff[3] = -2.6251e-14;
		rhoacoeff[4] = -5.5556e-19;
		
		dpdzcoeff[0] = -11.432;
		dpdzcoeff[1] = 0.0010267;
		dpdzcoeff[2] = -2.7622e-08;
		dpdzcoeff[3] = -5.2937e-13;
		dpdzcoeff[4] = 3.4713e-17;
		
	}

	if (refVariable == qvbhypref) {
		real qvbhyp = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i); 
			qvbhyp += qvbhypcoeff[i] * power;
		}
		if (qvbhyp < 0.) qvbhyp = 0.;
		return qvbhyp;
	} else if (refVariable == rhoaref) {
		real rhoa = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i); 
			rhoa += rhoacoeff[i] * power;
		}
		return rhoa;
	} else if (refVariable == rhoref) {
		real rho = 0.;
		real qvbhyp = 0.;
		real rhoa = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i); 
			rhoa += rhoacoeff[i] * power;
			qvbhyp += qvbhypcoeff[i] * power;
		}
		if (qvbhyp < 0.) qvbhyp = 0.;
		real qv = bhypInvTransform(qvbhyp);
		rho = rhoa*qv/1000. + rhoa;
		return rho;
	} else if ((refVariable == href) or (refVariable == tempref) or (refVariable == pressref)) {
		// Integrate hydrostatic equation to get pressure and/or solve for T or h
		real press = 0.;
		real temp = 0.;
		real rho = 0.;
		real qvbhyp = 0.;
		real rhoa = 0.;
		for (int i = 0; i < 5; i++) {
			real power = pow(heightm, i);
			real power1 = pow(heightm, i+1);
			press += dpdzcoeff[i] * power1 / (i+1);
			rhoa += rhoacoeff[i] * power;
			qvbhyp += qvbhypcoeff[i] * power;
		}
		if (qvbhyp < 0.) qvbhyp = 0.;
		real qv = bhypInvTransform(qvbhyp);
		rho = rhoa*qv/1000. + rhoa;
		press += 101510.0;
		temp = press/(286.9*rhoa + 461.5*rhoa*qv/1000.);
		real h = 1005.7*temp + 9.81*heightm + 2.5e3*qv;
		switch (refVariable) {
			case href:
				return h;
			case tempref:
				return temp;
			case pressref:
				return press;
			default:
				break;
		}
	}
	
	return 0;
}

real VarDriver::bhypTransform(real qv)
{
	
	real qvbhyp = 0.5*((qv + 1.e-7) - 1.e-14/(qv + 1.e-7));
	return qvbhyp;
	
}

real VarDriver::bhypInvTransform(real qvbhyp)
{
	real qv = 0.;
	if (qvbhyp > 0) {
		qv = sqrt(qvbhyp*qvbhyp + 1.e-14) + qvbhyp - 1.e-7;
	}
	return qv;
}


