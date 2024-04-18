/*
 *  Vardriver.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "FileList.h"
#include "VarDriver.h"
#include "Dorade.h"
#include "ReferenceState.h"
#include "LineSplit.h"
#include <fstream>
#include <cmath>
#include <euclid/GeographicLib/TransverseMercatorExact.hpp>
#include <euclid/GeographicLib/LambertConformalConic.hpp>

#include <Radx/Radx.hh>
#include <Radx/RadxField.hh>
#include <Radx/RadxGeoref.hh>
#include <Radx/RadxRay.hh>
#include <Radx/RadxVol.hh>
#include <Radx/NcfRadxFile.hh>

#include <algorithm>
// Constructor
VarDriver::VarDriver()
{
  // Constant for all drivers
  Pi = acos(-1);

  // Set up the datatype hash
  dataSuffix["cen"] = cen;
  dataSuffix["frd"] = frd;
  dataSuffix["cls"] = cls;
  dataSuffix["1sec"] = sec;
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
  dataSuffix["rad"] = rad;
  dataSuffix["cfrad"] = cfrad;
	dataSuffix["hgt"] = terrain;
	dataSuffix["model"] = model;
	dataSuffix["rf"] = crsim;
    dataSuffix["list"] = hrdradial;
    dataSuffix["hdob"] = hdob;

  // By default we have fixed grid dimensions coming from the config file
  fixedGrid = true;
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
  std::vector<std::string> filenames = FileList(dataPath);
  std::string centerFilename;
  for (std::size_t i = 0; i < filenames.size(); ++i) {
    std::string filename = filenames[i];
		if (filename.empty()) {
			continue;
		}
    std::string suffix = Extension(filename);
		suffix = suffix.substr(1);
    if (suffix == "cen") {
      // Match to centerfile
      centerFilename = filenames[i];
      break;
    }
  }

  // Open the file
	std::ifstream centerFile(dataPath + "/" + centerFilename, std::ifstream::in);
	if (!centerFile.is_open()) {
    std::cout << "Unable to open centerfile " << centerFilename << std::endl;
    return false;
  }

  // Get the date from the filename
  std::string datestr = centerFilename.substr(0,8);
  datetime startDate = ParseDate(datestr.c_str(), "%Y%m%d");

  std::cout << "DEBUG: center date = " << PrintDate(startDate) << ", time = " << Time(startDate) << std::endl;

  // Read the centers
  std::string line;
	float lat, lon, Vm, Um;
  while (std::getline(centerFile, line)) {
		std::istringstream iss(line);
    datetime date;
		std::string timestr;
 		iss >> timestr;
    int hours = std::stoi(timestr.substr(0,2));
    if (hours > 23) { // (FIXME : NCAR) The original code added a day here, then subtracted 24 from the hours, is this needed?
      //fixme date = date::make_zoned(startDate + date::days{1});
      // hours -= 24;
    } else {
      date = startDate;
    }
		int minutes = std::stoi(timestr.substr(2,2));
		int seconds = std::stoi(timestr.substr(4,2));
    datetime datetime_ = startDate + std::chrono::hours{hours} + std::chrono::minutes{minutes} + std::chrono::seconds{seconds};
		iss >> lat;
		iss >> lon;
		iss >> Vm;
		iss >> Um;
    FrameCenter center(datetime_, lat, lon, Um, Vm);
    frameVector.push_back(center);
  }

  return true;
}

// These are interface functions to allow caller manipulation of the center vector

// Empty the centers vector

void VarDriver::clearCenters()
{
  frameVector.clear();
}

// Append a new center to the centers vector

void VarDriver::appendCenter(std::string date, std::string time, float lat, float lon, float Vm, float Um)
{
  std::string in = date + time;
  datetime datetime_ = ParseDate(in.c_str(), "%Y%m%d%H%M%S");
  frameVector.push_back(FrameCenter(datetime_, lat, lon, Um, Vm));
}

// Pop the first center (Might add a date later if needed to remove a specific center)

void VarDriver::popCenter()
{
  frameVector.erase(frameVector.begin());
}

// Pop the first center (Might add a date later if needed to remove a specific center)

// bool VarDriver::readTerrain(std::string &filename, std::vector<MetObs>* metData)
// {
// 	std::string metFile = filename;
//   return read_terrain(metFile, metData);
// }

// interface function to read radar data

bool VarDriver::read_met_obs_file(int suffix, std::string &filename, std::vector<MetObs>* metData)
{

  std::string metFile = filename;
  switch (suffix) {
  case (frd):
    if (!read_frd(metFile, metData)) {
      cout << "Error reading frd file" << endl;
      return false;
    }
    break;
  case (cls):
    if (!read_cls(metFile, metData)) {
      cout << "Error reading cls file" << endl;
      return false;
    }
    break;
  case (sec):
    if (!read_sec(metFile, metData)) {
      cout << "Error reading 1sec file" << endl;
      return false;
    }
    break;
  case (ten):
    if (!read_ten(metFile, metData)) {
      cout << "Error reading ten file" << endl;
      return false;
    }
    break;
  case (swp):
    if (!read_dorade(metFile, metData)) {
      cout << "Error reading swp file" << endl;
      return false;
    }
    break;
  case (sfmr):
    if (!read_sfmr(metFile, metData)) {
      cout << "Error reading sfmr file" << endl;
      return false;
    }
    break;
  case (wwind):
    if (!read_wwind(metFile, metData)) {
      cout << "Error reading wwind file" << endl;
      return false;
    }
    break;
  case (eol):
    if (!read_eol(metFile, metData)) {
      cout << "Error reading eol file" << endl;
      return false;
    }
    break;
  case (qscat):
    if (!read_qscat(metFile, metData)) {
      cout << "Error reading wwind file" << endl;
      return false;
    }
    break;
  case (ascat):
    if (!read_ascat(metFile, metData)) {
      cout << "Error reading wwind file" << endl;
      return false;
    }
    break;
  case (nopp):
    if (!read_nopp(metFile, metData)) {
      cout << "Error reading wwind file" << endl;
      return false;
    }
    break;
  case (cimss):
    if (!read_cimss(metFile, metData)) {
      cout << "Error reading cimss file" << endl;
      return false;
    }
    break;
  case (dwl):
    if (!read_dwl(metFile, metData)) {
      cout << "Error reading dwl file" << endl;
      return false;
    }
    break;
  case (insitu):
    if (!read_insitu(metFile, metData)) {
      cout << "Error reading insitu file" << endl;
      return false;
    }
    break;
	case (terrain):
		if (!read_terrain(metFile, metData)) {
			cout << "Error reading terrain file" << endl;
			return false;
		}
		break;
	case (model):
		if (!read_model(metFile, metData)) {
			cout << "Error reading model file" << endl;
			return false;
		}
		break;
  case (mtp):
    if (!read_mtp(metFile, metData)) {
      cout << "Error reading mtp file" << endl;
      return false;
    }
  case (mesonet):
    if (!read_mesonet(metFile, metData)) {
      cout << "Error reading mesonet file" << endl;
      return false;
    }
    break;
  case (classnc):
    if (!read_classnc(metFile, metData)) {
      cout << "Error reading classnc file" << endl;
      return false;
    }
    break;
  case (qcf):
    if (!read_qcf(metFile, metData)) {
      cout << "Error reading classnc file" << endl;
      return false;
    }
    break;
  case(aeri):
    if (!read_aeri(metFile, metData)) {
      cout << "Error reading aeri file" << endl;
      return false;
    }
    break;
  case(rad):
    if (!read_rad(metFile, metData)) {
      cout << "Error reading rad file" << endl;
      return false;
    }
    break;
  case (cen):
    return false;
  case(cfrad):
    if (!read_cfrad(metFile, metData)) {
      cout << "Error reading cfrad file" << endl;
      return false;
    }
    break;
	case(crsim):
	  if (!read_crsim(metFile, metData)) {
	    cout << "Error reading rf file" << endl;
	    return false;
	  }
	    break;
  case(hrdradial):
    if (!read_hrdradial(metFile, metData)) {
      cout << "Error reading hrd radial file" << endl;
      return false;
    }
      break;
  case(hdob):
    if (!read_hdobs(metFile, metData)) {
      cout << "Error reading hdob file" << endl;
      return false;
    }
      break;
  default:
    cout << "Unknown data type, skipping..." << endl;
    return false;
  }

  return true;
}

/* This routine reads the FRD insitu format from NOAA/HRD */

bool VarDriver::read_frd(std::string& filename, std::vector<MetObs>* metObVector)
{
	std::ifstream metFile(filename);
  if (!metFile.is_open())
    return false;
/*
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
*/
  metFile.close();
  return true;
}

/* This routine reads the old CLASS dropsonde format from NCAR */

bool VarDriver::read_cls(std::string& filename, std::vector<MetObs>* metObVector)
{
	std::ifstream metFile(filename);
  if (!metFile.is_open())
    return false;

/*
  QTextStream in(&metFile);
  QString datestr, timestr, aircraft;
  QDateTime datetime;
  bool start = false;
  while (!in.atEnd()) {
    QString line = in.readLine();
    if (line.startsWith("Release Site Type")) {
      QStringList lineparts = line.split(":");
      aircraft = lineparts[1].trimmed();
      // } else if (line.startsWith("UTC") || line.startsWith("GMT")) {
    } else if (line.startsWith("UTC")) {
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
      int sec = lineparts[1].toFloat();
      ob.setTime(datetime.addSecs(sec));
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
*/
  metFile.close();
  return true;
}

/* This routine reads the modified CLASS dropsonde format from NCAR for the TPARC/TCS08 Dataset (and maybe TREX?) */

bool VarDriver::read_wwind(std::string& filename, std::vector<MetObs>* metObVector)
{
	std::ifstream metFile(filename);
  if (!metFile.is_open())
    return false;

/*

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
*/
  metFile.close();
  return true;
}

/* This routine reads the modified CLASS dropsonde format from NCAR for the PREDICT field project */

bool VarDriver::read_eol(std::string& filename, std::vector<MetObs>* metObVector)
{
	std::ifstream metFile(filename);
  if (!metFile.is_open())
    return false;

  std::string datestr, timestr, aircraft;
  datetime datetime_;
  bool start = false;

  std::string line;


  while (std::getline(metFile, line)) {
           //std::cout << "DEBUG: line = " << line << std::endl;
    if (line.find("Launch Site Type") == 0) {  // find of 0 means it's there, starting at position 0 (vs. npos)
      auto lineparts = LineSplit(line, ':');
      aircraft = lineparts[1];
    } else if (line.find("UTC") == 0) {
      datestr = line.substr(43,12);
      timestr = line.substr(57,8);

      std::string combined = datestr + " " + timestr;
      datetime_ = ParseDate(combined.c_str(), "%Y, %m, %d %H:%M:%S");
    } else if (line.find("------") == 0) {
      // Start reading data
      start = true;
    } else if (start) {
      MetObs ob;
      ob.setStationName(aircraft);
      //line = QString(" ") + line;
      line = "Blank " + line; // The blank is to force the +1 indexing here; original code used a space (see above)
      auto lineparts = LineSplit(line, ' ');
      int sec = (int)(std::stof(lineparts[1]));
      ob.setTime(datetime_ + std::chrono::seconds(sec));
      if (std::stof(lineparts[15]) != -999.) {
	      ob.setLon(std::stof(lineparts[15]));
      } else {
	      ob.setLon(-999.);
      }
      if (std::stof(lineparts[16]) != -999.) {
	      ob.setLat(std::stof(lineparts[16]));
      } else {
	      ob.setLat(-999.);
      }
      if (std::stof(lineparts[14]) != -999.0) {
	      ob.setAltitude(std::stof(lineparts[14]));
      } else {
	      ob.setAltitude(-999.);
      }
      if (std::stof(lineparts[5]) != -999.0) {
	      ob.setPressure(std::stof(lineparts[5]));
      } else {
	      ob.setPressure(-999.0);
      }
      if (std::stof(lineparts[6]) != -999.0) {
	      ob.setTemperature(std::stof(lineparts[6]) + 273.15);
      } else {
	      ob.setTemperature(-999.);
      }
      if (std::stof(lineparts[8]) != -999.0) {
	      ob.setRH(std::stof(lineparts[8]));
      } else {
	      ob.setRH(-999.);
      }
      if (std::stof(lineparts[12]) != -999.0) {
	      ob.setWindDirection(std::stof(lineparts[12]));
      } else {
	      ob.setWindDirection(-999.);
      }
      if (std::stof(lineparts[11]) != -999.0) {
	      ob.setWindSpeed(std::stof(lineparts[11]));
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

bool VarDriver::read_sec(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);
  if (!metFile.is_open()){
    return false;}
// Skip four lines
  std::string line;
  std::getline(metFile, line);
  auto parts = LineSplit(line, ' ');
  std::string datestr = parts[1].substr(0,8);
  datetime date = ParseDate(datestr.c_str(), "%Y%m%d") ;
  std::getline(metFile, line);std::getline(metFile, line);std::getline(metFile, line);

  MetObs ob;
  while (std::getline(metFile, line)) {
      auto parts = LineSplit(line, ' ');
      ob.setStationName("aircraft");
      std::string timestr;
      int hours = std::stoi(parts[0].substr(0,2));
      if (hours > 23) { // (FIXME : NCAR) The original code added a day here, then subtracted 24 from the hours, is this needed?
		    date += date::days{1};
		    hours -= 24;
      }
      int minutes = std::stoi(parts[0].substr(2,2));
	    int seconds = std::stoi(parts[0].substr(4,2));
      datetime datetime_ = date + std::chrono::hours{hours} + std::chrono::minutes{minutes} + std::chrono::seconds{seconds};
      std::cout << "Date: " << PrintDate(datetime_) << std::endl;
      ob.setTime(datetime_);
      // ob.setStationName(aircraft);
      ob.setLat(std::stof(parts[1]));
      // Note that the longitude is given in degrees West,
      ob.setLon(-std::stof(parts[2]));
      ob.setAltitude(std::stof(parts[7]));
      ob.setPressure(std::stof(parts[8]));
      if ((std::stof(parts[11]) != -999) or (std::stof(parts[11]) != -32767)) {
        ob.setTemperature(std::stof(parts[11]) + 273.15);
      } else {
        ob.setTemperature(-999.);
      }
      if ((std::stof(parts[12]) != -999) or (std::stof(parts[12]) != -32767)) {
        ob.setDewpoint(std::stof(parts[12]) + 273.15);
      } else {
        ob.setDewpoint(-999.);
      }
      if (std::stof(parts[10]) >= 0) {
        ob.setWindDirection(std::stof(parts[9]));
        ob.setWindSpeed(std::stof(parts[10]));
      }
      if ((std::stof(parts[16]) != -999) or (std::stof(parts[16]) != -32767)) {
        ob.setVerticalVelocity(std::stof(parts[16]));
      }
      ob.setObType(MetObs::flightlevel);
      metObVector->push_back(ob);
    }
    metFile.close();
    return true;

}
  

/*
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
*/
//   metFile.close();
//   return true;

// }

/* This routine reads a 10 sec flight level data file from the data provided by the USAF Hurricane Hunters */
/* GMT Time   AOA     BSP   CAS    CC     CSP    DPR     DVAL      GA     GPSA   GS      HSS     LAT      LON       PA PITCH       RA ROLL     SLP   SS     TA   TAS    TDA    TDD  THD   TRK    TT   V V     WD     WS Valid Flags Source Tags */

bool VarDriver::read_ten(std::string& filename, std::vector<MetObs>* metObVector)
{
	std::ifstream metFile(filename);
  if (!metFile.is_open())
    return false;

/*
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
*/
  metFile.close();
  return true;
}

/* This routine reads a Dorade Sweepfile
   This does not correctly read 'pure' Dorade files due to big-endian representation in that format,
   so it needs a byte swap which has not been implemented yet*/

bool VarDriver::read_dorade(std::string& filename, std::vector<MetObs>* metObVector)
{
  Dorade swpfile(filename);

  // Use a Transverse Mercator projection to map the radar gates to the grid
  //Geographiclib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();

  std::string radardbz = configHash["radar_dbz"]; // If radar_dbz isn't found, this will create it -- switch to find/check? (bdobbins)
  std::string radarvel = configHash["radar_vel"]; // See above (bdobbins)
  std::string radarsw = configHash["radar_sw"];  // See above (bdobbins)
  if(!swpfile.readSwpfile(radardbz, radarvel, radarsw))
    return false;

  int rayskip = std::stoi(configHash["radar_skip"]);
  int minstride = std::stoi(configHash["radar_stride"]);
  bool dynamicStride = std::stoi(configHash["dynamic_stride"]);
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
    datetime rayTime = swpfile.getRayTime(i);
    float* gatesp = swpfile.getGateSpacing();
    real gatelength = gatesp[1] - gatesp[0];
    real beamwidth = sin(swpfile.getBeamwidthDeg()*Pi/180.);

    for (int n=0; n < swpfile.getNumGates()-stride; n+=stride) {
      MetObs ob;
      real range = gatesp[n+stride/2];
      if (dynamicStride) {
	stride = (int)(range * beamwidth / gatelength);
	if (stride < minstride) stride = minstride;
      }
      real dz = 0.0;
      real vr = 0.0;
      real sw = 0.0;
      int dzcount = 0;
      int vrcount = 0;
      int swcount = 0;
      for (int g=n; g<(n+stride); g++) {
	if (gatesp[g] <= 0) continue;
	if (veldata[g] != -32768) {
	  vr += veldata[g];
	  vrcount++;
	}
	if (refdata[g] != -32768) {
	  dz += pow(10.0,(refdata[g]*0.1));
	  dzcount++;
	}
	if (swdata[g] != -32768) {
	  sw += swdata[g];
	  swcount++;
	}
      }
      if (dzcount > 0) {
	dz = dz/float(dzcount);
	dz = 10*log10(dz);
      } else {
	dz = -999.0;
      }
      if (vrcount > 0) {
	vr = vr/float(vrcount);
      } else {
	vr = -999.0;
      }
      if (swcount > 0) {
	sw = sw/float(swcount);
      } else {
	sw = -999.0;
      }
      if ((vr != -999.0) || (dz != -999.0)) {
	real relX = range*sin(az*Pi/180.)*cos(el*Pi/180.);
	real relY = range*cos(az*Pi/180.)*cos(el*Pi/180.);
	real rEarth = 6371000.0;

	// Take into account curvature of the earth for the height of the radar beam
	real relZ = sqrt(range*range + rEarth*rEarth + 2.0 * range * rEarth * sin(el*Pi/180.)) - rEarth;
	real radarX, radarY, gateLat, gateLon;
	projection.Forward(radarLon, radarLat, radarLon, radarX, radarY);
	projection.Reverse(radarLon, radarX + relX, radarY + relY, gateLat, gateLon);
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

bool VarDriver::read_sfmr(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);
  if (!metFile.is_open()) {
    return false;
	}
/*
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
    //cout << datestr.toLatin1().data() << endl;
    ob.setLat(lineparts[5].toFloat());
    ob.setLon(lineparts[6].toFloat());
    ob.setAltitude(10.0);
    ob.setWindSpeed(lineparts[7].toFloat());
    ob.setObType(MetObs::sfmr);
    metObVector->push_back(ob);
  }
*/

  metFile.close();
  return true;

}

/* This routine reads an ASCII dump of the Level II Quikscat data, see code for columns*/

bool VarDriver::read_qscat(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);
  if (!metFile.is_open()) {
    return false;
	}
/* fixme
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
*/

  metFile.close();
  return true;

}

/* This routine reads an ASCII dump of the ASCAT data, see code for columns*/

bool VarDriver::read_ascat(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);
  if (!metFile.is_open()) {
    return false;
	}

/*fixme
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
*/
  metFile.close();
  return true;

}

/* This routine reads an ASCII dump of the Quikscat data for NOPP, see code for columns*/

bool VarDriver::read_nopp(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);

  if (!metFile.is_open()) {
    return false;
	}
/*fixme
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
*/
  metFile.close();
  return true;

}

/* This routine reads the CIMSS atmospheric motion vector files.
   Pressure levels are converted to Heights via Newton's method and the background
   reference state */

bool VarDriver::read_cimss(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);

  if (!metFile.is_open()) {
    return false;
	}

  MetObs ob;
  datetime datetime_;
  // Skip first line
  std::string line;
  std::getline(metFile, line);

  while (std::getline(metFile, line)) {
    auto parts = LineSplit(line, ' ');
    ob.setStationName(parts[0] + " " + parts[1]);
    datetime datetime_ = ParseDate(parts[2] + parts[3], "%Y%m%d%H%M");
    //std::cout << "Date: " << PrintDate(datetime_) << std::endl;
    ob.setTime(datetime_);
    ob.setLat(std::stof(parts[4]));
    ob.setLon(-std::stof(parts[5]));
    //std::cout << "DEBUG: Lat/Lon : " << parts[4] << " / " << parts[5] << std::endl;
    real presslevel = 100*std::stof(parts[6]);
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
    ob.setWindSpeed(std::stof(parts[7]));
    ob.setWindDirection(std::stof(parts[8]));
    ob.setObType(MetObs::AMV);
    metObVector->push_back(ob);
  }
  metFile.close();
  return true;

}

/* This routine reads the LOS Doppler Wind Lidar data from TPARC, see code for columns*/

bool VarDriver::read_dwl(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);

  if (!metFile.is_open()) {
    return false;
	}
/*fixme
  //GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();

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
	  projection.Forward(radarLon, radarLat, radarLon, radarX, radarY);
	  projection.Reverse(radarLon, radarX + relX, radarY + relY, gateLat, gateLon);
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
	  //cout << datetime.toString(Qt::ISODate).toStdString() << "\t"
	  // << gateLat << "\t" << gateLon << "\t" << gateAlt << "\t"
	  //  << az << "\t" << el << "\t" << dz << "\t" << vr << "\t" << sw << endl;
	}
      }
    }
  }
*/
  return true;

}

/* This routine reads a generic insitu format suitable for many different observation types
   TIME       Lat       Lon    Altitude   Pressure  Temperature (C)  Dewpoint (C)  WindDir   WindSpeed  VerticalVelocity*/
bool VarDriver::read_insitu(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);

  if (!metFile.is_open()) {
    return false;
	}
/*fixme
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
*/
  metFile.close();
  return true;
}

bool VarDriver::read_mtp(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);

  if (!metFile.is_open()) {
    return false;
	}
/*fixme
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
*/
  metFile.close();
  return true;
}

/* This routine reads the Terrain text file.
   Pressure levels are converted to Heights via Newton's method and the background
   reference state */

bool VarDriver::read_terrain(std::string& filename, std::vector<MetObs>* metObVector)
{

  std::ifstream metFile(filename);
  if (!metFile.is_open()) {
    return false;
	}

	std::vector<std::string> filenames = FileList(dataPath);
	std::string centerFilename;
	for (std::size_t i = 0; i < filenames.size(); ++i) {
		std::string filename = filenames[i];
		if (filename.empty()) {
			continue;
		}
		std::string suffix = Extension(filename);
		suffix = suffix.substr(1);
		if (suffix == "cen") {
			// Match to centerfile
			centerFilename = filenames[i];
			break;
	}}
	// Open the file
	std::ifstream centerFile(dataPath + "/" + centerFilename, std::ifstream::in);
	// Get the date from the filename
	std::string datestr = centerFilename.substr(0,8);
	datetime startDate = ParseDate(datestr.c_str(), "%Y%m%d");

	std::string line;
	std::getline(centerFile, line);
	std::istringstream iss(line);
	datetime date;
	std::string timestr;
	iss >> timestr;
	int hours = std::stoi(timestr.substr(0,2));
	if (hours > 23) { // (FIXME : NCAR) The original code added a day here, then subtracted 24 from the hours, is this needed?
		//fixme date = date::make_zoned(startDate + date::days{1});
		// hours -= 24;
	} else {
		date = startDate;
	}
	int minutes = std::stoi(timestr.substr(2,2));
	int seconds = std::stoi(timestr.substr(4,2));
	datetime datetime_ = startDate + std::chrono::hours{hours} + std::chrono::minutes{minutes} + std::chrono::seconds{seconds};

  MetObs ob;
  std::string line2;

  while (std::getline(metFile, line2)) {
		ob.setTime(datetime_);
    auto parts = LineSplit(line2, ' ');
    ob.setLat(std::stof(parts[0]));
    ob.setLon(std::stof(parts[1]));
    ob.setAltitude(std::stof(parts[2]));
    ob.setTerrainDX(std::stof(parts[3]));
    ob.setTerrainDY(std::stof(parts[4]));
    ob.setTerrainX(std::stof(parts[5]));
    ob.setTerrainY(std::stof(parts[6]));
    ob.setObType(MetObs::terrain);
    metObVector->push_back(ob);
  }
	// std::cout << "Successfully read the terrain file" << std::endl;
  metFile.close();
  return true;

}

/* This routine reads a text file from WRF or other model output as pseudo-observations*/

bool VarDriver::read_model(std::string& filename, std::vector<MetObs>* metObVector)
{
	std::ifstream metFile(filename);
  if (!metFile.is_open()) {
    return false;
	}
	MetObs ob;
	datetime datetime_;
	std::string line;

  while (std::getline(metFile, line)) {
		auto parts = LineSplit(line, ' ');
		datetime datetime_ = ParseTime(parts[0].c_str(), "%Y-%m-%d_%H:%M:%S");
		ob.setTime(datetime_);
    ob.setLat(std::stof(parts[1]));
    ob.setLon(std::stof(parts[2]));
    ob.setAltitude(std::stof(parts[3]));
    ob.setZonalVelocity(std::stof(parts[4]));
    ob.setMeridionalVelocity(std::stof(parts[5]));
		ob.setVerticalVelocity(std::stof(parts[6]));
		ob.setTemperature(std::stof(parts[7]));
		ob.setModelQv(std::stof(parts[8]));
		ob.setModelAirDensity(std::stof(parts[9]));
		ob.setModelMoistDensity(std::stof(parts[10]));
		ob.setReflectivity(std::stof(parts[11]));
    ob.setObType(MetObs::model);
    metObVector->push_back(ob);
  }
	std::cout << "Successfully read the model file" << std::endl;
  metFile.close();
  return true;

}

/* This routine parses the supplied XML configuration
   and validates parameters that are common to all drivers
   Mode specific configs are validated in those drivers */

bool VarDriver::parseXMLconfig(const XMLNode& config)
{
  // Parse the nodes to a hash:
  // 1) Get the list (names) of elements for this node:
  std::vector<const XMLElement* > elementList = XMLGetElements(&config);

	// 2) Loop over all the entries in the node list:
	for (const auto element : elementList) {

		// 3) Assign the default 'iter' value, and check for the attribute:
		std::string iter = std::to_string(0);

		// 4) For each top-level node name, get its children:
		std::vector<const XMLAttribute *> attributeList = XMLGetAttributes(element);

		// 5) Loop over the children and check for 'iter':
		for (auto &attribute : attributeList) {
			if (std::string(attribute->Name()) == "iter" ) {
				iter = attribute->Value();
			}
		}

		// 6) Build the list of children (config options):
		std::vector<const XMLElement*> configList = XMLGetElements(element);

		// 7) Loop over the children and if 'iter' > 1 then create hash with iter appended:
		for (auto &config : configList) {
			std::string tag(config->Name());
			if (std::stoi(iter) > 1) {
				tag +=  "_" + iter;
			}
			configHash.insert(tag, config->GetText());
		}
	}

  // Validate the hash -- multiple passes are not validated currently
  std::set<std::string> configKeys;
	configKeys.insert("ref_state");
	configKeys.insert("qr_variable");
	configKeys.insert("i_background_roi");
	configKeys.insert("j_background_roi");
	configKeys.insert("i_reflectivity_roi");
	configKeys.insert("j_reflectivity_roi");
	configKeys.insert("k_reflectivity_roi");
	configKeys.insert("load_background");
	configKeys.insert("adjust_background");
	configKeys.insert("radar_dbz");
	configKeys.insert("radar_vel");
	configKeys.insert("radar_sw");
	configKeys.insert("radar_skip");
	configKeys.insert("radar_stride");
	configKeys.insert("dynamic_stride");
	configKeys.insert("dbz_pseudow_weight");
	configKeys.insert("melting_zone_width");
	configKeys.insert("mixed_phase_dbz");
	configKeys.insert("rain_dbz");
	configKeys.insert("num_iterations");
	configKeys.insert("output_mish");
	configKeys.insert("output_txt");
	configKeys.insert("output_qc");
	configKeys.insert("output_netcdf");
	configKeys.insert("output_asi");
	configKeys.insert("preprocess_obs");
	configKeys.insert("mask_reflectivity");
	configKeys.insert("dropsonde_rhou_error");
	configKeys.insert("dropsonde_rhov_error");
	configKeys.insert("dropsonde_rhow_error");
	configKeys.insert("dropsonde_tempk_error");
	configKeys.insert("dropsonde_qv_error");
	configKeys.insert("dropsonde_rhoa_error");
	configKeys.insert("flightlevel_rhou_error");
	configKeys.insert("flightlevel_rhov_error");
	configKeys.insert("flightlevel_rhow_error");
	configKeys.insert("flightlevel_tempk_error");
	configKeys.insert("flightlevel_qv_error");
	configKeys.insert("flightlevel_rhoa_error");
	configKeys.insert("insitu_rhou_error");
	configKeys.insert("insitu_rhov_error");
	configKeys.insert("insitu_rhow_error");
	configKeys.insert("insitu_tempk_error");
	configKeys.insert("insitu_qv_error");
	configKeys.insert("insitu_rhoa_error");
	configKeys.insert("sfmr_windspeed_error");
	configKeys.insert("qscat_rhou_error");
	configKeys.insert("qscat_rhov_error");
	configKeys.insert("ascat_rhou_error");
	configKeys.insert("ascat_rhov_error");
	configKeys.insert("amv_rhou_error");
	configKeys.insert("amv_rhov_error");
	configKeys.insert("lidar_sw_error");
	configKeys.insert("lidar_power_error");
	configKeys.insert("lidar_min_error");
	configKeys.insert("radar_sw_error");
	configKeys.insert("radar_fallspeed_error");
	configKeys.insert("radar_min_error");
	configKeys.insert("bg_rhou_error");
	configKeys.insert("bg_rhov_error");
	configKeys.insert("bg_rhow_error");
	configKeys.insert("bg_tempk_error");
	configKeys.insert("bg_qv_error");
	configKeys.insert("bg_rhoa_error");
	configKeys.insert("bg_qr_error");

  if (fixedGrid) {
    configKeys.insert("ref_time");
	}

  for (auto& key : configKeys) {
  if (configHash.exists(key) == false) {
      cout <<	"No configuration found for <" << key << "> aborting..." << endl;
      return false;
    }
  }

  // Set default for Bkgd Obs reader (only used if passed as arrays)
	if (configHash.exists("array_order") == false) {
  	configHash.insert("array_order", "column-major");
	}

  // Make sure we have a valid value for "array_order"
  if ( (configHash["array_order"] != "column-major") && (configHash["array_order"] != "row-major")) {
    cout << "Unsupported array_order" << configHash["array_order"] << ". Only column-major and row-major are supported" << ". Aborting..." << endl;
    return false;
  }

  return true;
}

// This routine does the same thing as parseXMLconfig, excepts from a samurai_config structure instead of a QDomElement

bool VarDriver::parseSamuraiConfig(const samurai_config &config)
{
  std::cout << "Parsing COAMPS structure ...\n";

  // Parse the values to a hash
  std::string tmpstr;
/*fixme
  configHash.insert("num_iterations",tmpstr.setNum(config.num_iterations));
  configHash.insert("radar_skip",tmpstr.setNum(config.radar_skip));
  configHash.insert("radar_stride",tmpstr.setNum(config.radar_stride));
  configHash.insert("dynamic_stride",tmpstr.setNum(config.dynamic_stride));
  configHash.insert("spline_approximation",tmpstr.setNum(config.spline_approximation));
*/
#if 0 /* old style call */
  configHash.insert("nx",tmpstr.setNum(config.nx));
  configHash.insert("ny",tmpstr.setNum(config.ny));
  configHash.insert("nz",tmpstr.setNum(config.nz));

  configHash.insert("i_min",tmpstr.setNum(config.i_min));
  configHash.insert("i_max",tmpstr.setNum(config.i_max));
  configHash.insert("i_incr",tmpstr.setNum(config.i_incr));
  configHash.insert("j_min",tmpstr.setNum(config.j_min));
  configHash.insert("j_max",tmpstr.setNum(config.j_max));
  configHash.insert("j_incr",tmpstr.setNum(config.j_incr));
  configHash.insert("k_min",tmpstr.setNum(config.k_min));
  configHash.insert("k_max",tmpstr.setNum(config.k_max));
  configHash.insert("k_incr",tmpstr.setNum(config.k_incr));
#endif
 /*fixme
  configHash.insert("i_background_roi",tmpstr.setNum(config.i_background_roi));
  configHash.insert("j_background_roi",tmpstr.setNum(config.j_background_roi));
  configHash.insert("i_reflectivity_roi",tmpstr.setNum(config.i_reflectivity_roi));
  configHash.insert("j_reflectivity_roi",tmpstr.setNum(config.j_reflectivity_roi));
  configHash.insert("k_reflectivity_roi",tmpstr.setNum(config.k_reflectivity_roi));
  configHash.insert("dbz_pseudow_weight",tmpstr.setNum(config.dbz_pseudow_weight));
  configHash.insert("melting_zone_width",tmpstr.setNum(config.melting_zone_width));
  configHash.insert("mixed_phase_dbz",tmpstr.setNum(config.mixed_phase_dbz));
  configHash.insert("rain_dbz",tmpstr.setNum(config.rain_dbz));
  configHash.insert("bg_rhou_error",tmpstr.setNum(config.bg_rhou_error));
  configHash.insert("bg_rhov_error",tmpstr.setNum(config.bg_rhov_error));
  configHash.insert("bg_rhow_error",tmpstr.setNum(config.bg_rhow_error));
  configHash.insert("bg_tempk_error",tmpstr.setNum(config.bg_tempk_error));
  configHash.insert("bg_qv_error",tmpstr.setNum(config.bg_qv_error));
  configHash.insert("bg_rhoa_error",tmpstr.setNum(config.bg_rhoa_error));
  configHash.insert("bg_qr_error",tmpstr.setNum(config.bg_qr_error));
  configHash.insert("mc_weight",tmpstr.setNum(config.mc_weight));
  configHash.insert("i_filter_length",tmpstr.setNum(config.i_filter_length));
  configHash.insert("j_filter_length",tmpstr.setNum(config.j_filter_length));
  configHash.insert("k_filter_length",tmpstr.setNum(config.k_filter_length));
  configHash.insert("i_spline_cutoff",tmpstr.setNum(config.i_spline_cutoff));
  configHash.insert("j_spline_cutoff",tmpstr.setNum(config.j_spline_cutoff));
  configHash.insert("k_spline_cutoff",tmpstr.setNum(config.k_spline_cutoff));
  configHash.insert("i_max_wavenumber",tmpstr.setNum(config.i_max_wavenumber));
  configHash.insert("j_max_wavenumber",tmpstr.setNum(config.j_max_wavenumber));
  configHash.insert("k_max_wavenumber",tmpstr.setNum(config.k_max_wavenumber));
  configHash.insert("dropsonde_rhou_error",tmpstr.setNum(config.dropsonde_rhou_error));
  configHash.insert("dropsonde_rhov_error",tmpstr.setNum(config.dropsonde_rhov_error));
  configHash.insert("dropsonde_rhow_error",tmpstr.setNum(config.dropsonde_rhow_error));
  configHash.insert("dropsonde_tempk_error",tmpstr.setNum(config.dropsonde_tempk_error));
  configHash.insert("dropsonde_qv_error",tmpstr.setNum(config.dropsonde_qv_error));
  configHash.insert("dropsonde_rhoa_error",tmpstr.setNum(config.dropsonde_rhoa_error));
  configHash.insert("flightlevel_rhou_error",tmpstr.setNum(config.flightlevel_rhou_error));
  configHash.insert("flightlevel_rhov_error",tmpstr.setNum(config.flightlevel_rhov_error));
  configHash.insert("flightlevel_rhow_error",tmpstr.setNum(config.flightlevel_rhow_error));
  configHash.insert("flightlevel_tempk_error",tmpstr.setNum(config.flightlevel_tempk_error));
  configHash.insert("flightlevel_qv_error",tmpstr.setNum(config.flightlevel_qv_error));
  configHash.insert("flightlevel_rhoa_error",tmpstr.setNum(config.flightlevel_rhoa_error));
  configHash.insert("insitu_rhou_error",tmpstr.setNum(config.insitu_rhou_error));
  configHash.insert("insitu_rhov_error",tmpstr.setNum(config.insitu_rhov_error));
  configHash.insert("insitu_rhow_error",tmpstr.setNum(config.insitu_rhow_error));
  configHash.insert("insitu_tempk_error",tmpstr.setNum(config.insitu_tempk_error));
  configHash.insert("insitu_qv_error",tmpstr.setNum(config.insitu_qv_error));
  configHash.insert("insitu_rhoa_error",tmpstr.setNum(config.insitu_rhoa_error));
  configHash.insert("sfmr_windspeed_error",tmpstr.setNum(config.sfmr_windspeed_error));
  configHash.insert("qscat_rhou_error",tmpstr.setNum(config.qscat_rhou_error));
  configHash.insert("qscat_rhov_error",tmpstr.setNum(config.qscat_rhov_error));
  configHash.insert("ascat_rhou_error",tmpstr.setNum(config.ascat_rhou_error));
  configHash.insert("ascat_rhov_error",tmpstr.setNum(config.ascat_rhov_error));
  configHash.insert("amv_rhou_error",tmpstr.setNum(config.amv_rhou_error));
  configHash.insert("amv_rhov_error",tmpstr.setNum(config.amv_rhov_error));
  configHash.insert("lidar_sw_error",tmpstr.setNum(config.lidar_sw_error));
  configHash.insert("lidar_power_error",tmpstr.setNum(config.lidar_power_error));
  configHash.insert("lidar_min_error",tmpstr.setNum(config.lidar_min_error));
  configHash.insert("radar_sw_error",tmpstr.setNum(config.radar_sw_error));
  configHash.insert("radar_fallspeed_error",tmpstr.setNum(config.radar_fallspeed_error));
  configHash.insert("radar_min_error",tmpstr.setNum(config.radar_min_error));
*/
#if 0 /* old style */
  configHash.insert("delx",tmpstr.setNum(config.delx));
  configHash.insert("dely",tmpstr.setNum(config.dely));
#endif

/*fixme
  if (config.load_background) configHash.insert("load_background","true");
  else configHash.insert("load_background","false");
  if (config.adjust_background) configHash.insert("adjust_background","true");
  else configHash.insert("adjust_background","false");
  if (config.preprocess_obs) configHash.insert("preprocess_obs","true");
  else configHash.insert("preprocess_obs","false");
  if (config.output_mish) configHash.insert("output_mish","true");
  else configHash.insert("output_mish","false");
  if (config.output_txt) configHash.insert("output_txt","true");
  else configHash.insert("output_txt","false");
  if (config.output_qc) configHash.insert("output_qc","true");
  else configHash.insert("output_qc","false");
  if (config.output_netcdf) configHash.insert("output_netcdf","true");
  else configHash.insert("output_netcdf","false");
  if (config.output_asi) configHash.insert("output_asi","true");
  else configHash.insert("output_asi","false");
  if (config.output_COAMPS) configHash.insert("output_COAMPS","true");
  else configHash.insert("output_COAMPS","false");

  configHash.insert("i_rhou_bcL",config.i_rhou_bcL);
  configHash.insert("i_rhou_bcR",config.i_rhou_bcR);
  configHash.insert("i_rhov_bcL",config.i_rhov_bcL);
  configHash.insert("i_rhov_bcR",config.i_rhov_bcR);
  configHash.insert("i_rhow_bcL",config.i_rhow_bcL);
  configHash.insert("i_rhow_bcR",config.i_rhow_bcR);
  configHash.insert("i_tempk_bcL",config.i_tempk_bcL);
  configHash.insert("i_tempk_bcR",config.i_tempk_bcR);
  configHash.insert("i_qv_bcL",config.i_qv_bcL);
  configHash.insert("i_qv_bcR",config.i_qv_bcR);
  configHash.insert("i_rhoa_bcL",config.i_rhoa_bcL);
  configHash.insert("i_rhoa_bcR",config.i_rhoa_bcR);
  configHash.insert("i_qr_bcL",config.i_qr_bcL);
  configHash.insert("i_qr_bcR",config.i_qr_bcR);
  configHash.insert("j_rhou_bcL",config.j_rhou_bcL);
  configHash.insert("j_rhou_bcR",config.j_rhou_bcR);
  configHash.insert("j_rhov_bcL",config.j_rhov_bcL);
  configHash.insert("j_rhov_bcR",config.j_rhov_bcR);
  configHash.insert("j_rhow_bcL",config.j_rhow_bcL);
  configHash.insert("j_rhow_bcR",config.j_rhow_bcR);
  configHash.insert("j_tempk_bcL",config.j_tempk_bcL);
  configHash.insert("j_tempk_bcR",config.j_tempk_bcR);
  configHash.insert("j_qv_bcL",config.j_qv_bcL);
  configHash.insert("j_qv_bcR",config.j_qv_bcR);
  configHash.insert("j_rhoa_bcL",config.j_rhoa_bcL);
  configHash.insert("j_rhoa_bcR",config.j_rhoa_bcR);
  configHash.insert("j_qr_bcL",config.j_qr_bcL);
  configHash.insert("j_qr_bcR",config.j_qr_bcR);
  configHash.insert("k_rhou_bcL",config.k_rhou_bcL);
  configHash.insert("k_rhou_bcR",config.k_rhou_bcR);
  configHash.insert("k_rhov_bcL",config.k_rhov_bcL);
  configHash.insert("k_rhov_bcR",config.k_rhov_bcR);
  configHash.insert("k_rhow_bcL",config.k_rhow_bcL);
  configHash.insert("k_rhow_bcR",config.k_rhow_bcR);
  configHash.insert("k_tempk_bcL",config.k_tempk_bcL);
  configHash.insert("k_tempk_bcR",config.k_tempk_bcR);
  configHash.insert("k_qv_bcL",config.k_qv_bcL);
  configHash.insert("k_qv_bcR",config.k_qv_bcR);
  configHash.insert("k_rhoa_bcL",config.k_rhoa_bcL);
  configHash.insert("k_rhoa_bcR",config.k_rhoa_bcR);
  configHash.insert("k_qr_bcL",config.k_qr_bcL);
  configHash.insert("k_qr_bcR",config.k_qr_bcR);
  configHash.insert("mode",config.mode);

  configHash.insert("qr_variable",config.qr_variable);
  configHash.insert("radar_dbz",config.radar_dbz);
  configHash.insert("radar_vel",config.radar_vel);
  configHash.insert("radar_sw",config.radar_sw);
  configHash.insert("mask_reflectivity",config.mask_reflectivity);
 */
#if 0 /* old style */
  configHash.insert("ref_time",config.ref_time);
#endif

/*fixme
  configHash.insert("ref_state",config.ref_state);
  configHash.insert("data_directory",config.data_directory);
  configHash.insert("output_directory",config.output_directory);

  configHash.insert("projection",config.projection);
*/
  // Validate the hash -- multiple passes are not validated currently
  std::vector<std::string> configKeys;
/*fixme
  configKeys << "ref_state" << //  << "ref_time" << // "reflat" << "reflon" are set by the VarDriver
    "qr_variable" << "i_background_roi" << "j_background_roi" <<
    "i_reflectivity_roi" << "j_reflectivity_roi" << "k_reflectivity_roi" <<
    "load_background" << "adjust_background" <<
    "radar_dbz" << "radar_vel" << "radar_sw" << "radar_skip" <<
    "radar_stride" << "dynamic_stride" << "dbz_pseudow_weight" <<
    "melting_zone_width" << "mixed_phase_dbz" << "rain_dbz" <<
    "num_iterations" << "output_mish" << "output_txt" << "output_qc" <<
    "output_netcdf" << "output_asi" << "output_COAMPS" << "preprocess_obs" << "mask_reflectivity" <<
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
    "bg_qv_error" << "bg_rhoa_error" << "bg_qr_error" << "projection";

  if (fixedGrid)
    configKeys << "ref_time";
*/

  for (int i = 0; i < configKeys.size(); i++) {
/*fixme    if (!configHash.contains(configKeys.at(i))) {
      std::cout << "No configuration found for <" << configKeys.at(i).toStdString() << "> aborting..." << std::endl;
      return false;
    }*/
  }
  return true;
}

Projection::ProjectionType VarDriver::projectionFromConfig()
{
  // default value
  Projection::ProjectionType retVal =
    Projection::TRANSVERSE_MERCATOR_EXACT;
  /*fixme
  if (configHash.contains("projection")) {
    if (configHash.value("projection") == "lambert_conformal_conic")
      retVal = Projection::LAMBERT_CONFORMAL_CONIC;
    else if (configHash.value("projection") != "transverse_mercator_exact")
      std::cerr << "Warning: Unrecognized projection type "
		<< configHash.value("projection").toLatin1().data()
		<< ". Defaulting to transverse_mercator_exact\n";
  }*/
  return retVal;
}

bool VarDriver::read_mesonet(std::string& filename, std::vector<MetObs>* metObVector)
{
  Nc3Error err(Nc3Error::verbose_nonfatal);

  // Read the file.
  Nc3File dataFile(filename.c_str(), Nc3File::ReadOnly);

  // Check to see if the file was read.
  if(!dataFile.is_valid())
    return false;

  // Get the number of records
  Nc3Dim* recnum;
  if (!(recnum = dataFile.get_dim("recNum")))
    return false;
  int NREC = recnum->size();

  Nc3Var *latVar, *lonVar, *altVar, *timeVar, *tempVar,
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

      datetime datetime_;
      datetime_ = std::chrono::seconds(long(obtime[rec])) + std::chrono::time_point<std::chrono::system_clock>{}; // bpd6 - casting issue?  What's 'rec'?

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

bool VarDriver::read_classnc(std::string& filename, std::vector<MetObs>* metObVector)
{
  Nc3Error err(Nc3Error::verbose_nonfatal);

  // Read the file.
  Nc3File dataFile(filename.c_str(), Nc3File::ReadOnly);

  // Check to see if the file was read.
  if(!dataFile.is_valid())
    return false;

  // Get the number of records
  // Nc3Dim* recnum;
  Nc3Dim* recnum;
  if (!(recnum = dataFile.get_dim("time")))
    return false;
  int NREC = recnum->size();

  Nc3Var *latVar, *lonVar, *altVar, *timeVar, *tempVar,
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

  datetime datetime_;
  std::vector<std::string> fileparts = LineSplit(filename, '.');

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

      auto combined = fileparts[1] + fileparts[2];
      datetime datetime_ = ParseDate(combined.c_str(), "%Y%m%d%H%M%S") + std::chrono::seconds(long(obtime[rec]));
      ob.setTime(datetime_);

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

bool VarDriver::read_qcf(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);

  if (!metFile.is_open())
    return false;

/*fixme
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
*/
  metFile.close();
  return true;

}

bool VarDriver::read_aeri(std::string& filename, std::vector<MetObs>* metObVector)
{
  Nc3Error err(Nc3Error::verbose_nonfatal);

  // Read the file.
  Nc3File dataFile(filename.c_str());

  // Check to see if the file was read.
  if(!dataFile.is_valid())
    return false;

  // Get the number of records
  Nc3Dim* recnum;
  if (!(recnum = dataFile.get_dim("time")))
    return false;
  int NREC = recnum->size();

  // Get the number of pressure levels
  Nc3Dim* pres_level;
  if (!(pres_level = dataFile.get_dim("pres_level")))
    return false;
  int NLVL = pres_level->size();

  Nc3Var *latVar, *lonVar, *altVarBase, *altVar, *timeVarBase, *timeVar,
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

  datetime datetime_;

  // Get the platform name
  MetObs ob;
  std::vector<std::string> fileparts = LineSplit(filename,'_');
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

/*fixme
      datetime = QDateTime::fromTime_t(basetime);
      QString obstring = datetime.toString(Qt::ISODate);
      QDateTime obdatetime = datetime.addSecs(obtime[rec]).toUTC();
      obstring = obdatetime.toString(Qt::ISODate);
      ob.setTime(obdatetime);
*/

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

/* This routine reads the ascii radar data format*/

bool VarDriver::read_rad(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);

  if (!metFile.is_open()) {
    return false;
	}

  // Use a Transverse Mercator projection to map the radar gates to the grid
  // GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();
/*fixme
  QTextStream in(&metFile);
  MetObs ob;
  QDate startDate;
  QDateTime datetime;
  QString line;
  while (!in.atEnd()) {
    line = in.readLine();
    QStringList lineparts = line.split(QRegExp(","));
    QDate date = QDate::fromString(lineparts[0].mid(0,10), "yyyy-MM-dd");
    QTime time = QTime::fromString(lineparts[0].mid(11,8), "HH:mm:ss");
    datetime = QDateTime(date, time, Qt::UTC);
    ob.setTime(datetime);
    real radarLat = lineparts[1].toFloat();
    real radarLon = lineparts[2].toFloat();
    real radarAlt = lineparts[3].toFloat();
    real az = lineparts[4].toFloat();
    real el = lineparts[5].toFloat();
    real range = lineparts[6].toFloat();
    real relX = range*sin(az*Pi/180.)*cos(el*Pi/180.);
    real relY = range*cos(az*Pi/180.)*cos(el*Pi/180.);
    real rEarth = 6371000.0;

    // Take into account curvature of the earth for the height of the radar beam
    real relZ = sqrt(range*range + rEarth*rEarth + 2.0 * range * rEarth * sin(el*Pi/180.)) - rEarth;
    real radarX, radarY, gateLat, gateLon;
    projection.Forward(radarLon, radarLat, radarLon, radarX, radarY);
    projection.Reverse(radarLon, radarX + relX, radarY + relY, gateLat, gateLon);
    real gateAlt = relZ + radarAlt;

    ob.setLat(gateLat);
    ob.setLon(gateLon);
    ob.setAltitude(gateAlt);
    ob.setAzimuth(az);
    ob.setElevation(el);
    if(lineparts[7].toFloat() != -32768.0) {
      ob.setReflectivity(lineparts[7].toFloat());
    } else {
      ob.setReflectivity(-999.0);
    }
    if(lineparts[8].toFloat() != -32768.0) {
      ob.setRadialVelocity(lineparts[8].toFloat());
    } else {
      ob.setRadialVelocity(-999.0);
    }
    if(lineparts[9].toFloat() != -32768.0) {
      ob.setSpectrumWidth(lineparts[9].toFloat());
    } else {
      ob.setSpectrumWidth(-999.0);
    }
    ob.setObType(MetObs::radar);
    metObVector->push_back(ob);
  }
*/

  metFile.close();
  return true;
}

// This routing reads the Lrose Radx format

bool VarDriver::read_cfrad(std::string &fileName, std::vector<MetObs>* metObVector)
{
  RadxFile rxFile;
  RadxVol rxVol;

  if ( ! rxFile.isSupported(fileName) ) {
    std::cerr << "ERROR - File '" << fileName
             << "' is not in a Radx supported format." << std::endl;
    return false;
  }

  rxFile.setReadPreserveSweeps(true); // prevent Radx from tossing away long sweeps

  if (rxFile.readFromPath(fileName, rxVol) ) {
    std::cerr << "ERROR - reading file: " << fileName << std::endl;
    std::cerr << rxFile.getErrStr() << std::endl;
    return false;
  }

  // Name of the fields to use
  // TODO Should we use cfradial defaults if not specified?
  // DBZ, VEL, and WIDTH?

  std::string radarDbzStr = configHash["radar_dbz"];
  std::string radarVelStr = configHash["radar_vel"];
  std::string radarSwStr  = configHash["radar_sw"];

  // Use a Transverse Mercator projection to map the radar gates to the grid
  // GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();

#if 0
  double tmStart = rxVol.getStartTimeSecs() + 1.e-9 * rxVol.getStartNanoSecs();
  double tmEnd = rxVol.getEndTimeSecs() + 1.e-9 * rxVol.getEndNanoSecs();
#endif

  vector<RadxRay *> rays = rxVol.getRays();

  for (size_t index = 0; index < rays.size(); index++) {
    RadxRay *ray = rays[index];
    if (ray == NULL) {
      std::cout << "Error: Failed to read ray " << index
               << " from file " << fileName << std::endl;
      return false;
    }

    double radarLat = numeric_limits<double>::quiet_NaN();     // degrees
    double radarLon = numeric_limits<double>::quiet_NaN();     // degrees
    double radarAlt = numeric_limits<double>::quiet_NaN();     // km

    double beamWidth = sin(ray->getFixedAngleDeg() * Pi / 180.0);

    const RadxGeoref *gref = ray->getGeoreference();
    if (gref == NULL ) { // gref is NULL for some ground-based stations
      radarLat = rxVol.getLatitudeDeg();
      radarLon = rxVol.getLongitudeDeg();
      radarAlt = rxVol.getAltitudeKm();
    } else {
      radarLat = gref->getLatitude();
      radarLon = gref->getLongitude();
      radarAlt = gref->getAltitudeKmMsl();
    }
    if (std::isnan(radarLat)
       || std::isnan(radarLon)
       || std::isnan(radarAlt))
    {
      std::cout << "Error: incomplete file spec or radx Georeference" << std::endl;
      return false;
    }

    // Qt 5.8 //  QDateTime rayTime = QDateTime::fromSecsSinceEpoch(ray->getTimeSecs());
    datetime rayTime = date::sys_seconds(std::chrono::seconds(ray->getTimeSecs()));
    double az = ray->getAzimuthDeg();
    double el = ray->getElevationDeg();

    // Get the ref, vel, and swdata

    RadxField *radarDbz = ray->getField(radarDbzStr);
    if (radarDbz == NULL) {
      std::cout << "Failed to get variable " << radarDbzStr << " from " << fileName << std::endl;
      return false;
    }

    RadxField *radarVel = ray->getField(radarVelStr);
    if (radarVel == NULL) {
      std::cout << "Failed to get variable " << radarVelStr << " from " << fileName << std::endl;
      return false;
    }
    RadxField *radarSw  = ray->getField(radarSwStr);
    if (radarSw == NULL) {
      std::cout << "Failed to get variable " << radarSwStr << " from " << fileName << std::endl;
      return false;
    }

    double dbzMissingVal = radarDbz->getMissingFl64();
    double velMissingVal = radarVel->getMissingFl64();
    double swMissingVal =  radarSw->getMissingFl64();

    size_t nGates = ray->getNGates();
    float gatelength = ray->getGateSpacingKm() * 1000;

    //    int rayskip = configHash.value("radar_skip").toInt();
    int minstride = std::stoi(configHash["radar_stride"]);
    bool dynamicStride = std::stoi(configHash["dynamic_stride"]);
    int stride = minstride;

    // std::cout << "-I- Gates: " << nGates << ", stride: " << stride
    // << ", length: " << gatelength << std::endl;

    for (size_t gateIndex = 0; gateIndex < nGates - stride; gateIndex += stride) {
      MetObs ob;
      float range = gatelength * (gateIndex + stride / 2);
      if (dynamicStride) {
	stride = (int) (range * beamWidth / gatelength);
	if (stride < minstride) stride = minstride;
      }
      real dz = 0.0;
      real vr = 0.0;
      real sw = 0.0;
      int dzCount = 0;
      int vrCount = 0;
      int swCount = 0;

      for (size_t idx = gateIndex; idx < (gateIndex + stride); idx++) {

	double valDbz = radarDbz->getDoubleValue(idx);
	double valVel = radarVel->getDoubleValue(idx);
	double valSw  = radarSw->getDoubleValue(idx);

	if (valVel != velMissingVal) {
	  vr += valVel;
	  vrCount++;
	}
	if (valDbz != dbzMissingVal) {
	  dz += pow(10.0, valDbz * 0.1);
	  dzCount++;
	}
	if(valSw != swMissingVal) {
	  sw += valSw;
	  swCount++;
	}
      }

      if (dzCount > 0) {
	dz = dz / float(dzCount);
	dz = 10 * log10(dz);
      } else
	dz = -999.0;
      if (vrCount > 0) {
	vr = vr / float(vrCount);
      } else {
	vr = -999.0;
      }
      if (swCount > 0) {
	sw = sw / float(swCount);
      } else {
	sw = -999.0;
      }

      if ((vr != -999.0) || (dz != -999.0)) {
	real relX = range * sin(az * Pi / 180.0) * cos(el * Pi / 180.0);
	real relY = range * cos(az * Pi / 180.0) * cos(el * Pi / 180.0);
	real rEarth = 6371000.0;

	// Take into account curvature of the earth for the height of the radar beam
	real relZ = sqrt(range * range + rEarth * rEarth + 2.0 * range
			 * rEarth * sin(el * Pi / 180.0)) - rEarth;

	real radarX, radarY, gateLat, gateLon;
	projection.Forward(radarLon, radarLat, radarLon, radarX, radarY);
	projection.Reverse(radarLon, radarX + relX, radarY + relY, gateLat, gateLon);
	real gateAlt = relZ + radarAlt * 1000;

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
    } // gates
  } // rays
  return true;
}

/* This routine reads a text file from a crsim radar-filter output*/

bool VarDriver::read_crsim(std::string& filename, std::vector<MetObs>* metObVector)
{
	std::ifstream metFile(filename);
  if (!metFile.is_open()) {
    return false;
	}
	MetObs ob;
	std::string line;
	std::getline(metFile, line);
	auto parts = LineSplit(line, ' ');
	datetime datetime_ = ParseTime(parts[0].c_str(), "%Y-%m-%d_%H:%M:%S");
	real radarLat = std::stof(parts[1]);
	real radarLon = std::stof(parts[2]);
	real radarAlt = std::stof(parts[3]);

  while (std::getline(metFile, line)) {
		auto parts = LineSplit(line, ' ');
		real az = std::stof(parts[0]);
		real el = std::stof(parts[1]);
		real range = std::stof(parts[2]);
		real dz = std::stof(parts[3]);
		real vr = std::stof(parts[4]);
		if ((vr != -999.0) || (dz != -999.0)) {
		real relX = range * sin(az * Pi / 180.0) * cos(el * Pi / 180.0);
		real relY = range * cos(az * Pi / 180.0) * cos(el * Pi / 180.0);
		real rEarth = 6371000.0;

		// Take into account curvature of the earth for the height of the radar beam
		real relZ = sqrt(range * range + rEarth * rEarth + 2.0 * range
				 * rEarth * sin(el * Pi / 180.0)) - rEarth;

		real radarX, radarY, gateLat, gateLon;
		projection.Forward(radarLon, radarLat, radarLon, radarX, radarY);
		projection.Reverse(radarLon, radarX + relX, radarY + relY, gateLat, gateLon);
		real gateAlt = relZ + radarAlt * 1000;

		ob.setTime(datetime_);
		ob.setLat(gateLat);
		ob.setLon(gateLon);
		ob.setAltitude(gateAlt);
		ob.setAzimuth(az);
		ob.setElevation(el);
		ob.setRadialVelocity(vr);
		ob.setReflectivity(dz);
    ob.setObType(MetObs::crsim);
    metObVector->push_back(ob);
  }}
	std::cout << "Successfully read the crsim file" << std::endl;
  metFile.close();
  return true;
}

bool VarDriver::read_hrdradial(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);
  if (!metFile.is_open()) {
    return false;
    }
    MetObs ob;
    std::string line;

    // These lines are not needed data, but must be read through
    std::getline(metFile, line); // Date
    std::getline(metFile, line); // Aircraft number NXXRF (42, 43, or 49)
    std::getline(metFile, line); // Storm number
    
    std::getline(metFile, line); // How many radial bins
    int num_gates = std::stoi(line);
    
    std::getline(metFile, line); // The radius of the first bin
    real radius_first_bin = std::stof(line);

    std::getline(metFile, line); // The radial resolution of the data
    real range_bin_resolution = std::stof(line);
    
    std::getline(metFile, line); // Octal for type of data.
    std::getline(metFile, line); // Type of antenna.
    
    /* Epoch time, aircraft latitude, aircraft longitude,
     aircraft altitude (m), mathematical earth-relative azimuth of
     pointing angle (0 east, 90 north, 180 west, 270 south),
     earth-relative elevation angle */
    std::string radial;
    while (std::getline(metFile, radial)) {
        auto header = LineSplit(radial, ' ');
        std::time_t tt = std::stof(header[0]);
        std::gmtime(&tt);
        datetime datetime_ = std::chrono::system_clock::from_time_t(tt);
        real radarLat = std::stof(header[1]);
        real radarLon = std::stof(header[2]);
        real radarAlt = std::stof(header[3]);
        real az = 90.0 - std::stof(header[4]);
        real el = std::stof(header[5]);
        
        int num_lines = num_gates / 10;
        if ( num_gates % 10 ) {
            num_lines += 1;
        }
        real range = radius_first_bin;
        for (size_t idx = 0; idx < num_lines; idx++) {
            std::getline(metFile, line);
            auto parts = LineSplit(line, ' ');
            std::vector<std::string>::iterator ptr;
            for (ptr = parts.begin(); ptr < parts.end(); ptr++) {
                real vr = std::stof(*ptr);
                range += range_bin_resolution;
                if (vr != -8888.0) {
                    real vr_scale = vr / 10.0;
                    real relX = range * sin(az * Pi / 180.0) * cos(el * Pi / 180.0);
                    real relY = range * cos(az * Pi / 180.0) * cos(el * Pi / 180.0);
                    real rEarth = 6371000.0;
                    
                    // Take into account curvature of the earth for the height of the radar beam
                    real relZ = sqrt(range * range + rEarth * rEarth + 2.0 * range
                                     * rEarth * sin(el * Pi / 180.0)) - rEarth;
                    
                    real radarX, radarY, gateLat, gateLon;
                    projection.Forward(radarLon, radarLat, radarLon, radarX, radarY);
                    projection.Reverse(radarLon, radarX + relX, radarY + relY, gateLat, gateLon);
                    real gateAlt = relZ + radarAlt;
                    
                    ob.setTime(datetime_);
                    ob.setLat(gateLat);
                    ob.setLon(gateLon);
                    ob.setAltitude(gateAlt);
                    ob.setAzimuth(az);
                    ob.setElevation(el);
                    ob.setRadialVelocity(vr_scale);
                    // Set a fake reflectivity
                    ob.setReflectivity(10.0);
                    ob.setObType(MetObs::radar);
                    metObVector->push_back(ob);
                }
            }
        }
    }
    std::cout << "Successfully read the hrd radial file" << std::endl;
    metFile.close();
    return true;
}

bool VarDriver::read_hdobs(std::string& filename, std::vector<MetObs>* metObVector)
{
  std::ifstream metFile(filename);
  if (!metFile.is_open()) {
    return false;
    }
    MetObs ob;
    std::string line;

    // These lines are not needed data, but must be read through
    std::getline(metFile, line); // Blank
    std::getline(metFile, line); // Not sure what this is
    std::getline(metFile, line); // HDOB header
    
    std::getline(metFile, line); // Storm, aircraft, and date info
    auto header = LineSplit(line, ' ');
    string date = header[5];
    int year = std::stoi(date.substr(0,4));
    int mon = std::stoi(date.substr(4,2));
    int day = std::stoi(date.substr(6,2));
    
    std::string hdob;
    while (std::getline(metFile, hdob)) {
        if (hdob.empty()) {
            break;
        }
        if (hdob.compare("$$") == 0) {
            break;
        }
        auto obdata = LineSplit(hdob, ' ');
        string obtime = obdata[0];
        int hr = std::stoi(obtime.substr(0,2));
        int min = std::stoi(obtime.substr(2,2));
        int sec = std::stoi(obtime.substr(4,2));
        
        std::time_t tt;
        time ( &tt );
        struct tm *timeinfo;
        timeinfo = gmtime ( &tt );
        timeinfo->tm_year = year - 1900;
        timeinfo->tm_mon = mon - 1;
        timeinfo->tm_mday = day;
        timeinfo->tm_hour = hr;
        timeinfo->tm_min = min;
        timeinfo->tm_sec = sec;
        tt = timegm ( timeinfo );
        datetime datetime_ = std::chrono::system_clock::from_time_t(tt);
        
        real obLatDeg = std::stof(obdata[1].substr(0,2));
        real obLatMin = std::stof(obdata[1].substr(2,2));
        string obLatHem = obdata[1].substr(4,1);
        real obLat = obLatDeg + obLatMin/60.0;
        if (obLatHem.compare("S") == 0) {
            obLat = -obLat;
        }
        
        real obLonDeg = std::stof(obdata[2].substr(0,3));
        real obLonMin = std::stof(obdata[2].substr(3,2));
        string obLonHem = obdata[2].substr(5,1);
        real obLon = obLonDeg + obLonMin/60.0;
        if (obLonHem.compare("W") == 0) {
            obLon = -obLon;
        }
        
        real obPressure = std::stof(obdata[3]);
        // Check to see if pressure is above 1000 hPa
        if (obPressure < 1000.0) {
            obPressure = obPressure + 10000.0;
        }
        obPressure = obPressure / 10.0;
        
        real obAlt = std::stof(obdata[4]);
        
        string press = obdata[5];
        real obSfcpress = -999.0;
        if (press.compare("////") != 0) {
            obSfcpress = std::stof(press);
            // Check to see if pressure is above 1000 hPa
            if (obSfcpress < 1000.0) {
                obSfcpress = obSfcpress + 10000.0;
            }
            obSfcpress = obSfcpress / 10.0;
        }
        
        string temp = obdata[6].substr(1,3);
        real obTemp = -999.0;
        if (temp.compare("///") != 0) {
            obTemp = std::stof(temp);
        }
        if (obdata[6].substr(1,3).compare("-") == 0) {
            obTemp = -obTemp;
        }
        obTemp = (obTemp/10.0) + 273.15;
        
        string dewp = obdata[7].substr(1,3);
        real obDewp = -999.0;
        if (dewp.compare("///") != 0) {
            obDewp = std::stof(dewp);
        }
        if (obdata[7].substr(1,3).compare("-") == 0) {
            obDewp = -obDewp;
        }
        obDewp = (obDewp/10.0) + 273.15;
        
        real obWdir = std::stof(obdata[8].substr(0,3));
        //real obWspd30s = std::stof(obdata[8].substr(3,3)) * 0.514;
        string wspd = obdata[9];
        real obWspd10s = -999.0;
        if (wspd.compare("///") != 0) {
            obWspd10s = std::stof(wspd) * 0.514;
        }
        wspd = obdata[10];
        real obSFMR10s = -999.0;
        if (wspd.compare("///") != 0) {
            obSFMR10s = std::stof(wspd) * 0.514;
        }
        
        // Ignoring SFMR rain rate for now, but could be useful
        
        int qualityflag = std::stoi(obdata[12]);
        
        ob.setTime(datetime_);
        ob.setLat(obLat);
        ob.setLon(obLon);
        ob.setAltitude(obAlt);
        ob.setPressure(obPressure);
        ob.setTemperature(obTemp);
        ob.setDewpoint(obDewp);
        ob.setWindSpeed(obWspd10s);
        ob.setWindDirection(obWdir);
        ob.setObType(MetObs::flightlevel);
        metObVector->push_back(ob);
        
        /*  Pseudo-surface ob
        ob.setAltitude(10.0);
        ob.setPressure(obSfcpress);
        // Crudely estimate the surface temperature with a 8 K/km lapse rate
        real sfcTemp = obTemp + (8.0 * obAlt/1000.0);
        ob.setTemperature(sfcTemp);
        ob.setDewpoint(sfcTemp - 1.0);
        ob.setWindSpeed(obSFMR10s);
        ob.setWindDirection(obWdir);
        ob.setObType(MetObs::sfmr);
        //metObVector->push_back(ob); */

    }
        
    std::cout << "Successfully read the hdob file" << std::endl;
    metFile.close();
    return true;
}
