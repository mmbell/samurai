/*
 *  MetObs.cpp
 *  tcvar
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "MetObs.h"
#include <QTextStream>
#include "math.h"
#include <iostream>

MetObs::MetObs()
{
	latitude = -999;
	longitude = -999;
	altitude = -999;
	pressure = -999;
	windSpeed = -999;
	windDirection = -999;
	verticalVelocity = -999;
	dewpoint = -999;
	temperature = -999;
	radialVelocity = -999;
	reflectivity = -999;
	spectrumWidth = -999;
	azimuth = -999;
	elevation = -999;
	stationName = QString();
	time = QDateTime();	
	obType = -1;

}

MetObs::MetObs(const MetObs& other)
{
	latitude = other.latitude;
	longitude = other.longitude;
	altitude = other.altitude;
	pressure = other.pressure;
	windSpeed = other.windSpeed;
	windDirection = other.windDirection;
	verticalVelocity = other.verticalVelocity;
	dewpoint = other.dewpoint;
	temperature = other.temperature;
	radialVelocity = other.radialVelocity;
	reflectivity = other.reflectivity;
	spectrumWidth = other.spectrumWidth;
	azimuth = other.azimuth;
	elevation = other.elevation;
	stationName = other.stationName;
	time = other.time;
	obType = other.obType;
	
}

MetObs::~MetObs()
{
}

bool MetObs::readObs()
{
	// Virtual function
	return false;
}

QString MetObs::getStationName() const
{
	return stationName;
}

void MetObs::setStationName(const QString& name)
{
	stationName = name;
}

float MetObs::getLat() const
{
	return latitude;
}

void MetObs::setLat(const float& lat)
{
	latitude = lat;
}


float MetObs::getLon() const
{
	return longitude;
}

void MetObs::setLon(const float& lon)
{
	longitude = lon;
}

float MetObs::getAltitude() const
{
	return altitude;
}

void MetObs::setAltitude(const float& alt)
{
	altitude = alt;
}

QDateTime MetObs::getTime() const
{
	return time;
}

void MetObs::setTime(const QDateTime& obTime)
{
	time = QDateTime(obTime);
}

float  MetObs::getPressure() const
{
	return pressure;
}

void  MetObs::setPressure(const float& press)
{
	pressure = press;
}

float  MetObs::getWindSpeed() const
{
	return windSpeed;
}

void  MetObs::setWindSpeed(const float& speed)
{
	windSpeed = speed;
}

float  MetObs::getWindDirection() const 
{
	return windDirection;
}

void  MetObs::setWindDirection(const float& dir)
{
	windDirection = dir;
}


float MetObs::getVerticalVelocity() const
{
	return verticalVelocity;
}

void MetObs::setVerticalVelocity(const float& w)
{
	verticalVelocity = w;
}

float MetObs::getTemperature() const
{
	return temperature;
}

void MetObs::setTemperature(const float& T)
{
	temperature = T;
}

float MetObs::getDewpoint() const
{
	return dewpoint;
}

void MetObs::setDewpoint(const float& D)
{
	dewpoint = D;
}

void MetObs::setRH(const float& RH)
{
	if ((RH > 0) and (RH < 101) and (temperature != -999)) { 
		float t = temperature;
		float es = (E_3 * exp (_A_ * log (T_3 / t)) * 
					exp ((_A_ + _B_) * (1 - T_3 / t)));
		float e = es*RH/100;
		float u = log (e / E_3);	
		dewpoint = (237.3 * u / (17.2694 - u) + T_3);
	} else {
		dewpoint = -999;
	}
}

float MetObs::getRadialVelocity() const
{
	return radialVelocity;
}

void MetObs::setRadialVelocity(const float& vr)
{
	radialVelocity = vr;
}

float MetObs::getReflectivity() const
{
	return reflectivity;
}

void MetObs::setReflectivity(const float& dz)
{
	reflectivity = dz;
}

float MetObs::getSpectrumWidth() const
{
	return spectrumWidth;
}

void MetObs::setSpectrumWidth(const float& sw)
{
	spectrumWidth = sw;
}

float MetObs::getAzimuth() const
{
	return azimuth;
}

void MetObs::setAzimuth(const float& az)
{
	azimuth = az;
}

float MetObs::getElevation() const
{
	return elevation;
}

void MetObs::setElevation(const float& el)
{
	elevation = el;
}

int MetObs::getObType() const
{
	return obType;
}

void MetObs::setObType(const int& type)
{
	obType = type;
}

// Derived variables
float MetObs::getQv() const
{	
	if ((dewpoint != -999) and (pressure != -999)) {
		float e = (E_3 * exp (17.2694 * (1.0 - 237.3 / (dewpoint - T_3 + 237.3))));
		return (1000.0 * EPSILON * e / (pressure - e));
	} else {
		return -999;
	}
}

float MetObs::getQvSaturation() const
{	
	if ((temperature != -999) and (pressure != -999)) {
		float e = (E_3 * exp (17.2694 * (1.0 - 237.3 / (temperature - T_3 + 237.3))));
		return (1000.0 * EPSILON * e / (pressure - e));
	} else {
		return -999;
	}
}


float  MetObs::getCartesianUwind() const
{
	if (windSpeed != -999) {
		float uWind = -windSpeed * sin(windDirection * acos(-1) / 180.);
		return uWind;
	} else {
		return -999;
	}
}

float  MetObs::getCartesianVwind() const
{
	if (windSpeed != -999) {
		float vWind = -windSpeed * cos(windDirection * acos(-1) / 180.);
		return vWind;
	} else {
		return -999;
	}
}

float MetObs::getPolarUwind(const float& centerLat, const float& centerLon) const
{
	
	float centerLatRadians = centerLat * acos(-1.0)/180.0;
	float fac_lat = 111.13209 - 0.56605 * cos(2.0 * centerLatRadians)
	+ 0.00012 * cos(4.0 * centerLatRadians) - 0.000002 * cos(6.0 * centerLatRadians);
	float fac_lon = 111.41513 * cos(centerLatRadians)
	- 0.09455 * cos(3.0 * centerLatRadians) + 0.00012 * cos(5.0 * centerLatRadians);
	float relX = (longitude - centerLon) * fac_lon;
	float relY = (latitude - centerLat) * fac_lat;
	float dist = sqrt(relX*relX + relY*relY);
	float cs = relY/dist;
	float sn = -relX/dist;
	float vr = -(getCartesianUwind())*sn + getCartesianVwind()*cs;
	return vr;
}

float MetObs::getPolarVwind(const float& centerLat, const float& centerLon) const
{
	float centerLatRadians = centerLat * acos(-1.0)/180.0;
	float fac_lat = 111.13209 - 0.56605 * cos(2.0 * centerLatRadians)
	+ 0.00012 * cos(4.0 * centerLatRadians) - 0.000002 * cos(6.0 * centerLatRadians);
	float fac_lon = 111.41513 * cos(centerLatRadians)
	- 0.09455 * cos(3.0 * centerLatRadians) + 0.00012 * cos(5.0 * centerLatRadians);
	float relX = (longitude - centerLon) * fac_lon;
	float relY = (latitude - centerLat) * fac_lat;
	float dist = sqrt(relX*relX + relY*relY);
	float cs = relY/dist;
	float sn = -relX/dist;
	float vt = -(getCartesianUwind())*cs - getCartesianVwind()*sn;
	return vt;
}

float MetObs::getDryDensity() const
{
	if ((temperature != -999) and (pressure != -999)) {
		return (pressure*100)/(R_D*temperature);
	} else {
		return -999;
	}
}

float MetObs::getMoistDensity() const
{
	if ((dewpoint != -999) and (pressure != -999)) {
		return (pressure*100)/(R_D*getVirtualTemp());
	} else {
		return -999;
	}
}

float MetObs::getVirtualTemp() const
{
	float t = temperature;
	float e = (E_3 * exp (_A_ * log (T_3 / t)) * 
			   exp ((_A_ + _B_) * (1 - T_3 / t)));	
	return (t / (1 - (e / pressure) * (1 - EPSILON)));
}

float MetObs::getMoistStaticEnergy() const
{
	float h = C_P * temperature + 2.5e3 * getQv() + 9.81 * altitude;
	return h;
}

float MetObs::getMoistSaturationStaticEnergy() const
{
	if ((temperature != -999) and (pressure != -999)) {
		float h = C_P * temperature + 2.5e3 * getQvSaturation() + 9.81 * altitude;
		return h;
	} else {
		return -999;
	}
}

bool MetObs::operator ==(const MetObs &other)
{
    
	if((this->time.time() == other.time.time())
	   &&(this->stationName==other.stationName))
		return true;
	return false;
}

bool MetObs::operator < (const MetObs &other)
{
	if(this->time < other.time)
		return true;
	return false;
}

bool MetObs::operator > (const MetObs &other)
{
	if(this->time > other.time)
		return true;
	return false;
}

void MetObs::printString()
{
	
	QString printMessage(getStationName()+"_"+getTime().toString()+"_"+QString().setNum(getPressure()));
	//std::cout << printMessage.toAscii();
	
}

