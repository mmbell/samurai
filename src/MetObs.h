/*
 *  MetObs.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef METOBS_H
#define METOBS_H

#include "datetime.h"

/*
 * Constants
 */
# define _A_	5.0065
# define _B_	19.83923
# define C_P	1005.7
# define R_D	287.
# define E_3	6.1078
# define T_3	273.15
# define EPSILON	0.622

class MetObs
{

 public:
  MetObs();
  MetObs(const MetObs& other);
  virtual ~MetObs();
  virtual bool readObs();

  void printString();

  std::string getStationName() const;
  void setStationName(const std::string& name);

  float getLat() const;
  void setLat(const float& lat);

  float getLon() const;
  void setLon(const float& lon);

  float getAltitude() const;
  void setAltitude(const float& alt);

  datetime getTime() const;
  void setTime(const datetime& obTime);

  float getPressure() const;
  void setPressure(const float& press);

  float getWindSpeed() const;
  void setWindSpeed(const float& speed);
  float getWindDirection() const;
  void setWindDirection(const float& dir);

  float getVerticalVelocity() const;
  void setVerticalVelocity(const float& w);

  float getMeridionalVelocity() const;
  void setMeridionalVelocity(const float& v);

  float getZonalVelocity() const;
  void setZonalVelocity(const float& u);

  float getTemperature() const;
  void setTemperature(const float& T);
  float getTemperatureError() const;
  void setTemperatureError(const float& T);

  float getDewpoint() const;
  void setDewpoint(const float& D);
  void setRH(const float& RH);

  float getRadialVelocity() const;
  void setRadialVelocity(const float& vr);

  float getReflectivity() const;
  void setReflectivity(const float& dz);

  float getSpectrumWidth() const;
  void setSpectrumWidth(const float& sw);

  float getAzimuth() const;
  void setAzimuth(const float& az);

  float getElevation() const;
  void setElevation(const float& el);

  float getTerrainDX() const;
  void setTerrainDX(const float& dhdx);

  float getTerrainDY() const;
  void setTerrainDY(const float& dhdy);

  float getTerrainX() const;
  void setTerrainX(const float& x);

  float getTerrainY() const;
  void setTerrainY(const float& y);

  float getModelQv() const;
  void setModelQv(const float& qv);

  float getModelMoistDensity() const;
  void setModelMoistDensity(const float& rho);

  float getModelAirDensity() const;
  void setModelAirDensity(const float& rhoa);

  int getObType() const;
  void setObType(const int& type);

  // Derived variables
  float getQv() const;
  float getQvSaturation() const;
  float getCartesianUwind() const;
  float getCartesianVwind() const;
  float getPolarUwind(const float& centerLat, const float& centerLon) const;
  float getPolarVwind(const float& centerLat, const float& centerLon) const;
  float getDryDensity() const;
  float getAirDensity() const;
  float getVaporDensity() const;
  float getMoistDensity() const;
  float getVaporPressure() const;
  float getSatVaporPressure() const;
  float getVirtualTemp() const;
  float getMoistStaticEnergy() const;
  float getMoistSaturationStaticEnergy() const;
  float getTotalEnergy() const;

  bool operator ==(const MetObs &other);
  bool operator < (const MetObs &other);
  bool operator > (const MetObs &other);

  enum MetObTypes {
    dropsonde,
    flightlevel,
    radar,
    sfmr,
    qscat,
    ascat,
    AMV,
    lidar,
    insitu,
    mtp,
    mesonet,
    aeri,
    terrain,
    model,
    crsim,
    hrdradial
    // add terrain type
  };

 protected:

  float latitude;
  float longitude;
  float altitude;
  float pressure;
  datetime time;
  std::string stationName;
  float windSpeed;
  float windDirection;
  float verticalVelocity;
  float meridionalVelocity;
  float zonalVelocity;
  float dewpoint;
  float temperature;
  float temperatureError;
  float radialVelocity;
  float reflectivity;
  float spectrumWidth;
  float azimuth;
  float elevation;
  int obType;
  float terrain_dx;
  float terrain_dy;
  float terrain_x;
  float terrain_y;
  float moistDensity;
  float airDensity;
  float mixingRatio;
  // float differentialReflectivity; add Zdr
  // float
  // add terrain slope here


};

#endif
