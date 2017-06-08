#include <iostream>
#include "BkgdAdapter.h"

BkgdArray::BkgdArray(int nx, int ny, int nsigma,
		     float dx, float dy,
		     float *sigmas,
		     float *latitude,	// these 2 are 2D
		     float *longitude,
		     float *u1,		// these are 3D
		     float *v1,
		     float *w1,
		     float *th1,
		     float *p1)
{
  _x_incr = dx;
  _y_incr = dy;

  _xd = nx;
  _yd = ny;
  _zd = nsigma;

  _sigmas = sigmas;
  _lat_2d  = latitude;
  _long_2d = longitude;

  _u1_3d = u1;
  _v1_3d = v1;
  _w1_3d = w1;
  _th1_3d = th1;
  _p1_3d = p1;

  _curr_x = 0;
  _curr_y = 0;
  _curr_z = 0;
}

// Return the item at the given position
// a3d is a 3d array.

float BkgdArray::item3d(float *a3d, int x, int y, int z)
{
  return *(a3d + x + _xd * (y + _yd * z));
}

float BkgdArray::item2d(float *a2d, int x, int y)
{
  return *(a2d + x + _xd * y);
}

BkgdArray::~BkgdArray()
{
}

bool BkgdArray::next(int &time, real &lat, real &lon, real &alt, real &u,
	    real &v, real &w, real &t, real &qv, real &rhoa, real &qr)
{
  if(_curr_x >= _xd) return false;

  // Don't care about time
  
  lat = item2d(_lat_2d, _curr_x, _curr_y);
  lon = item2d(_long_2d, _curr_x, _curr_y);
    
  u = item3d(_u1_3d, _curr_x, _curr_y, _curr_z);
  v = item3d(_v1_3d, _curr_x, _curr_y, _curr_z);
  w = item3d(_w1_3d, _curr_x, _curr_y, _curr_z);

  // TODO
  // th1 is background pertubation potential temp on current nest
  // p1  is background pertubation exner velocity on current nest
  //
  t = item3d(_th1_3d, _curr_x, _curr_y, _curr_z);
  qv = item3d(_p1_3d, _curr_x, _curr_y, _curr_z);
  // Need a map to do a linear transform to go from these to samuari_background.in format

  // Need to get moisture: humidity(nx, ny, nsigma)

  // qv = ?
  // rhoa = ?
  qr = 0; // for now

  std::cout << "(" << _curr_x << ", " << _curr_y << ", " << _curr_z;
  std::cout << ")\t";
  std::cout << "lat: " << lat << ", long: " << lon << ", alt: " << alt;
  std::cout << " u: " << u << ", v: " << v << ", w: " << w;
  std::cout << " t: " << t << ", qv: " << qv << ", rhoa: " << rhoa << ", qr: " << qr << std::endl;
  // TODO Double (and triple) check this.
  // From the Samurai doc:
  // The data does not have to be evenly spaced, but the ordering of the positions is important.
  // Data should be sorted by altitude, such that each successive column of data is grouped together
  // with altitude increasing.

  _curr_z += 1;
  if(_curr_z >= _zd) {
    _curr_z %= _zd;
    _curr_y += 1;
    if(_curr_y >= _yd) {
      _curr_y %= _yd;
      _curr_x += 1;
    }
  }
  
  return true;
}
