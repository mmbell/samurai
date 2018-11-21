#ifndef BKGD_ADAPTER_H
#define BKGD_ADAPTER_H

#include <fstream>
#include <QDateTime>

#include "precision.h"	// for 'real' typedef

// This is an abstract class to support the multiple ways Background
// observations can be loaded into samurai.
//
// Right now we support loading from
//   a file (samurai_Background.in) The format is documented on the
//       Samurai github wiki
//   a set of arrays. This is a new format to support COAMPS.
//       This interface is a work in progress.
//   Both the methods above expect to iterate over each iteration (line per line, or array rows)
//   a Fractl generated netcdf file
//       This one doesn't need to act as an iterator since no interpolation is done.
//       So it really is just a placeholder

class BkgdAdapter {

 public:
  BkgdAdapter() {};
  virtual ~BkgdAdapter() {};
  
  virtual bool next(int &time, real &lat, real &lon, real &alt, real &u,
		    real &v, real &w, real &t, real &qv, real &rhoa, real &qr) = 0;
  virtual bool checkTime() = 0;
};

class BkgdStream : public BkgdAdapter {

 public:
  
  BkgdStream(const char *fname);
  ~BkgdStream();

  bool next(int &time, real &lat, real &lon, real &alt, real &u,
	    real &v, real &w, real &t, real &qv, real &rhoa, real &qr);

  bool checkTime() { return true; };
  
 private:

  std::ifstream *_stream;
  bool _valid;
};

// Abstract class.
// You need to pick a BkgdFArray, or BkgdCArray depending
// on whether your arrays are row-major or column-major

class BkgdArray : public BkgdAdapter {

 public:
  
  BkgdArray(int nx, int ny, int nsigma,
	    char *ctdg,		// current time group string
	    int delta,		// model time step on coarse grid
	    int iter,		// current iteration
	    float *sigmas,
	    float *latitude,	// 2d
	    float *longitude,
	    float *u1,		// 3d
	    float *v1,
	    float *w1,
	    float *th1,
	    float *p1);
  ~BkgdArray();

  bool next(int &time, real &lat, real &lon, real &alt, real &u,
	    real &v, real &w, real &t, real &qv, real &rhoa, real &qr);
  bool checkTime() { return false; };

 protected:
  
  virtual float item3d(float *a3d, int x, int y, int z) = 0;
  virtual float item2d(float *a2d, int x, int y)= 0;
    
  // TODO Do we need these 2?

  float _x_incr, _y_incr;

  // Time elements
  
  QDateTime _obTime;
  
  // 3d array dimentions

  int _xd, _yd, _zd;

  // 1D array

  float *_sigmas;
  
  // These are 2D

  float *_lat_2d;
  float *_long_2d;

  // These are 3D
  
  float* _u1_3d;
  float* _v1_3d;
  float* _w1_3d;
  float* _th1_3d;
  float* _p1_3d;

  // Iterator indices.
  
  int _curr_x, _curr_y, _curr_z;
};

// Redefine the pointer arithmetic to deal with column-major arrays (Fortran, Julia, ...)

class BkgdFArray : public BkgdArray {
  
 public:

 BkgdFArray(int nx, int ny, int nsigma,
	    char *ctdg,		// current time group string
	    int delta,		// model time step on coarse grid
	    int iter,		// current iteration
	    float *sigmas,
	    float *latitude,	// 2d
	    float *longitude,
	    float *u1,		// 3d
	    float *v1,
	    float *w1,
	    float *th1,
	    float *p1) : BkgdArray(nx, ny, nsigma,
				   ctdg,
				   delta,
				   iter,
				   sigmas,
				   latitude,
				   longitude,
				   u1,
				   v1,
				   w1,
				   th1,
				   p1) {};

  ~BkgdFArray();
  

 protected:
  
  float item3d(float *a3d, int x, int y, int z);
  float item2d(float *a2d, int x, int y);
};

// Redefine the pointer arithmetic to deal with row-major arrays (C, C++, c_float arrays from Python...)

class BkgdCArray : public BkgdArray {
  
 public:

 BkgdCArray(int nx, int ny, int nsigma,
	    char *ctdg,		// current time group string
	    int delta,		// model time step on coarse grid
	    int iter,		// current iteration
	    float *sigmas,
	    float *latitude,	// 2d
	    float *longitude,
	    float *u1,		// 3d
	    float *v1,
	    float *w1,
	    float *th1,
	    float *p1) : BkgdArray(nx, ny, nsigma,
				   ctdg,
				   delta,
				   iter,
				   sigmas,
				   latitude,
				   longitude,
				   u1,
				   v1,
				   w1,
				   th1,
				   p1) {};

  ~BkgdCArray();
  

 protected:
  
  float item3d(float *a3d, int x, int y, int z);
  float item2d(float *a2d, int x, int y);
};

// Place holder for Fractl input

class BkgdFractl : public BkgdAdapter {
  bool next(int &time, real &lat, real &lon, real &alt, real &u,
	    real &v, real &w, real &t, real &qv, real &rhoa, real &qr) { return true;};
  bool checkTime() { return true; };
};

#endif
