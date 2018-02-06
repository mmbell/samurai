#ifndef BKGD_ADAPTER_H
#define BKGD_ADAPTER_H

#include <fstream>
#include <QDateTime>

#include "precision.h"	// for 'real' typedef

// This is an abstract class to support the multiple ways Background observations can
// be loaded into samurai.
//
// Right now we support loading from
//   a file (samurai_Background.in) The format is documented on the Samurai github wiki
//   a set of arrays. This is a new format to support COAMPS. Interface is a work in progress.
//

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

 private:
  
  float item3d(float *a3d, int x, int y, int z);
  float item2d(float *a2d, int x, int y);
    
  // TODO Do we need these 2?

  float _x_incr, _y_incr;

  // Time elements
  
  QDateTime _obTime;
  
  // 3d array dimentions

  int _xd, _yd, _zd;

  // 1D array from Fortran

  float *_sigmas;
  
  // These are 2D from Fortran (column dominant)

  float *_lat_2d;
  float *_long_2d;

  // These are 3D arrays from Fortran (column dominant)
  
  float* _u1_3d;
  float* _v1_3d;
  float* _w1_3d;
  float* _th1_3d;
  float* _p1_3d;

  // Iterator indices.
  
  int _curr_x, _curr_y, _curr_z;
};

#endif
