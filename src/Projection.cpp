#include "Projection.h"

// We only use one of the projections, but it isn't worth optimizing
// and creating only one.

Projection::Projection(ProjectionType t) :
  tm1(GeographicLib::TransverseMercatorExact::UTM()),
  tm2(GeographicLib::LambertConformalConic::Mercator())
{
  projection_type = t;
}

void Projection::setProjection(ProjectionType t)
{
  projection_type = t;
}

void Projection::Forward (real lon0, real lat, real lon,
			  real &x, real &y) const 
{
  switch(projection_type) {
  case LAMBERT_CONFORMAL_CONIC:
    tm2.Forward(lon0, lat, lon, x, y);
    break;
  case TRANSVERSE_MERCATOR_EXACT:
  default:
    tm1.Forward(lon0, lat, lon, x, y);
    break;
  }
}

void Projection::Reverse (real lon0, real x, real y,
			  real &lat, real &lon) const
{
  switch(projection_type) {
  case LAMBERT_CONFORMAL_CONIC:
    tm2.Reverse(lon0, x, y, lat, lon);
    break;
  case TRANSVERSE_MERCATOR_EXACT:
  default:
    tm1.Reverse(lon0, x, y, lat, lon);
    break;
  }
}

