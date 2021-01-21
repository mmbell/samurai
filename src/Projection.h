/*
 * Projection.h
 * samurai
 *
 * Copyright?
 *
 * This is an abstaction of the various (2 so far) of the
 * geographiclib projections we use in this program.
 *
 * Only one version of Forward and Reverse are supported, since
 * these are the only ones we use.
 *
 * Unfortunately, geographiclib doesn't offer an abstract
 * Projection class, which makes it hard to use one or the
 * other depending on configuration.
 */

#ifndef PROJECTION_H
#define PROJECTION_H

#include <euclid/GeographicLib/TransverseMercatorExact.hpp>
#include <euclid/GeographicLib/LambertConformalConic.hpp>
#include "precision.h"

class Projection
{

 public:

  // The kind of projections we currently support
  enum ProjectionType {
    UNSET,
    TRANSVERSE_MERCATOR_EXACT,
    LAMBERT_CONFORMAL_CONIC,
  };

  // constructor
  Projection(ProjectionType proj = TRANSVERSE_MERCATOR_EXACT);

  // Set the projection type
  void setProjection(ProjectionType t);
  void Forward (real lon0, real lat, real lon, real &x, real &y) const;
  void Reverse (real lon0, real x, real y, real &lat, real &lon) const;

 private:

  ProjectionType projection_type;
  GeographicLib::TransverseMercatorExact	tm1;
  GeographicLib::LambertConformalConic		tm2;
};

#endif /* PROJECTION_H */

