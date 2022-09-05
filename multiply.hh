#pragma once

#include <functional>

#include <geometry.hh>

Geometry::BSSurface
multiplyBSplines(const Geometry::BSBasis &basis1,
                 const Geometry::BSBasis &basis2,
                 const std::function<Geometry::Point3D(double)> &point,
                 const std::function<Geometry::Vector3D(double)> &derivative);

// Point & CrossDer are in the same basis, Scaling is in another
// This computes a surface that is linear in V, interpolates Point in V = 0,
// and interpolates Point + CrossDer * Scaling[0] + Point' * Scaling[1] in V = 1.
Geometry::BSSurface
multiplyBSplinesByBezier(Geometry::BSCurve point,
                         Geometry::BSCurve crossder,
                         Geometry::BSCurve scaling);
