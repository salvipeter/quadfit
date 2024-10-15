#pragma once

#include <geometry.hh>

Geometry::Point3D
closestPoint(const Geometry::BSCurve &curve,
             const Geometry::Point3D &point, double &u,
             size_t max_iteration, double distance_tol, double cosine_tol);

Geometry::Point3D
closestPoint(const Geometry::BSSurface &surface,
             const Geometry::Point3D &point, double &u, double &v,
             size_t max_iteration, double distance_tol, double cosine_tol);
