#pragma once

#include <functional>

#include <geometry.hh>

Geometry::BSSurface
multiplyBSplines(const Geometry::BSBasis &basis1,
                 const Geometry::BSBasis &basis2,
                 const std::function<Geometry::Point3D(double)> &point,
                 const std::function<Geometry::Vector3D(double)> &derivative);
