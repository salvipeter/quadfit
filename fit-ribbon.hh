#pragma once

#include <geometry.hh>

std::vector<Geometry::BSSurface> fitSlices(const Geometry::BSSurface &ribbon,
                                           const Geometry::Point2DVector &sh);
