#pragma once

#include <functional>

#include <geometry.hh>

void bsplineFit(Geometry::BSSurface &surface, size_t resolution,
                const Geometry::PointVector &samples,
                const std::function<bool(size_t,size_t)> &isFixed, double smoothness);
