#pragma once

#include <geometry.hh>

void writeCurves(const std::vector<Geometry::BSCurve> &curve,
                 std::string filename, size_t resolution);

void writeSTL(const std::vector<Geometry::BSSurface> &surfaces,
              std::string filename, size_t resolution);

void writeControlNet(const std::vector<Geometry::BSSurface> &surfaces,
                     std::string filename);

void writeBoundaries(const std::vector<Geometry::BSSurface> &surfaces,
                     std::string filename, size_t resolution);
