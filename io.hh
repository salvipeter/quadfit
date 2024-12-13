#pragma once

#include <geometry.hh>

void writeCurves(const std::vector<Geometry::BSCurve> &curve,
                 std::string filename, size_t resolution);

void writeControlPolygon(const std::vector<Geometry::BSCurve> &curves, std::string filename);

void writeQDS(const std::vector<Geometry::BSSurface> &surfaces, std::string filename);

void writeSTL(const std::vector<Geometry::BSSurface> &surfaces,
              std::string filename, size_t resolution);

void writeControlNet(const std::vector<Geometry::BSSurface> &surfaces,
                     std::string filename);

void writeBoundaries(const std::vector<Geometry::BSSurface> &surfaces,
                     std::string filename, size_t resolution);
