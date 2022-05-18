#pragma once

#include <geometry.hh>

void writeSTL(const std::vector<Geometry::BSSurface> &surfaces,
              std::string filename, size_t resolution);

void writeControlNet(const std::vector<Geometry::BSSurface> &surfaces,
                     std::string filename);

void writeBoundaries(const std::vector<Geometry::BSSurface> &surfaces,
                     std::string filename, size_t resolution);
