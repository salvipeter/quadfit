#include "fit-ribbon.hh"
#include "io.hh"

using namespace Geometry;

int main(int argc, char **argv) {
  BSSurface ribbon(3, 1, { 0, 0, 0, 0, 0.3, 0.5, 0.6, 1, 1, 1, 1 }, { 0, 0, 1, 1 },
                   { { 0, 0, 0 }, { 0, 4, 1 },
                     { 3, 1, 1 }, { 4, 5, 2 },
                     { 7, 2, 1 }, { 8, 5, 1 },
                     {10, 2, 2 }, {11, 6, 1 },
                     {14, 0, 1 }, {16, 4, 2 },
                     {17, 0, 0 }, {18, 4, 1 },
                     {20, 2, 1 }, {21, 5, 2 } });
  Point2DVector sh = { { 0, 0.5 }, { 0.3, 0.6 }, { 0.4, 0.7 }, { 1, 0.9 } };
  size_t resolution = 50;
  auto surfaces = fitSlices(ribbon, sh);
  writeSTL({ ribbon }, "/tmp/ribbon.stl", resolution);
  writeControlNet({ ribbon }, "/tmp/ribbon-controls.obj");
  writeSTL(surfaces, "/tmp/surfaces.stl", resolution);
  writeControlNet(surfaces, "/tmp/controls.obj");
  writeBoundaries(surfaces, "/tmp/boundaries.obj", resolution);
}
