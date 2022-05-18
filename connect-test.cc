#include "connect-g1.hh"
#include "io.hh"

using namespace Geometry;

int main(int argc, char **argv) {
  DoubleVector knots = { 0, 0, 0, 0, 1, 3, 4, 4, 4, 4 };
  BSCurve c(3, knots, {{0,0,0},{4,1,0},{8,1,0},{12,0,1},{16,0,1},{19,1,0}});
  BSCurve c1(3, knots, {{-2,3,0},{3,5,1},{8,5,2},{12,4,1},{16,4,0},{20,4,0}});
  BSCurve c2(3, knots, {{-1,-4,0},{4,-3,1},{8,-3,2},{13,-4,1},{17,-3,1},{21,-2,0}});
  size_t resolution = 50;
  std::vector<BSSurface> surfaces = { connectG1(c, c1, c2), connectG1(c, c2, c1) };
  surfaces[1].reverseV();
  writeSTL(surfaces, "/tmp/surfaces.stl", resolution);
  writeControlNet(surfaces, "/tmp/controls.obj");
}
