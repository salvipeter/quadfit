#include "connect-g1.hh"
#include "io.hh"
#include "multiply.hh"

using namespace Geometry;

BSSurface createSurface(const BSCurve &p, const BSCurve &d) {
  return multiplyBSplines(p.basis(), p.basis(),
                          [&](double u) { return p.eval(u); },
                          [&](double u) { return d.eval(u) - p.eval(u); });
}

int main(int argc, char **argv) {
  DoubleVector knots = { 0, 0, 0, 0, 1, 3, 4, 4, 4, 4 };
  BSCurve c(3, knots, {{0,0,0},{4,1,0},{8,1,0},{12,0,1},{16,0,0},{19,1,0}});
  BSCurve c1(3, knots, {{-2,3,0},{3,5,-2.8},{8,5,4},{12,4,1},{16,4,-18.0/11.0},{20,4,0}});
  BSCurve c2(3, knots, {{-1,-4,0},{4,-3,3},{8,-3,2},{13,-4,2},{17,-3,2},{21,-2,0}});
  size_t resolution = 50;
  std::vector<BSSurface> surfaces = { connectG1(c, c1, c2), connectG1(c, c2, c1) };
  surfaces[1].reverseV();
  writeSTL(surfaces, "/tmp/surfaces.stl", resolution);
  writeControlNet(surfaces, "/tmp/controls.obj");

  std::vector<BSSurface> surfaces_old = { createSurface(c, c1), createSurface(c, c2) };
  surfaces_old[1].reverseV();
  writeSTL(surfaces_old, "/tmp/surfaces-old.stl", resolution);
  writeControlNet(surfaces_old, "/tmp/controls_old.obj");
}
