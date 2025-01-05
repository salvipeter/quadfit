#include <fstream>

#include <geometry.hh>

using namespace Geometry;

Point3D readPoint(std::istream &is) {
  Point3D p;
  is >> p[0] >> p[1] >> p[2];
  return p;
}

BSCurve readBSpline(std::istream &is) {
  size_t degree, n_knots, n_cpts;
  DoubleVector knots;
  PointVector cpts;
  is >> degree;
  is >> n_knots;
  knots.resize(n_knots);
  for (size_t i = 0; i < n_knots; ++i)
    is >> knots[i];
  is >> n_cpts;
  cpts.resize(n_cpts);
  for (size_t i = 0; i < n_cpts; ++i)
    cpts[i] = readPoint(is);
  return { degree, knots, cpts };
}

void writeCurves(const std::vector<BSCurve> &curves, std::string filename, size_t resolution) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  size_t index = 0;
  for (const auto &c : curves) {
    double lo = c.basis().low(), hi = c.basis().high();
    for (size_t i = 0; i < resolution; ++i) {
      double u = lo + (double)i / (resolution - 1) * (hi - lo);
      f << "v " << c.eval(u) << std::endl;
    }
    f << 'l';
    for (size_t i = 0; i < resolution; ++i)
      f << ' ' << ++index;
    f << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc != 3 && argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input.bcs> <output.obj> [resolution]" << std::endl;
    return 1;
  }
  size_t resolution = 100;
  if (argc == 4)
    resolution = std::atoi(argv[3]);
  std::ifstream f(argv[1]);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  std::vector<BSCurve> curves;
  size_t n;
  f >> n;
  for (size_t i = 0; i < n; ++i)
    curves.push_back(readBSpline(f));
  f.close();
  writeCurves(curves, argv[2], resolution);
}
