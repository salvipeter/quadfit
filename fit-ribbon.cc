#include <algorithm>
#include <fstream>
#include <iterator>

#include <geometry.hh>

#include <Eigen/Dense>

using namespace Geometry;

BSCurve interpolateWithZeroTangents(const Point2DVector &values) {
  size_t n = values.size() - 1;
  size_t m = n + 2;
  DoubleVector u; u.reserve(n + 1);
  for (const auto &v : values)
    u.push_back(v[0]);

  DoubleVector knots;

  std::fill_n(std::back_inserter(knots), 3, u.front());
  for (const auto &v : values)
    knots.push_back(v[0]);
  std::fill_n(std::back_inserter(knots), 3, u.back());

  BSBasis basis(3, knots);

  Eigen::MatrixXd A(n + 3, n + 3); A.setZero();
  for (size_t i = 0; i <= n; ++i) {
    size_t span = basis.findSpan(u[i]);
    DoubleVector coeff; basis.basisFunctions(span, u[i], coeff);
    for (size_t j = 0; j <= 3; ++j)
      A(i, span - 3 + j) = coeff[j];
  }
  A(n + 1, 1) = 1.0;
  A(n + 2, m - 1) = 1.0;

  Eigen::VectorXd b(n + 3);
  for (size_t i = 0; i <= n; ++i)
    b(i) = values[i][1];
  b(n + 1) = values.front()[1];
  b(n + 2) = values.back()[1];

  Eigen::VectorXd x;
  x = A.fullPivLu().solve(b);

  PointVector cpts;
  for (size_t i = 0; i <= m; ++i)
    cpts.emplace_back(x(i), 0, 0);

  return { 3, knots, cpts };
}

BSSurface splitU(BSSurface surface, double from, double to) {
  size_t p = surface.basisU().degree();
  surface = surface.insertKnotU(from, p);
  surface = surface.insertKnotU(to, p);
  const auto &knots = surface.basisU().knots();
  DoubleVector slice_knots;
  std::fill_n(std::back_inserter(slice_knots), p + 1, from);
  std::copy(std::upper_bound(knots.begin(), knots.end(), from),
            std::lower_bound(knots.begin(), knots.end(), to),
            std::back_inserter(slice_knots));
  std::fill_n(std::back_inserter(slice_knots), p + 1, to);
  PointVector cpts;
  size_t ncol = surface.numControlPoints()[1];
  std::copy_n(surface.controlPoints().begin() + (surface.basisU().findSpan(from) - p) * ncol,
              (slice_knots.size() - p - 1) * ncol, std::back_inserter(cpts));
  return { p, surface.basisV().degree(), slice_knots, surface.basisV().knots(), cpts };
}

BSSurface multiplySurfaceWithFunction(const BSSurface &surface, const BSCurve &function) {
  // TODO
  return surface;
}

std::vector<BSSurface> fitSlices(const BSSurface &ribbon, const Point2DVector &sh) {
  auto alpha = interpolateWithZeroTangents(sh);
  BSSurface sextic = multiplySurfaceWithFunction(ribbon, alpha);
  std::vector<BSSurface> result;
  for (size_t i = 1; i < sh.size(); ++i)
    result.push_back(splitU(sextic, sh[i-1][0], sh[i][0]));
  
  return result;
}


// I/O

template<typename T>
void writeType(std::ostream &os, T x) {
  os.write(reinterpret_cast<const char *>(&x), sizeof(T));
}

void writeVector(std::ostream &os, const Vector3D &v) {
  writeType<float>(os, v[0]);
  writeType<float>(os, v[1]);
  writeType<float>(os, v[2]);
}

void writeSTL(const std::vector<BSSurface> &surfaces, std::string filename, size_t resolution) {
  size_t n = surfaces.size();
  std::ofstream f(filename, std::ios::binary);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  std::string comment("Sampled B-spline surfaces.");
  comment.resize(80, ' ');
  f.write(comment.c_str(), 80);
  writeType<uint32_t>(f, n * resolution * resolution * 2); // # of all faces
  PointVector points;
  auto addTriangle = [&](size_t i0, size_t i1, size_t i2) {
                       const auto &p0 = points[i0], &p1 = points[i1], &p2 = points[i2];
                       auto normal = ((p1 - p0) ^ (p2 - p0)).normalize();
                       writeVector(f, normal);
                       writeVector(f, p0);
                       writeVector(f, p1);
                       writeVector(f, p2);
                       writeType<uint16_t>(f, 0);
                     };
  for (const auto &s : surfaces) {
    points.clear();
    for (size_t i = 0; i <= resolution; ++i) {
      double u = (double)i / resolution;
      u = s.basisU().low() * (1 - u) + s.basisU().high() * u;
      for (size_t j = 0; j <= resolution; ++j) {
        double v = (double)j / resolution;
        v = s.basisV().low() * (1 - v) + s.basisV().high() * v;
        points.push_back(s.eval(u, v));
      }
    }
    for (size_t i = 0; i < resolution; ++i) {
      size_t base = i * (resolution + 1);
      for (size_t j = 0; j < resolution; ++j) {
        addTriangle(base + j, base + j + 1, base + j + resolution + 1);
        addTriangle(base + j + resolution + 1, base + j + 1, base + j + resolution + 2);
      }
    }
  }
}

void writeControlNet(const std::vector<BSSurface> &surfaces, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  std::vector<size_t> start_indices;
  size_t index = 1;
  for (const auto &s : surfaces) {
    start_indices.push_back(index);
    for (const auto &p : s.controlPoints())
      f << "v " << p << std::endl;
    index += s.controlPoints().size();
  }
  for (size_t i = 0; i < surfaces.size(); ++i) {
    const auto &s = surfaces[i];
    size_t base = start_indices[i];
    auto [nu, nv] = s.numControlPoints();
    for (size_t j = 0; j < nu; ++j)
      for (size_t k = 1; k < nv; ++k)
        f << "l " << base + j * nv + k - 1 << ' ' << base + j * nv + k << std::endl;
    for (size_t j = 1; j < nu; ++j)
      for (size_t k = 0; k < nv; ++k)
        f << "l " << base + (j - 1) * nv + k << ' ' << base + j * nv + k << std::endl;
  }
}

void writeBoundaries(const std::vector<BSSurface> &surfaces, std::string filename,
                     size_t resolution) {
  std::ofstream f(filename);
  for (size_t n = 0; n < surfaces.size(); ++n) {
    const auto &s = surfaces[n];
    double u_lo = s.basisU().low(), u_hi = s.basisU().high(), u_d = u_hi - u_lo;
    double v_lo = s.basisV().low(), v_hi = s.basisV().high(), v_d = v_hi - v_lo;
    for (size_t i = 0; i <= resolution; ++i) {
      double t = (double)i / resolution;
      double u = u_lo + u_d * t, v = v_lo + v_d * t;
      f << "v " << s.eval(u_lo, v) << std::endl;
      f << "v " << s.eval(u, v_lo) << std::endl;
      f << "v " << s.eval(u_hi, v) << std::endl;
      f << "v " << s.eval(u, v_hi) << std::endl;
    }
    for (size_t j = 0; j < 4; ++j) {
      f << "l";
      for (size_t k = 0; k <= resolution; ++k)
        f << ' ' << 1 + (resolution + 1) * 4 * n + 4 * k + j;
      f << std::endl;
    }
  }
}


// Testing

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
