#include <fstream>

#include "io.hh"

using namespace Geometry;

void writeCurves(const std::vector<BSCurve> &curves, std::string filename, size_t resolution) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  for (const auto &curve : curves) {
    double lo = curve.basis().low(), hi = curve.basis().high();
    for (size_t i = 0; i <= resolution; ++i) {
      double u = lo + (double)i / resolution * (hi - lo);
      f << "v " << curve.eval(u) << std::endl;
    }
  }
  for (size_t j = 0; j < curves.size(); ++j) {
    f << "l";
    for (size_t i = 0; i <= resolution; ++i)
      f << ' ' << j * (resolution + 1) + i + 1;
    f << std::endl;
  }
}

template<typename T>
static void writeType(std::ostream &os, T x) {
  os.write(reinterpret_cast<const char *>(&x), sizeof(T));
}

static void writeVector(std::ostream &os, const Vector3D &v) {
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
