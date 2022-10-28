#include <cmath>
#include <fstream>
#include <iterator>

#include "io.hh"

using namespace Geometry;

void writeCurves(const std::vector<BSCurve> &curves, std::string filename, size_t resolution) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  f << "# vtk DataFile Version 2.0" << std::endl;
  f << "Curves with  curvature values & normal directions" << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET POLYDATA" << std::endl;

  f << "POINTS " << curves.size() * (resolution + 1) << " float" << std::endl;
  for (const auto &curve : curves) {
    double lo = curve.basis().low(), hi = curve.basis().high();
    for (size_t i = 0; i <= resolution; ++i) {
      double u = lo + (double)i / resolution * (hi - lo);
      f << curve.eval(u) << std::endl;
    }
  }

  f << "LINES " << curves.size() << " " << curves.size() * (resolution + 2) << std::endl;
  for (size_t j = 0; j < curves.size(); ++j) {
    f << resolution + 1;
    for (size_t i = 0; i <= resolution; ++i)
      f << ' ' << j * (resolution + 1) + i;
    f << std::endl;
  }

  f << "POINT_DATA " << curves.size() * (resolution + 1) << std::endl;
  f << "NORMALS normal float" << std::endl;
  for (const auto &curve : curves) {
    double lo = curve.basis().low(), hi = curve.basis().high();
    for (size_t i = 0; i <= resolution; ++i) {
      double u = lo + (double)i / resolution * (hi - lo);
      VectorVector der;
      curve.eval(u, 2, der);
      f << (der[0] ^ (der[0] ^ der[1])).normalize() << std::endl;
    }
  }

  f << "SCALARS curvature float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for (const auto &curve : curves) {
    double lo = curve.basis().low(), hi = curve.basis().high();
    for (size_t i = 0; i <= resolution; ++i) {
      double u = lo + (double)i / resolution * (hi - lo);
      VectorVector der;
      curve.eval(u, 2, der);
      f << (der[0] ^ der[1]).norm() / std::pow(der[0].norm(), 3) << std::endl;
    }
  }
}

void writeQDS(const std::vector<Geometry::BSSurface> &surfaces, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  f << surfaces.size() << std::endl;
  for (const auto &s : surfaces) {
    f << s.basisU().degree() << ' ' << s.basisV().degree() << std::endl;
    f << s.basisU().knots().size() << ' ';
    std::copy(s.basisU().knots().begin(), s.basisU().knots().end(),
              std::ostream_iterator<double>(f, " "));
    f << std::endl;
    f << s.basisV().knots().size() << ' ';
    std::copy(s.basisV().knots().begin(), s.basisV().knots().end(),
              std::ostream_iterator<double>(f, " "));
    f << std::endl;
    std::copy(s.controlPoints().begin(), s.controlPoints().end(),
              std::ostream_iterator<Point3D>(f, "\n"));
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
