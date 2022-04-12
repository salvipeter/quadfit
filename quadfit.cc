#include <cassert>
#include <fstream>

#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include "quadfit.hh"

using namespace Geometry;

static Point3D readPoint(std::istream &is) {
  Point3D p;
  is >> p[0] >> p[1] >> p[2];
  return p;
}

static BSCurve readBSpline(std::istream &is) {
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

std::string QuadFit::readPWGB(std::string filename) {
  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  std::string description;
  std::getline(f, description);

  // Part 1 : Ribbons
  size_t n_ribbons;
  f >> n_ribbons;
  ribbons.resize(n_ribbons);
  for (size_t i = 0; i < n_ribbons; ++i) {
    ribbons[i][0] = readBSpline(f);
    ribbons[i][1] = readBSpline(f);
  }

  // Part 2 : Curvenet segments
  size_t n_segments;
  f >> n_segments;
  segments.resize(n_segments);
  for (size_t i = 0; i < n_segments; ++i)
    segments[i] = readBSpline(f);

  // Part 3 : Quads
  size_t n_quads;
  f >> n_quads;
  quads.resize(n_quads);
  for (size_t i = 0; i < n_quads; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      auto &b = quads[i].boundaries[j];
      int s, r;
      f >> s >> r;
      b.segment = std::abs(s) - 1;
      b.reversed = s < 0;
      b.on_ribbon = r != 0;
      if (r > 0) {
        b.ribbon = r - 1;
        f >> b.s_min >> b.s_max;
      }
    }
  }

  // Part 4 : Sampled points
  for (size_t i = 0; i < n_quads; ++i) {
    size_t res;
    f >> res;
    quads[i].resolution = res;
    size_t n_samples = (res + 1) * (res + 1);
    quads[i].samples.resize(n_samples);
    for (size_t j = 0; j < n_samples; ++j)
      quads[i].samples[j] = readPoint(f);
  }

  return description;
}

static Matrix3x3 weingarten(const std::array<Vector3D, 5> &derivatives) {
  auto &[Su, Sv, Suv, Suu, Svv] = derivatives;
  double E = Su * Su, F = Su * Sv, G = Sv * Sv;
  Vector3D n = (Su ^ Sv).normalize();
  double L = n * Suu, M = n * Suv, N = n * Svv;
  Eigen::Matrix2d I, II;
  I << E, F, F, G;
  II << L, M, M, N;
  Eigen::Matrix<double, 3, 2> J;
  J.col(0) = Eigen::Map<const Eigen::Vector3d>(Su.data());
  J.col(1) = Eigen::Map<const Eigen::Vector3d>(Sv.data());
  Eigen::Matrix<double, 2, 3> J1 = I.inverse() * J.transpose();
  Eigen::Matrix3d W = J1.transpose() * II * J1;
  return Matrix3x3(W.data());
}

Matrix3x3 QuadFit::fitCurvature(const Point3D &p, const std::set<size_t> &quad_indices) const {
  Vector3D Su, Sv, Suv, Suu, Svv;

  if (quad_indices.size() == 1) {
    // Corner - everything is known
    const auto &b = quads[*quad_indices.begin()].boundaries;
    size_t i = 0;
    while (i < 4) {
      const auto &cp = segments[b[i].segment].controlPoints();
      if ((b[i].reversed && (cp.front() - p).norm() < epsilon) ||
          (!b[i].reversed && (cp.back() - p).norm() < epsilon))
        break;
      i++;
    }
    assert(i < 4 && "Invalid corner quad");
    size_t ip = (i + 1) % 4;
    assert(b[i].on_ribbon && b[ip].on_ribbon && "Corner quad without two ribbons");
    const auto &rib_left = ribbons[b[i].ribbon];
    const auto &rib_right = ribbons[b[ip].ribbon];
    double s_left = b[i].s_max, s_right = b[ip].s_min;
    VectorVector der;
    rib_left[0].eval(s_left, 2, der);
    Sv = -der[1];
    Svv = der[2];
    rib_right[0].eval(s_right, 2, der);
    Su = der[1];
    Suu = der[2];
    rib_left[1].eval(s_left, 1, der);
    Suv = (-der[1] - Sv) * 3; // assumes cubic Bezier in the cross direction ???
    rib_right[1].eval(s_right, 1, der);
    Suv = (Suv + (der[1] - Su) * 3) / 2; // take the mean ???

  } else if (quad_indices.size() == 2) {
    // Boundary - Svv is unknown

  } else {
    // Simple fit
  }

  return weingarten({ Su, Sv, Suv, Suu, Svv });
}

[[maybe_unused]]
static Vector3D weingartenNormal(const Matrix3x3 &m) {
  Eigen::EigenSolver<Eigen::Matrix3d>
    es(Eigen::Map<const Eigen::Matrix3d>(m.data()), true);
  auto n = es.eigenvectors().col(0);
  return { n(0).real(), n(1).real(), n(2).real() };
}

[[maybe_unused]]
static double weingartenPrincipalCurvature(const Matrix3x3 &m, size_t index) {
  Eigen::EigenSolver<Eigen::Matrix3d>
    es(Eigen::Map<const Eigen::Matrix3d>(m.data()), false);
  return es.eigenvalues()(index).real();
}

[[maybe_unused]]
static Vector3D weingartenPrincipalDirection(const Matrix3x3 &m, size_t index) {
  Eigen::EigenSolver<Eigen::Matrix3d>
    es(Eigen::Map<const Eigen::Matrix3d>(m.data()), true);
  auto d = es.eigenvectors().col(index);
  return { d(0).real(), d(1).real(), d(2).real() };
}

[[maybe_unused]]
static void writeVertexCurvatures(const std::vector<Point3D> &vertices,
                                  const std::vector<Matrix3x3> &curvatures,
                                  std::string filename) {
  std::ofstream f(filename);
  f << "# vtk DataFile Version 2.0" << std::endl;
  f << "Vertices with principal curvature values & directions" << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET POLYDATA" << std::endl;
  f << "POINTS " << vertices.size() << " float" << std::endl;
  for (const auto &v : vertices) {
    f << v << std::endl;
  }
  f << "POINT_DATA " << curvatures.size() << std::endl;
  f << "NORMALS normal float" << std::endl;
  for (const auto &m : curvatures) {

    f << weingartenNormal(m) << std::endl;
  }
  f << "SCALARS k1 float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for (const auto &m : curvatures) {
    f << weingartenPrincipalCurvature(m, 1) << std::endl;
  }
  f << "SCALARS k2 float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for (const auto &m : curvatures) {
    f << weingartenPrincipalCurvature(m, 2) << std::endl;
  }
  f << "NORMALS d1 float" << std::endl;
  for (const auto &m : curvatures) {
    f << weingartenPrincipalDirection(m, 1) << std::endl;
  }
  f << "NORMALS d2 float" << std::endl;
  for (const auto &m : curvatures) {
    f << weingartenPrincipalDirection(m, 2) << std::endl;
  }
}

std::vector<BSSurface> QuadFit::fit() {
  std::vector<BSSurface> result(quads.size());

  // Topology & geometry structures
  std::vector<Point3D> vertices;
  std::vector<std::set<size_t>> quad_indices;       // vertex -> quads
  std::vector<std::pair<size_t, size_t>> endpoints; // segment -> from_vertex, to_vertex
  std::vector<Matrix3x3> curvatures;                // vertex -> embedded Weingarten matrix

  // Build topology

  auto vertexIndex = [&](const Point3D &p) {
    for (size_t i = 0; i < vertices.size(); ++i)
      if ((p - vertices[i]).norm() < epsilon)
        return i;
    vertices.push_back(p);
    return vertices.size() - 1;
  };
  for (size_t i = 0; i < segments.size(); ++i) {
    size_t j0 = vertexIndex(segments[i].controlPoints().front());
    size_t j1 = vertexIndex(segments[i].controlPoints().back());
    endpoints.emplace_back(j0, j1);
  }
  quad_indices.resize(vertices.size());
  for (size_t i = 0; i < quads.size(); ++i)
    for (size_t j = 0; j < 4; ++j) {
      size_t s = quads[i].boundaries[j].segment;
      quad_indices[endpoints[s].first].insert(i);
      quad_indices[endpoints[s].second].insert(i);
    }

  // Fit vertex curvatures

  curvatures.resize(vertices.size());
  for (size_t i = 0; i < vertices.size(); ++i)
    curvatures[i] = fitCurvature(vertices[i], quad_indices[i]);
  writeVertexCurvatures(vertices, curvatures, "/tmp/curvatures.vtk");

  // Simple C0 fit

  for (size_t i = 0; i < quads.size(); ++i) {
    auto &b = quads[i].boundaries;
    PointVector cpts(16);
    for (size_t j = 0; j < 4; ++j) {
      auto getCP = [&](size_t k) {
        size_t j1 = b[k].reversed ? 3 - j : j;
        return segments[b[k].segment].controlPoints()[j1];
      };
      cpts[j] = getCP(0);
      cpts[j*4] = getCP(1);
      cpts[12+j] = getCP(2);
      cpts[3+j*4] = getCP(3);
    }
    auto setTwist = [&](size_t p, size_t d1, size_t d2, size_t t) {
      cpts[t] = cpts[d1] + (cpts[d2] - cpts[p]);
    };
    setTwist(0, 1, 4, 5);
    setTwist(3, 2, 7, 6);
    setTwist(12, 8, 13, 9);
    setTwist(15, 11, 14, 10);
    result[i] = BSSurface(3, 3, cpts);
  }

  return result;
}
