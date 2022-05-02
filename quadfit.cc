#include <cassert>
#include <fstream>

#include "jet-wrapper.hh"

#include "quadfit.hh"

using namespace Geometry;
using JetWrapper::JetData;

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

[[maybe_unused]]
static void writeVertexCurvatures(const std::vector<Point3D> &vertices,
                                  const std::vector<JetData> &jet,
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
  f << "POINT_DATA " << jet.size() << std::endl;
  f << "NORMALS normal float" << std::endl;
  for (const auto &jd : jet) {

    f << jd.normal << std::endl;
  }
  f << "SCALARS k1 float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for (const auto &jd : jet) {
    f << jd.k_min << std::endl;
  }
  f << "SCALARS k2 float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for (const auto &jd : jet) {
    f << jd.k_max << std::endl;
  }
  f << "NORMALS d1 float" << std::endl;
  for (const auto &jd : jet) {
    f << jd.d_min << std::endl;
  }
  f << "NORMALS d2 float" << std::endl;
  for (const auto &jd : jet) {
    f << jd.d_max << std::endl;
  }
}

std::vector<BSSurface> QuadFit::initialFit() const {
  std::vector<BSSurface> result(quads.size());

  for (size_t i = 0; i < quads.size(); ++i) {
    auto &b = quads[i].boundaries;
    PointVector cpts(16);
    for (size_t j = 0; j < 4; ++j) {
      auto getCP = [&](size_t k) {
        const auto &curve = segments[b[k].segment];
        bool rev = (b[k].reversed && j < 2) || (!b[k].reversed && j >= 2);
        bool deriv = j == 1 || j == 2;
        double u = rev ? curve.basis().high() : curve.basis().low();
        VectorVector der;
        auto p = curve.eval(u, deriv ? 1 : 0, der);
        if (deriv) {
          double interval = curve.basis().high() - curve.basis().low();
          return p + der[1] / 3 * (rev ? -1 : 1) * interval;
        }
        return p;
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

static BSSurface elevateBezierU(const BSSurface &surface) {
  size_t du = surface.basisU().degree();
  size_t dv = surface.basisV().degree();
  const auto &old = surface.controlPoints();
  PointVector cpts((du + 2) * (dv + 1));
  size_t index = 0;
  for (size_t i = 0; i <= du + 1; ++i)
    for (size_t j = 0; j <= dv; ++j) {
      if (i == 0)
        cpts[index] = old[index];
      else if (i > du)
        cpts[index] = old[index-dv-1];
      else {
        double alpha = (double)i / (du + 1);
        cpts[index] = old[index] * (1 - alpha) + old[index-dv-1] * alpha;
      }
      index++;
    }
  return { du + 1, dv, cpts };
}

static BSSurface elevateBezierV(const BSSurface &surface) {
  size_t du = surface.basisU().degree();
  size_t dv = surface.basisV().degree();
  const auto &old = surface.controlPoints();
  PointVector cpts((du + 1) * (dv + 2));
  size_t index = 0;
  for (size_t i = 0; i <= du; ++i)
    for (size_t j = 0; j <= dv + 1; ++j) {
      if (j == 0)
        cpts[index] = old[index-i];
      else if (j > dv)
        cpts[index] = old[index-i-1];
      else {
        double alpha = (double)j / (dv + 1);
        cpts[index] = old[index-i] * (1 - alpha) + old[index-i-1] * alpha;
      }
      index++;
    }
  return { du, dv + 1, cpts };
}

static BSSurface elevateBezier(const BSSurface &surface, size_t target) {
  BSSurface result = surface;
  while (result.basisU().degree() < target)
    result = elevateBezierU(result);
  while (result.basisV().degree() < target)
    result = elevateBezierV(result);
  return result;
}

std::vector<BSSurface> QuadFit::fit() {
  // Topology & geometry structures
  std::vector<Point3D> vertices;
  std::vector<std::set<size_t>> quad_indices;       // vertex -> quads
  std::vector<std::pair<size_t, size_t>> endpoints; // segment -> from_vertex, to_vertex

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

  // 0. Fit vertex curvatures
  PointVector all_points;
  for (const auto &q : quads)
    all_points.insert(all_points.end(), q.samples.begin(), q.samples.end());
  auto jet = JetWrapper::fit(vertices, JetWrapper::Nearest(all_points));
  // TODO: these should be fixed at non-corner ribbon vertices
  writeVertexCurvatures(vertices, jet, "/tmp/curvatures.vtk");

  // 1. Simple C0 fit
  // - cubic Bezier surfaces
  // - boundaries by the boundary curve endpoints & derivatives
  // - twists by parallelogram rule
  auto result = initialFit();

  // 2. Compute better tangents
  // - by the ribbons where available
  // - projected to the vertices' normal plane
  constexpr std::array<size_t, 8> vertex_cps = { 0, 3, 0, 12, 12, 15, 3, 15 };
  constexpr std::array<size_t, 8> tangent_cps = { 1, 2, 4, 8, 13, 14, 7, 11 };
  for (size_t i = 0; i < quads.size(); ++i) {
    auto &cpts = result[i].controlPoints();
    for (size_t side = 0; side < 4; ++side) {
      const auto &b = quads[i].boundaries[side];
      if (b.on_ribbon) {
        const auto &curve = ribbons[b.ribbon][0];
        auto matchRibbon = [&](size_t v, size_t t, double interval, const Vector3D &d) {
          cpts[t] = cpts[v] + d / 3 * interval;
        };
        VectorVector der;
        curve.eval(b.s_min, 1, der);
        matchRibbon(vertex_cps[2*side], tangent_cps[2*side], b.s_max - b.s_min, der[1]);
        curve.eval(b.s_max, 1, der);
        matchRibbon(vertex_cps[2*side+1], tangent_cps[2*side+1], b.s_max - b.s_min, -der[1]);
      } else {
        auto projectToPlane = [&](size_t v, size_t t, const Vector3D &n) {
          cpts[t] += n * (n * (cpts[v] - cpts[t]));
        };
        const auto &n1 = jet[endpoints[b.segment].first].normal;
        const auto &n2 = jet[endpoints[b.segment].second].normal;
        projectToPlane(vertex_cps[2*side], tangent_cps[2*side], b.reversed ? n2 : n1);
        projectToPlane(vertex_cps[2*side+1], tangent_cps[2*side+1], b.reversed ? n1 : n2);
      }
    }
  }

  // 3. Compute better twists

  // 4. Degree elevation to sextic
  for (auto &s : result)
    s = elevateBezier(s, 6);

  // 5. Compute better curvatures

  // 6. Fit sampled points using inner control points

  // 7. Set G1 continuity by direction blends

  // 8. Insert knots & refine boundaries

  // 9. Re-fit inner control points

  return result;
}
