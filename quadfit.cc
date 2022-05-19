#include <cassert>
#include <fstream>

#include "jet-wrapper.hh"

#include "quadfit.hh"

#include "bspline-fit.hh"

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
        f >> b.s0 >> b.s1 >> b.h0 >> b.h1;
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

[[maybe_unused]]
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

  writeVertexCurvatures(vertices, jet, "/tmp/curvatures.vtk"); // (just for debugging)


  // 1. Simple C0 fit
  // - cubic Bezier surfaces
  // - boundaries by the boundary curve endpoints & derivatives
  // - twists by parallelogram rule
  auto result = initialFit();


  // 2. Fit ribbons


  // 3. Compute better first derivatives for the quad boundary curves
  // - by the fitted ribbons where available
  // - in the interior: project to the vertices' normal plane
  constexpr std::array<size_t, 8> vertex_cps = { 0, 3, 0, 12, 12, 15, 3, 15 };
  constexpr std::array<size_t, 8> tangent_cps = { 1, 2, 4, 8, 13, 14, 7, 11 };
  for (size_t i = 0; i < quads.size(); ++i) {
    auto &cpts = result[i].controlPoints();
    for (size_t side = 0; side < 4; ++side) {
      const auto &b = quads[i].boundaries[side];
      const auto &bp = quads[i].boundaries[(side+3)%4];
      const auto &bn = quads[i].boundaries[(side+1)%4];
      const auto &curve = ribbons[b.ribbon][0];
      auto matchRibbon = [&](size_t v, size_t t, double interval, const Vector3D &d) {
        cpts[t] = cpts[v] + d / 3 * interval;
      };
      VectorVector der;
      if (b.on_ribbon) {
        curve.eval(b.s0, 1, der);
        matchRibbon(vertex_cps[2*side], tangent_cps[2*side], b.s1 - b.s0, der[1]);
        curve.eval(b.s1, 1, der);
        matchRibbon(vertex_cps[2*side+1], tangent_cps[2*side+1], b.s1 - b.s0, -der[1]);
      } else {
        auto projectToPlane = [&](size_t v, size_t t, const Vector3D &n) {
          cpts[t] += n * (n * (cpts[v] - cpts[t]));
        };
        const auto &n1 = jet[endpoints[b.segment].first].normal;
        const auto &n2 = jet[endpoints[b.segment].second].normal;
        bool replace_start = ((bp.on_ribbon && (side == 1 || side == 2)) ||
                              (bn.on_ribbon && (side == 0 || side == 3)));
        bool replace_end   = ((bp.on_ribbon && (side == 0 || side == 3)) ||
                              (bn.on_ribbon && (side == 1 || side == 2)));
        if (replace_start) {
          auto v = cpts[vertex_cps[2*side]], t = cpts[tangent_cps[2*side]];
          Point3D p;
          if (side < 2) {
            const auto &B = side == 0 ? bn : bp;
            p = ribbons[B.ribbon][1].eval(B.s0);
          } else {
            const auto &B = side == 3 ? bn : bp;
            p = ribbons[B.ribbon][1].eval(B.s1);
          }
          cpts[tangent_cps[2*side]] = v + (p - v).normalized() * (t - v).norm(); // or: p
        } else
          projectToPlane(vertex_cps[2*side], tangent_cps[2*side], b.reversed ? n2 : n1);
        if (replace_end) {
          auto v = cpts[vertex_cps[2*side+1]], t = cpts[tangent_cps[2*side+1]];
          Point3D p;
          if (side < 2) {
            const auto &B = side == 1 ? bn : bp;
            p = ribbons[B.ribbon][1].eval(B.s0);
          } else {
            const auto &B = side == 2 ? bn : bp;
            p = ribbons[B.ribbon][1].eval(B.s1);
          }
          cpts[tangent_cps[2*side+1]] = v + (p - v).normalized() * (t - v).norm(); // or: p
        } else
          projectToPlane(vertex_cps[2*side+1], tangent_cps[2*side+1], b.reversed ? n1 : n2);
      }
    }
  }


  // 4. Compute better second derivatives for the quad boundary curves


  // 5. Compute twist vectors
  constexpr std::array<size_t, 16> corner_cps =
    { 0, 1, 4, 5, 12, 8, 13, 9, 15, 14, 11, 10, 3, 7, 2, 6 };
  for (size_t i = 0; i < quads.size(); ++i) {
    auto &cpts = result[i].controlPoints();
    for (size_t side = 0; side < 4; ++side) {
      size_t p = corner_cps[4*side], t1 = corner_cps[4*side+1];
      size_t t2 = corner_cps[4*side+2], tw = corner_cps[4*side+3];
      const auto &b = quads[i].boundaries[side];
      const auto &bn = quads[i].boundaries[(side+1)%4];
      if (b.on_ribbon && bn.on_ribbon) {
        VectorVector der0, der1;
        if (side == 0 || side == 3) {
          ribbons[b.ribbon][0].eval(b.s0, 1, der0);
          ribbons[b.ribbon][1].eval(b.s0, 1, der1);
        } else {
          ribbons[b.ribbon][0].eval(b.s1, 1, der0);
          ribbons[b.ribbon][1].eval(b.s1, 1, der1);
        }
        auto twist = (der1[1] - der0[1]) * (b.s1 - b.s0);
        cpts[tw] = twist / 9 - cpts[p] + cpts[t1] + cpts[t2];
      } else {
        // TODO
      }
    }
  }


  // 6. Compute inner boundary ribbons


  // 7. Fill sextic patches


  // 8. Fit sampled points using inner control points
  // for (size_t i = 0; i < quads.size(); ++i)
  //   bsplineFit(result[i], quads[i].resolution, quads[i].samples, fixC0Twist, 0.1);

  return result;
}
