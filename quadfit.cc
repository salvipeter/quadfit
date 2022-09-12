#include <cassert>
#include <cmath>
#include <fstream>

#include <jet-wrapper.hh>

#include "quadfit.hh"

#include "bspline-fit.hh"
#include "connect-g1.hh"
#include "fit-ribbon.hh"
#include "io.hh"

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

std::vector<std::pair<size_t,size_t>> QuadFit::ribbonSegments(size_t index) const {
  std::vector<std::pair<size_t,size_t>> result;
  const auto &r = ribbons[index][0];
  auto p = r.controlPoints().front(), last = r.controlPoints().back();
  auto samePoint = [](const Point3D &a, const Point3D &b) {
    return (a - b).normSqr() == 0; // assumes numerical equivalence!
  };
  while (!samePoint(p, last)) {
    bool found = false;
    for (size_t i = 0; i < quads.size(); ++i) {
      const auto &q = quads[i];
      for (size_t j = 0; j < 4; ++j) {
        const auto &b = q.boundaries[j];
        if (b.on_ribbon && b.ribbon == index) {
          const auto &cp = segments[b.segment].controlPoints();
          if (samePoint(cp.front(), p) &&
              (result.empty() || result.back().first != i)) {
            result.push_back({ i, j });
            p = cp.back();
            found = true;
            break;
          } else if (samePoint(cp.back(), p) &&
                     (result.empty() || result.back().first != i)) {
            result.push_back({ i, j });
            p = cp.front();
            found = true;
            break;
          }
        }
      }
      if (found)
        break;
    }
    assert(found);
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

static BSSurface ribbonToSurface(const std::array<BSCurve, 2> &ribbon) {
  const auto &basis = ribbon[0].basis();
  PointVector cpts;
  const auto &cp0 = ribbon[0].controlPoints(), &cp1 = ribbon[1].controlPoints();
  for (size_t i = 0; i < cp0.size(); ++i) {
    cpts.push_back(cp0[i]);
    cpts.push_back(cp0[i] + (cp1[i] - cp0[i]) * 3);
  }
  return { basis.degree(), 1, basis.knots(), { 0, 0, 1, 1 }, cpts };
}

// Returns { S, Su, Suu, Sv, Suv }
// - `u` is the side direction away from the specified corner
// - `v` is the cross-direction, into the surface interior
static VectorVector extractDerivatives(const BSSurface &surface, size_t side, bool at_end) {
  Point2DVector uvs = { {0,0},{0,1}, {0,0},{1,0}, {1,0},{1,1}, {0,1},{1,1} };
  auto uv = uvs[2 * side + (at_end ? 1 : 0)];
  bool swap_uv = side == 0 || side == 2;
  double u_sign = at_end ? -1 : 1;
  double v_sign = side == 2 || side == 3 ? -1 : 1;
  VectorMatrix der;
  surface.eval(uv[0], uv[1], 2, der);
  if (swap_uv)
    return { der[0][0], der[0][1] * u_sign, der[0][2],
             der[1][0] * v_sign, der[1][1] * u_sign * v_sign };
  return { der[0][0], der[1][0] * u_sign, der[2][0],
           der[0][1] * v_sign, der[1][1] * u_sign * v_sign };
}

std::vector<BSSurface> QuadFit::fit() {
  // Topology & geometry structures
  std::vector<Point3D> vertices;
  std::vector<std::set<size_t>> quad_indices;       // vertex -> quads
  std::vector<std::pair<size_t, size_t>> endpoints; // segment -> from_vertex, to_vertex
  std::vector<std::vector<std::pair<size_t, size_t>>> adjacent; // segment -> (quad,side)*

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
  adjacent.resize(vertices.size());
  for (size_t i = 0; i < quads.size(); ++i)
    for (size_t j = 0; j < 4; ++j) {
      size_t s = quads[i].boundaries[j].segment;
      quad_indices[endpoints[s].first].insert(i);
      quad_indices[endpoints[s].second].insert(i);
      adjacent[s].emplace_back(i, j);
    }


  // 0. Fit vertex curvatures
  PointVector all_points;
  for (const auto &q : quads)
    all_points.insert(all_points.end(), q.samples.begin(), q.samples.end());
  auto jet = JetWrapper::fit(vertices, JetWrapper::Nearest(all_points));
  // Note: only those at inner vertices should be used

  // writeVertexCurvatures(vertices, jet, "/tmp/curvatures.vtk");


  // 1. Simple C0 fit
  // - cubic Bezier surfaces
  // - boundaries by the boundary curve endpoints & derivatives
  // - twists by parallelogram rule
  auto result = initialFit();


  // 2. Fit ribbons
  for (size_t i = 0; i < ribbons.size(); ++i) {
    auto rs = ribbonSegments(i);
    Point2DVector sh;
    for (auto s : rs) {
      const auto &b = quads[s.first].boundaries[s.second];
      if (b.s0 < b.s1) {
        if (sh.empty())
          sh.emplace_back(b.s0, b.h0);
        sh.emplace_back(b.s1, b.h1);
      } else {
        if (sh.empty())
          sh.emplace_back(b.s1, b.h1);
        sh.emplace_back(b.s0, b.h0);
      }
    }
    auto slices = fitSlices(ribbonToSurface(ribbons[i]), sh);
    for (size_t j = 0; j < slices.size(); ++j) {
      auto s = rs[j];
      auto &b = quads[s.first].boundaries[s.second];
      b.sextic = slices[j];
      if (b.reversed)
        b.sextic.reverseU();
      b.sextic.normalize();
    }
  }


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
      if (b.on_ribbon) {
        VectorMatrix der;
        const auto &basis = b.sextic.basisU();
        double interval = basis.high() - basis.low();
        b.sextic.eval(basis.low(), 0, 1, der);
        cpts[tangent_cps[2*side]] = cpts[vertex_cps[2*side]] + der[1][0] / 3 * interval;
        b.sextic.eval(basis.high(), 0, 1, der);
        cpts[tangent_cps[2*side+1]] = cpts[vertex_cps[2*side+1]] - der[1][0] / 3 * interval;
      } else {
        auto projectToPlane = [&](size_t v, size_t t, const Vector3D &n) {
          cpts[t] += n * (n * (cpts[v] - cpts[t]));
        };
        const auto &n1 = jet[endpoints[b.segment].first].normal;
        const auto &n2 = jet[endpoints[b.segment].second].normal;
        // Replace start/end tangent based on the ribbons where applicable
        VectorMatrix der;
        if (bp.on_ribbon) {
          if (side == 2 || side == 3)
            bp.sextic.eval(bp.sextic.basisU().high(), 0, 1, der);
          else
            bp.sextic.eval(bp.sextic.basisU().low(), 0, 1, der);
        } else if (bn.on_ribbon) {
          if (side == 2 || side == 3)
            bn.sextic.eval(bn.sextic.basisU().high(), 0, 1, der);
          else
            bn.sextic.eval(bn.sextic.basisU().low(), 0, 1, der);
        }
        bool replace_start = ((bp.on_ribbon && (side == 1 || side == 2)) ||
                              (bn.on_ribbon && (side == 0 || side == 3)));
        bool replace_end   = ((bp.on_ribbon && (side == 0 || side == 3)) ||
                              (bn.on_ribbon && (side == 1 || side == 2)));
        if (replace_start)
          cpts[tangent_cps[2*side]] = cpts[vertex_cps[2*side]] + der[0][1] / 3;
        else
          projectToPlane(vertex_cps[2*side], tangent_cps[2*side], b.reversed ? n2 : n1);
        if (replace_end)
          cpts[tangent_cps[2*side+1]] = cpts[vertex_cps[2*side+1]] + der[0][1] / 3;
        else
          projectToPlane(vertex_cps[2*side+1], tangent_cps[2*side+1], b.reversed ? n1 : n2);
      }
    }
  }


  // 4. Compute better second derivatives for the quad boundary curves

  // First degree elevate to quintic, just for storing the second derivatives
  for (auto &s : result)
    s = elevateBezier(s, 5);

  constexpr std::array<size_t, 8> vertex_cps5 = { 0, 5, 0, 30, 30, 35, 5, 35 };
  constexpr std::array<size_t, 8> tangent_cps5 = { 1, 4, 6, 24, 31, 34, 11, 29 };

  for (size_t i = 0; i < quads.size(); ++i) {
    const auto &q = quads[i];
    auto &cpts = result[i].controlPoints();
    for (size_t side = 0; side < 4; ++side) {
      const auto &b = q.boundaries[side];
      if (b.on_ribbon)          // we already have everything we need here
        continue;
      auto fixSecondDerivative = [&](bool at_end) {
        size_t j = at_end ? 2 * side + 1 : 2 * side;
        const auto &ends = endpoints[b.segment];
        bool seg_end = b.reversed ? !at_end : at_end;
        size_t v = seg_end ? ends.second : ends.first;
        const auto &n = jet[v].normal;
        auto der = (cpts[tangent_cps5[j]] - cpts[vertex_cps5[j]]) * 5 * (at_end ? -1 : 1);
        size_t index = vertex_cps5[j] + ((int)tangent_cps5[j] - (int)vertex_cps5[j]) * 2;
        double cos_theta = jet[v].d_max * der.normalized();
        double c2 = cos_theta * cos_theta, s2 = 1 - c2;
        double k = jet[v].k_max * c2 + jet[v].k_min * s2;
        double h = k * der.normSqr();
        cpts[index] += n * (n * (cpts[vertex_cps5[j]] - cpts[index]) + h / 20);
      };
      // TODO: we may need better 2nd derivatives at the outer vertices, as well
      if (!q.boundaries[(side+3)%4].on_ribbon)
        fixSecondDerivative(side == 0 || side == 3);
      if (!q.boundaries[(side+1)%4].on_ribbon)
        fixSecondDerivative(side == 1 || side == 2);
    }
  }


  // 5. Compute twist vectors
  constexpr std::array<size_t, 24> corner_cps = // S, Su, Sv, Suu, Svv, Suv
    {
       0,  6,  1, 12,  2,  7, // u = 0, v = 0
      30, 24, 31, 18, 32, 25, // u = 1, v = 0
      35, 29, 34, 23, 33, 28, // u = 1, v = 1
       5, 11,  4, 17,  3, 10  // u = 0, v = 1
    };
  for (size_t i = 0; i < quads.size(); ++i) {
    auto &cpts = result[i].controlPoints();
    for (size_t side = 0; side < 4; ++side) {
      const auto &b = quads[i].boundaries[side];
      const auto &bn = quads[i].boundaries[(side+1)%4];

      // If there are 2 ribbons here, we don't need to compute the twist
      if (b.on_ribbon && bn.on_ribbon)
        continue;

      size_t CP = corner_cps[6*side], CPu = corner_cps[6*side+1], CPv = corner_cps[6*side+2],
        CPuu = corner_cps[6*side+3], CPvv = corner_cps[6*side+4], CPuv = corner_cps[6*side+5];

      // If there is 1 fixed ribbon, take the twist from there
      if (b.on_ribbon) {
        VectorMatrix der;
        b.sextic.eval(b.sextic.basisU().high(), 0, 1, der);
        cpts[CPuv] = der[1][1] / 25 + cpts[CPu] + cpts[CPv] - cpts[CP];
      } else if (bn.on_ribbon) {
        VectorMatrix der;
        bn.sextic.eval(bn.sextic.basisU().low(), 0, 1, der);
        cpts[CPuv] = der[1][1] / 25 + cpts[CPu] + cpts[CPv] - cpts[CP];
      } else {
        // Otherwise use the surface curvature to tweak the CP height
        auto Su = (cpts[CPu] - cpts[CP]) * 5 * (side == 1 || side == 2 ? -1 : 1);
        auto Sv = (cpts[CPv] - cpts[CP]) * 5 * (side == 2 || side == 3 ? -1 : -1);
        auto Suu = (cpts[CPuu] - cpts[CPu] * 2 + cpts[CP]) * 20;
        auto Svv = (cpts[CPvv] - cpts[CPv] * 2 + cpts[CP]) * 20;
        bool at_end = side == 1 || side == 2;
        const auto &ends = endpoints[b.segment];
        bool seg_end = b.reversed ? !at_end : at_end;
        size_t v = seg_end ? ends.second : ends.first;
        const auto &n = jet[v].normal;
        double K = jet[v].k_max * jet[v].k_min;
        double h = std::sqrt((Suu * n) * (Svv * n) - K * (Su ^ Sv).normSqr());
        double sign = (side == 1 || side == 3) ? -1 : 1;
        cpts[CPuv] += n * (n * (cpts[CP] - cpts[CPuv]) + h * sign / 25);
      }
    }
  }


  // 6. Compute inner boundary ribbons
  for (size_t i = 0; i < quads.size(); ++i) {
    for (size_t side = 0; side < 4; ++side) {
      auto &b = quads[i].boundaries[side];
      if (b.on_ribbon)
        continue;
      const auto &adj = adjacent[b.segment];
      const auto &opp = adj[0].first == i ? adj[1] : adj[0];
      const auto &b_opp = quads[opp.first].boundaries[opp.second];
      
      auto der1_0 = extractDerivatives(result[i], side, false);
      auto der1_1 = extractDerivatives(result[i], side, true);
      auto der2_0 = extractDerivatives(result[opp.first], opp.second, false);
      auto der2_1 = extractDerivatives(result[opp.first], opp.second, true);
      BSCurve c(4, { 0, 0, 0, 0, 0, 0.5, 1, 1, 1, 1, 1 },
                {  der1_0[0],
                   der1_0[0] + der1_0[1] / 8,
                   der1_0[0] + der1_0[1] * 3 / 8 + der1_0[2] / 24,
                   der1_1[0] + der1_1[1] * 3 / 8 + der1_1[2] / 24,
                   der1_1[0] + der1_1[1] / 8,
                   der1_1[0] });
      BSCurve c1({
          der1_0[0] + der1_0[3],
          der1_0[0] + der1_0[1] / 3 + der1_0[3] + der1_0[4] / 3,
          der1_1[0] + der1_1[1] / 3 + der1_1[3] + der1_1[4] / 3,
          der1_1[0] + der1_1[3]});
      c1.insertKnot(0.5, 1);
      BSCurve c2({
          der2_0[0] + der2_0[3],
          der2_0[0] + der2_0[1] / 3 + der2_0[3] + der2_0[4] / 3,
          der2_1[0] + der2_1[1] / 3 + der2_1[3] + der2_1[4] / 3,
          der2_1[0] + der2_1[3]});
      if (b.reversed != b_opp.reversed)
        c2.reverse();
      c2.insertKnot(0.5, 1);
      b.sextic = connectG1(c, c1, c2);
      // writeSTL({b.sextic}, std::string("/tmp/sextic") + std::to_string(i) + "_" + std::to_string(side) + ".stl", 50);
    }
  }


  // 7. Fill sextic patches
  for (size_t i = 0; i < quads.size(); ++i) {
    auto &quad = quads[i];
    result[i] = elevateBezier(result[i], 6);
    {
      auto &s1 = quad.boundaries[0].sextic;
      auto &s2 = quad.boundaries[2].sextic;
      auto p1 = s1.basisU().knots().begin();
      auto p2 = s2.basisU().knots().begin();
      while (*p1 < 1 || *p2 < 1) {
        if (*p1 < *p2) {
          result[i].insertKnotV(*p1, 1);
          s2.insertKnotU(*p1, 1);
          p1++;
        } else if (*p2 < *p1) {
          result[i].insertKnotV(*p2, 1);
          s1.insertKnotU(*p2, 1);
          p2++;
        } else {
          result[i].insertKnotV(*p1, 1);
          p1++;
          p2++;
        }
      }
    }
    {
      auto &s1 = quad.boundaries[1].sextic;
      auto &s2 = quad.boundaries[3].sextic;
      auto p1 = s1.basisU().knots().begin();
      auto p2 = s2.basisU().knots().begin();
      while (*p1 < 1 || *p2 < 1) {
        if (*p1 < *p2) {
          result[i].insertKnotU(*p1, 1);
          s2.insertKnotU(*p1, 1);
          p1++;
        } else if (*p2 < *p1) {
          result[i].insertKnotU(*p2, 1);
          s1.insertKnotU(*p2, 1);
          p2++;
        } else {
          result[i].insertKnotU(*p1, 1);
          p1++;
          p2++;
        }
      }
    }
    auto n_cpts = result[i].numControlPoints();
    auto &s = result[i];
    const auto &r0 = quad.boundaries[0].sextic;
    const auto &r1 = quad.boundaries[1].sextic;
    const auto &r2 = quad.boundaries[2].sextic;
    const auto &r3 = quad.boundaries[3].sextic;
    for (size_t j = 0; j < n_cpts[1]; ++j) {
      s.controlPoint(0, j) = r0.controlPoint(j, 0);
      s.controlPoint(1, j) =
        r0.controlPoint(j, 0) + (r0.controlPoint(j, 1) - r0.controlPoint(j, 0)) / 6;
      s.controlPoint(n_cpts[0] - 1, j) = r2.controlPoint(j, 0);
      s.controlPoint(n_cpts[0] - 2, j) =
        r2.controlPoint(j, 0) + (r2.controlPoint(j, 1) - r2.controlPoint(j, 0)) / 6;
    }
    for (size_t j = 0; j < n_cpts[0]; ++j) {
      s.controlPoint(j, 0) = r1.controlPoint(j, 0);
      s.controlPoint(j, 1) =
        r1.controlPoint(j, 0) + (r1.controlPoint(j, 1) - r1.controlPoint(j, 0)) / 6;
      s.controlPoint(j, n_cpts[1] - 1) = r3.controlPoint(j, 0);
      s.controlPoint(j, n_cpts[1] - 2) =
        r3.controlPoint(j, 0) + (r3.controlPoint(j, 1) - r3.controlPoint(j, 0)) / 6;
    }
  }


  // 8. Fit sampled points using inner control points
  // for (size_t i = 0; i < quads.size(); ++i) {
  //   auto ncp = result[i].numControlPoints();
  //   auto fix2 = [&](size_t i, size_t j) {
  //     return i < 2 || j < 2 || i >= ncp[0] - 2 || j >= ncp[1] - 2;
  //   };
  //   bsplineFit(result[i], quads[i].resolution, quads[i].samples, fix2, 0.1);
  // }

  return result;
}
