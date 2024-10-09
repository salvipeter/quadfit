#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>

#include "quadfit.hh"

#include "bspline-fit.hh"
#include "connect-g1.hh"
#include "discrete-mask.hh"
#include "fit-ribbon.hh"
#include "io.hh"
#include "switches.hh"

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

  // Part 4 : Sampled points/normals
  for (size_t i = 0; i < n_quads; ++i) {
    size_t res;
    f >> res;
    quads[i].resolution = res;
    size_t n_samples = (res + 1) * (res + 1);
    quads[i].samples.resize(n_samples);
    quads[i].normals.resize(n_samples);
    for (size_t j = 0; j < n_samples; ++j) {
      quads[i].samples[j] = readPoint(f);
      quads[i].normals[j] = readPoint(f);
      quads[i].normals[j].normalize();
    }
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

void QuadFit::update(const std::vector<std::string> &switches) {
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

  adjacency.resize(segments.size());
  for (size_t i = 0; i < quads.size(); ++i)
    for (size_t j = 0; j < 4; ++j) {
      size_t s = quads[i].boundaries[j].segment;
      adjacency[s].emplace_back(i, j);
    }

  PointVector all_points;
  std::string mesh_filename;
  parseSwitch<std::string>(switches, "mesh", &mesh_filename);
  if (mesh_filename.empty())
    for (const auto &q : quads)
      all_points.insert(all_points.end(), q.samples.begin(), q.samples.end());
  else
    all_points = TriMesh::readOBJ(mesh_filename).points();

  jet = JetWrapper::fit(vertices, JetWrapper::Nearest(all_points));
  // writeVertexCurvatures(vertices, jet, "/tmp/curvatures.vtk");
}

BSSurface QuadFit::preliminaryFit(size_t i) const {
  BSSurface result(3, 3, PointVector(16));
  for (size_t j = 1; j <= 15; ++j) {
    result = result.insertKnotU(j / 16.0, 1);
    result = result.insertKnotV(j / 16.0, 1);
  }
  auto constraint = [&](size_t i, size_t j) -> MoveConstraint {
    return MoveType::Free();
  };
  bsplineFit(result, quads[i].resolution, quads[i].samples, constraint, 0, false);
  return result;
}

// Fit cubic Bezier surfaces
// - boundaries by the boundary curve endpoints & derivatives
// - twists by parallelogram rule
std::vector<BSSurface> QuadFit::initialFit(bool fit_tangents, bool fit_twists) const {
  std::vector<BSSurface> result(quads.size());

  for (size_t i = 0; i < quads.size(); ++i) {
    auto &quad = quads[i];
    auto &b = quad.boundaries;

    PointVector cpts(16);
    for (size_t j = 0; j < 4; ++j) {
      auto getCP = [&](size_t k) {
        const auto &curve = segments[b[k].segment];
        bool rev = (b[k].reversed && j < 2) || (!b[k].reversed && j >= 2);
        bool deriv = j == 1 || j == 2;
        double u = rev ? curve.basis().high() : curve.basis().low();
        VectorVector der;
        auto p = curve.eval(u, deriv ? 1 : 0, der);
        if (deriv)
          return p + der[1] * (rev ? -1 : 1);
        return p;
      };
      cpts[j] = getCP(0);
      cpts[j*4] = getCP(1);
      cpts[12+j] = getCP(2);
      cpts[3+j*4] = getCP(3);
    }

    if ((fit_tangents || fit_twists) && !quad.preliminary_fit)
      quad.preliminary_fit = preliminaryFit(i);

    if (fit_tangents) {
      auto setTangentLength = [&](bool u, bool u_end, bool v_end) {
        size_t i = u_end ? (v_end ? 15 : 12) : (v_end ? 3 : 0);
        size_t j = i + (u ? (u_end ? -4 : 4) : (v_end ? -1 : 1));
        VectorMatrix der;
        quad.preliminary_fit->eval(u_end ? 1 : 0, v_end ? 1 : 0, 1, der);
        Vector3D d = u ? der[1][0] : der[0][1];
        if ((u && u_end) || (!u && v_end))
          d *= -1;
        cpts[j] = cpts[i] + d / 3;
      };
      setTangentLength(false, false, false);
      setTangentLength(false, false, true);
      setTangentLength(false, true, false);
      setTangentLength(false, true, true);
      setTangentLength(true, false, false);
      setTangentLength(true, false, true);
      setTangentLength(true, true, false);
      setTangentLength(true, true, true);
    } else {
      auto setTangentLength = [&](size_t p1, size_t t1, size_t t2, size_t p2) {
        double d = (cpts[p2] - cpts[p1]).norm() / 3;
        cpts[t1] = cpts[p1] + (cpts[t1] - cpts[p1]).normalize() * d;
        cpts[t2] = cpts[p2] + (cpts[t2] - cpts[p2]).normalize() * d;
      };
      setTangentLength(0, 1, 2, 3);
      setTangentLength(0, 4, 8, 12);
      setTangentLength(12, 13, 14, 15);
      setTangentLength(3, 7, 11, 15);
    }

    if (fit_twists) {
      auto setTwist = [&](bool u_end, bool v_end) {
        size_t CP = u_end ? (v_end ? 15 : 12) : (v_end ? 3 : 0);
        size_t CPu = CP + (u_end ? -4 : 4);
        size_t CPv = CP + (v_end ? -1 : 1);
        size_t CPuv = CPu + (v_end ? -1 : 1);
        VectorMatrix der;
        quad.preliminary_fit->eval(u_end ? 1 : 0, v_end ? 1 : 0, 2, der);
        bool at_end = (u_end && !v_end) || (!u_end && v_end);
        cpts[CPuv] = der[1][1] * (at_end ? -1 : 1) / 9 + cpts[CPu] + cpts[CPv] - cpts[CP];
      };
      setTwist(false, false);
      setTwist(false, true);
      setTwist(true, false);
      setTwist(true, true);
    } else {
      auto setTwist = [&](size_t p, size_t d1, size_t d2, size_t t) {
        cpts[t] = cpts[d1] + (cpts[d2] - cpts[p]);
      };
      setTwist(0, 1, 4, 5);
      setTwist(3, 2, 7, 6);
      setTwist(12, 8, 13, 9);
      setTwist(15, 11, 14, 10);
    }

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

void QuadFit::correctFirstDerivatives(BSSurface &cubic, size_t quad_index) const {
  static constexpr std::array<size_t, 8> vertex_cps = { 0, 3, 0, 12, 12, 15, 3, 15 };
  static constexpr std::array<size_t, 8> tangent_cps = { 1, 2, 4, 8, 13, 14, 7, 11 };

  auto &cpts = cubic.controlPoints();
  for (size_t side = 0; side < 4; ++side) {
    const auto &quad = quads[quad_index];
    const auto &b = quad.boundaries[side];
    const auto &bp = quad.boundaries[(side+3)%4];
    const auto &bn = quad.boundaries[(side+1)%4];
    if (b.on_ribbon) {
      // Compute first derivative from the fitted ribbon
      VectorMatrix der;
      const auto &knots = b.sextic.basisU().knots();
      double lo = knots[7] - knots[1], hi = knots.rbegin()[1] - knots.rbegin()[7];
      b.sextic.eval(knots[6], 0, 1, der);
      cpts[tangent_cps[2*side]] = cpts[vertex_cps[2*side]] + der[1][0] / 3 * lo;
      b.sextic.eval(knots.rbegin()[6], 0, 1, der);
      cpts[tangent_cps[2*side+1]] = cpts[vertex_cps[2*side+1]] - der[1][0] / 3 * hi;
    } else {
      // In the interior project to the vertices' normal plane
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

static Point3D intersectLines(const Point3D &ap, const Vector3D &ad,
                              const Point3D &bp, const Vector3D &bd) {
  double a = ad * ad, b = ad * bd, c = bd * bd;
  double d = ad * (ap - bp), e = bd * (ap - bp);
  if (a * c - b * b < 1.0e-7)
    return ap;
  double D = a * c - b * b;
  double s = (b * e - c * d) / D;
  double t = (a * e - b * d) / D;
  return ((ap + ad * s) + (bp + bd * t)) / 2;
}

// At the corner between side `side` and `side+1`
void QuadFit::correctCubicTwists(BSSurface &cubic, size_t quad_index, bool heuristic) const {
  static constexpr std::array<size_t, 16> corner_cps = // S, Su, Sv, Suv
    {
       0,  4,  1,  5, // u = 0, v = 0
      12,  8, 13,  9, // u = 1, v = 0
      15, 11, 14, 10, // u = 1, v = 1
       3,  7,  2,  6  // u = 0, v = 1
    };
  auto &cpts = cubic.controlPoints();

  std::vector<size_t> not_fixed;

  for (size_t side = 0; side < 4; ++side) {
    const auto &b = quads[quad_index].boundaries[side];
    const auto &bn = quads[quad_index].boundaries[(side+1)%4];

    size_t CP = corner_cps[4*side], CPu = corner_cps[4*side+1],
      CPv = corner_cps[4*side+2], CPuv = corner_cps[4*side+3];

    // If there is a fixed ribbon, take the twist from there
    // (assuming that 2 fixed ribbons share a common twist);
    // otherwise we will need the 2nd derivatives, see correctTwists()
    auto fixByRibbon = [&](const Quad::Boundary &b, bool at_end) {
      double u = at_end ? b.sextic.basisU().high() : b.sextic.basisU().low();
      VectorMatrix der;
      b.sextic.eval(u, 0, 2, der);
      cpts[CPuv] = der[1][1] * (at_end ? -1 : 1) / 9 + cpts[CPu] + cpts[CPv] - cpts[CP];
    };
    if (b.on_ribbon)
      fixByRibbon(b, side == 1 || side == 2);
    else if (bn.on_ribbon)
      fixByRibbon(bn, side == 2 || side == 3);
    else
      not_fixed.push_back(CPuv);
  }

  if (!heuristic)
    return;

  // For the time being, set the rest to the mean of the adjacent CPs,
  // or at the quasi-intersection of two fixed segments
  if (not_fixed.size() == 1) {
    size_t i = not_fixed[0];
    cpts[i] = (cpts[i-1] + cpts[i+1] + cpts[i-4] + cpts[i+4]) / 4;
    // cpts[i] = intersectLines(cpts[i-1], cpts[i+1] - cpts[i-1], cpts[i-4], cpts[i+4] - cpts[i-4]);
  } else {
    size_t i = not_fixed[0], j = not_fixed[1];
    if (i > j)
      std::swap(i, j);
    size_t a = i - (j - i), b = j + (j - i);
    if (j - i == 1) {
      cpts[i] = intersectLines(cpts[a], cpts[b] - cpts[a], cpts[i-4], cpts[i+4] - cpts[i-4]);
      cpts[j] = intersectLines(cpts[a], cpts[b] - cpts[a], cpts[j-4], cpts[j+4] - cpts[j-4]);
    } else {
      cpts[i] = intersectLines(cpts[a], cpts[b] - cpts[a], cpts[i-1], cpts[i+1] - cpts[i-1]);
      cpts[j] = intersectLines(cpts[a], cpts[b] - cpts[a], cpts[j-1], cpts[j+1] - cpts[j-1]);
    }
  }
}

void QuadFit::correctSecondDerivatives(BSSurface &quintic, size_t quad_index) const {
  static constexpr std::array<size_t, 8> vertex_cps = { 0, 5, 0, 30, 30, 35, 5, 35 };
  static constexpr std::array<size_t, 8> tangent_cps = { 1, 4, 6, 24, 31, 34, 11, 29 };

  const auto &q = quads[quad_index];
  auto &cpts = quintic.controlPoints();
  for (size_t side = 0; side < 4; ++side) {
    const auto &b = q.boundaries[side];
    if (b.on_ribbon)          // we already have everything we need here
      continue;
    // Project the second derivative into the plane at the constrained height
    auto fix = [&](bool at_end) {
      size_t j = at_end ? 2 * side + 1 : 2 * side;
      const auto &ends = endpoints[b.segment];
      bool seg_end = b.reversed ? !at_end : at_end;
      size_t v = seg_end ? ends.second : ends.first;
      const auto &n = jet[v].normal;
      auto der = (cpts[tangent_cps[j]] - cpts[vertex_cps[j]]) * 5;
      size_t index = vertex_cps[j] + ((int)tangent_cps[j] - (int)vertex_cps[j]) * 2;
      double cos_theta = jet[v].d_max * der.normalized();
      double c2 = cos_theta * cos_theta, s2 = 1 - c2;
      double k = jet[v].k_max * c2 + jet[v].k_min * s2;
      double h = k * der.normSqr();
      cpts[index] += n * (n * (cpts[vertex_cps[j]] - cpts[index]) + h / 20);
    };
    // TODO: we may need better 2nd derivatives at the outer vertices, as well
    if (!q.boundaries[(side+3)%4].on_ribbon)
      fix(side == 0 || side == 3);
    if (!q.boundaries[(side+1)%4].on_ribbon)
      fix(side == 1 || side == 2);
  }
}

// At the corner between side `side` and `side+1`
void QuadFit::correctTwists(BSSurface &quintic, size_t quad_index) const {
  static constexpr std::array<size_t, 8> corner_cps = { 0, 7, 30, 25, 35, 28, 5, 10 };
  static constexpr std::array<double, 8> uv = { 0,0, 1,0, 1,1, 0,1 };
  auto &cpts = quintic.controlPoints();

  for (size_t side = 0; side < 4; ++side) {
    const auto &b = quads[quad_index].boundaries[side];
    const auto &bn = quads[quad_index].boundaries[(side+1)%4];

    // If it is near a ribbon, the twist is already computed
    if (b.on_ribbon || bn.on_ribbon)
      continue;

    // Otherwise use the surface curvature to tweak the CP height
    size_t CP = corner_cps[2*side], CPuv = corner_cps[2*side+1];
    VectorMatrix der;
    quintic.eval(uv[2*side], uv[2*side+1], 2, der);
    const auto &Su = der[1][0];
    const auto &Sv = der[0][1];
    const auto &Suu = der[2][0];
    const auto &Svv = der[0][2];
    bool at_end = side == 1 || side == 2;
    const auto &ends = endpoints[b.segment];
    bool seg_end = b.reversed ? !at_end : at_end;
    size_t v = seg_end ? ends.second : ends.first;
    const auto &n = jet[v].normal;
    double K = jet[v].k_max * jet[v].k_min;
    double h = std::sqrt((Suu * n) * (Svv * n) - K * (Su ^ Sv).normSqr());
    double sign = (side == 1 || side == 3) ? -1 : 1;
    {
      auto E = Su.normSqr(), F = Su * Sv, G = Sv.normSqr();
      auto L = Suu * n, M = h, N = Svv * n;
      auto H = jet[v].k_min + jet[v].k_max;
      auto H1 = (L * G - 2 * M * F + N * E) / (E * G - F * F);
      if (std::abs(H - H1) > epsilon) // dirty hack
        sign *= -1;
    }
    cpts[CPuv] += n * (n * (cpts[CP] - cpts[CPuv]) + h * sign / 25);
  }
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

static auto edgeSamples(const PointVector &samples, const VectorVector &normals,
                         size_t side, size_t res) {
  std::vector<std::pair<Point3D, Vector3D>> result;
  for (size_t i = 0; i <= res; ++i) {
    size_t index = i;
    if (side == 1)
      index = index * (res + 1);
    else if (side == 2)
      index = res * (res + 1) + index;
    else if (side == 3)
      index = index * (res + 1) + res;
    const auto &p = samples[index];
    const auto &n = normals[index];
    result.emplace_back(p, n);
  }
  return result;
}

[[maybe_unused]]
static Point3D closestPoint(const BSCurve &curve, const Point3D &point, double &u,
                            size_t max_iteration, double distance_tol, double cosine_tol) {
  VectorVector der;
  auto lo = curve.basis().low();
  auto hi = curve.basis().high();

  for (size_t iteration = 0; iteration < max_iteration; ++iteration) {
    auto deviation = curve.eval(u, 2, der) - point;
    auto distance = deviation.norm();
    if (distance < distance_tol)
      break;

    double scaled_error = der[1] * deviation;
    double cosine_err = std::abs(scaled_error) / (der[1].norm() * distance);
    if (cosine_err < cosine_tol)
      break;

    double old = u;
    u -= scaled_error / (der[2] * deviation + der[1] * der[1]);
    u = std::min(std::max(u, lo), hi);

    if ((der[1] * (u - old)).norm() < distance_tol)
      break;
  }

  return der[0];
}

static Vector3D quadNormal(const BSSurface &surface, size_t side, double t) {
  double u = 0, v = t;
  if (side == 1)
    std::swap(u, v);
  else if (side == 2)
    u = 1;
  else if (side == 3) {
    u = v;
    v = 1;
  }
  VectorMatrix der;
  surface.eval(u, v, 1, der);
  return (der[1][0] ^ der[0][1]).normalize();
}

BSSurface QuadFit::innerBoundaryRibbon(const std::vector<BSSurface> &quintic_patches,
                                       size_t quad_index, size_t side, size_t extra_knots,
                                       bool prelim_normals, bool fitC0, bool fitG1) const {
  const auto &result = quintic_patches;
  const auto &quad = quads[quad_index];
  const auto &b = quad.boundaries[side];
  const auto &adj = adjacency[b.segment];
  const auto &opp = adj[0].first == quad_index ? adj[1] : adj[0];
  const auto &b_opp = quads[opp.first].boundaries[opp.second];

  auto der1_0 = extractDerivatives(result[quad_index], side, false);
  auto der1_1 = extractDerivatives(result[quad_index], side, true);
  auto der2_0 = extractDerivatives(result[opp.first], opp.second, false);
  auto der2_1 = extractDerivatives(result[opp.first], opp.second, true);
  BSCurve c(4, { 0, 0, 0, 0, 0, 0.5, 1, 1, 1, 1, 1 },
            {  der1_0[0],
               der1_0[0] + der1_0[1] / 8,
               der1_0[0] + der1_0[1] * 3 / 8 + der1_0[2] / 24,
               der1_1[0] + der1_1[1] * 3 / 8 + der1_1[2] / 24,
               der1_1[0] + der1_1[1] / 8,
               der1_1[0] });
  for (size_t i = 1; i <= extra_knots; ++i) {
    double u = (double)i / (extra_knots + 1);
    c = c.insertKnot(u / 2, 1);
    c = c.insertKnot(0.5 + u / 2, 1);
  }

  if (prelim_normals) {
    if (!quad.preliminary_fit)
      quad.preliminary_fit = preliminaryFit(quad_index);
    if (!quads[opp.first].preliminary_fit)
      quads[opp.first].preliminary_fit = preliminaryFit(opp.first);
  }

  // Prepare sampled points & normals
  size_t res = quad.resolution;
  PointVector points;
  VectorVector normals;
  if (fitC0 || fitG1) {
    auto cross = edgeSamples(quad.samples, quad.normals, side, res);
    for (size_t i = 0; i <= res; ++i) {
      const auto &[p, normal] = cross[i];
      points.push_back(p);
      if (prelim_normals) {
        double t = (double)i / res;
        auto n1 = quadNormal(*quad.preliminary_fit, side, t);
        if (b.reversed != b_opp.reversed)
          t = 1 - t;
        auto n2 = quadNormal(*quads[opp.first].preliminary_fit, opp.second, t);
        normals.push_back((n1 + n2).normalize());
      } else
        normals.push_back(normal);
    }
  }

  if (fitC0) {
    // Better curve approximation
    // Only fix the endpoints and the tangent near the ribbon
    bool start_on_ribbon, end_on_ribbon;
    if (side == 0 || side == 3) {
      start_on_ribbon = quad.boundaries[(side+1)%4].on_ribbon;
      end_on_ribbon = quad.boundaries[(side+3)%4].on_ribbon;
    } else {
      start_on_ribbon = quad.boundaries[(side+3)%4].on_ribbon;
      end_on_ribbon = quad.boundaries[(side+1)%4].on_ribbon;
    }
    auto constraint = [=](size_t i) -> MoveConstraint {
      size_t max = c.controlPoints().size() - 1;
      if ((i == 0 || (start_on_ribbon && i == 1)) ||
          (i == max || (end_on_ribbon && i == max - 1)))
        return MoveType::Fixed();
      return MoveType::Free();
      // if (i == 2)
      //   return MoveType::Tangent({(der1_0[1] ^ der1_0[3]).normalize()});
      // if (i > 2 && i < max - 2)
      //   return MoveType::Free();
      // if (i == max - 2)
      //   return MoveType::Tangent({(der1_1[1] ^ der1_1[3]).normalize()});
      // return MoveType::Fixed();
    };
    bsplineFit(c, points, /*normals,*/ constraint, 0, 0);
  }

  BSCurve c1({
      der1_0[0] + der1_0[3],
      der1_0[0] + der1_0[1] / 3 + der1_0[3] + der1_0[4] / 3,
      der1_1[0] + der1_1[1] / 3 + der1_1[3] + der1_1[4] / 3,
      der1_1[0] + der1_1[3]});
  BSCurve c2({
      der2_0[0] + der2_0[3],
      der2_0[0] + der2_0[1] / 3 + der2_0[3] + der2_0[4] / 3,
      der2_1[0] + der2_1[1] / 3 + der2_1[3] + der2_1[4] / 3,
      der2_1[0] + der2_1[3]});
  if (b.reversed != b_opp.reversed)
    c2.reverse();
  c1 = c1.insertKnot(0.5, 1);
  c2 = c2.insertKnot(0.5, 1);
  for (size_t i = 1; i <= extra_knots; ++i) {
    double u = (double)i / (extra_knots + 1);
    c1 = c1.insertKnot(u / 2, 1);
    c1 = c1.insertKnot(0.5 + u / 2, 1);
    c2 = c2.insertKnot(u / 2, 1);
    c2 = c2.insertKnot(0.5 + u / 2, 1);
  }

  if (fitG1) {
    // Better normal approximation
    points.clear();
    for (size_t i = 0; i <= res; ++i)
      points.push_back(c.eval((double)i / res));

    auto fixCenter = [&](BSCurve &curve) {
      auto constraint = [=](size_t i) -> MoveConstraint {
        if (i >= 2 && i < curve.controlPoints().size() - 2)
          return MoveType::Normal({normals[res/2]}); // Note: this assumes that resolution is even!
        return MoveType::Fixed();
      };
      bsplineFit(curve, points, normals, constraint, 0);
    };

    fixCenter(c1);
    fixCenter(c2);
  }

  return connectG1(c, c1, c2);
}

static void unifyKnots(BSSurface &sextic, BSSurface &s1, BSSurface &s2, bool u_dir) {
  auto p1 = s1.basisU().knots().begin();
  auto p2 = s2.basisU().knots().begin();
  auto t1 = s1, t2 = s2;
  while (*p1 < 1 || *p2 < 1) {
    double to_insert;
    if (*p1 < *p2) {
      to_insert = *p1;
      t2 = t2.insertKnotU(*p1, 1);
      p1++;
    } else if (*p2 < *p1) {
      to_insert = *p2;
      t1 = t1.insertKnotU(*p2, 1);
      p2++;
    } else {
      to_insert = *p1;
      p1++;
      p2++;
    }
    if (u_dir)
      sextic = sextic.insertKnotU(to_insert, 1);
    else
      sextic = sextic.insertKnotV(to_insert, 1);
  }
  s1 = t1;
  s2 = t2;
}

static void fillSextic(BSSurface &s, const BSSurface &r0, const BSSurface &r1,
                       const BSSurface &r2, const BSSurface &r3) {
  auto n_cpts = s.numControlPoints();
  const auto &knots_u = s.basisU().knots();
  double lo = knots_u[7] - knots_u[1], hi = knots_u.rbegin()[1] - knots_u.rbegin()[7];
  for (size_t j = 0; j < n_cpts[1]; ++j) {
    s.controlPoint(0, j) = r0.controlPoint(j, 0);
    s.controlPoint(n_cpts[0] - 1, j) = r2.controlPoint(j, 0);
    if (j != 0 && j != n_cpts[1] - 1) {
      s.controlPoint(1, j) =
        r0.controlPoint(j, 0) + (r0.controlPoint(j, 1) - r0.controlPoint(j, 0)) / 6 * lo;
      s.controlPoint(n_cpts[0] - 2, j) =
        r2.controlPoint(j, 0) + (r2.controlPoint(j, 1) - r2.controlPoint(j, 0)) / 6 * hi;
    }
  }
  const auto &knots_v = s.basisV().knots();
  lo = knots_v[7] - knots_v[1]; hi = knots_v.rbegin()[1] - knots_v.rbegin()[7];
  for (size_t j = 0; j < n_cpts[0]; ++j) {
    s.controlPoint(j, 0) = r1.controlPoint(j, 0);
    s.controlPoint(j, n_cpts[1] - 1) = r3.controlPoint(j, 0);
    if (j != 0 && j != n_cpts[0] - 1) {
      s.controlPoint(j, 1) =
        r1.controlPoint(j, 0) + (r1.controlPoint(j, 1) - r1.controlPoint(j, 0)) / 6 * lo;
      s.controlPoint(j, n_cpts[1] - 2) =
        r3.controlPoint(j, 0) + (r3.controlPoint(j, 1) - r3.controlPoint(j, 0)) / 6 * hi;
    }
  }
}

[[maybe_unused]]
static BSCurve baseCurve(const BSSurface &quad, size_t side) {
  PointVector cpts;
  if (side == 0 || side == 2) {
    // v-side
    size_t m = side == 0 ? 0 : quad.numControlPoints()[0] - 1;
    size_t n = quad.numControlPoints()[1];
    cpts.resize(n);
    for (size_t i = 0; i < n; ++i)
      cpts[i] = quad.controlPoint(m, i);
    return { quad.basisV().degree(), quad.basisV().knots(), cpts };
  }

  // u-side
  double m = side == 1 ? 0 : quad.numControlPoints()[1] - 1;
  size_t n = quad.numControlPoints()[0];
  cpts.resize(n);
  for (size_t i = 0; i < n; ++i)
    cpts[i] = quad.controlPoint(i, m);
  return { quad.basisU().degree(), quad.basisU().knots(), cpts };
}

[[maybe_unused]]
static std::pair<Point3D, Vector3D> evalNormal(const BSSurface &s, size_t side, double u) {
  VectorMatrix der;
  if (side == 0 || side == 2) // v-side
    s.eval(side == 0 ? 0 : 1, u, 1, der);
  else // u-side
    s.eval(u, side == 1 ? 0 : 1, 1, der);
  return { der[0][0], (der[1][0] ^ der[0][1]).normalize() };
}

void QuadFit::printContinuityErrors(const std::vector<BSSurface> &result) const {
  std::cout << "\nMaximal errors:\tC0\tG1 (degrees)" << std::endl;
  for (const auto &adj : adjacency) {
    if (adj.size() < 2) {
      size_t side = adj[0].second;
      const auto &b = quads[adj[0].first].boundaries[side];
      const auto &q = result[adj[0].first];
      size_t resolution = 100;
      const auto &r = ribbons[b.ribbon];
      std::cout << "Ribbon #" << b.ribbon + 1 << ":\t";
      double max_p_error = 0, max_t_error = 0;
      for (size_t i = 0; i <= resolution; ++i) {
        double u = (double)i / resolution;
        double v = b.s0 * (1 - u) + b.s0 * u;
        auto [p1, n1] = evalNormal(q, side, u);
        closestPoint(r[0], p1, v, 20, 0, 0);
        VectorVector der;
        Point3D p2 = r[0].eval(v, 1, der);
        Vector3D n2 = (der[1] ^ (r[1].eval(v) - p2)).normalize();
        double p_error = (p1 - p2).norm();
        double t_error = std::acos(std::min(std::max(n1 * n2, -1.0), 1.0)) * 180 / M_PI;
        max_p_error = std::max(max_p_error, p_error);
        max_t_error = std::max(max_t_error, t_error);
      }
      std::cout << max_p_error << " \t" << max_t_error << std::endl;
      continue;
    }
    std::cout << "Quads #" << adj[0].first + 1 << " - #" << adj[1].first + 1 << ":\t";
    const auto &q1 = result[adj[0].first];
    const auto &q2 = result[adj[1].first];
    size_t side1 = adj[0].second, side2 = adj[1].second;
    bool reversed = quads[adj[0].first].boundaries[side1].reversed !=
      quads[adj[1].first].boundaries[side2].reversed;
    double max_p_error = 0, max_t_error = 0;
    size_t resolution = 100;
    for (size_t i = 0; i <= resolution; ++i) {
      double u = (double)i / resolution;
      auto [p1, n1] = evalNormal(q1, side1, u);
      auto [p2, n2] = evalNormal(q2, side2, reversed ? 1 - u : u);
      double p_error = (p1 - p2).norm();
      double t_error = std::acos(std::min(std::max(n1 * n2, -1.0), 1.0)) * 180 / M_PI;
      max_p_error = std::max(max_p_error, p_error);
      max_t_error = std::max(max_t_error, t_error);
    }
    std::cout << max_p_error << " \t" << max_t_error << std::endl;
  }
  std::cout << std::endl;
}

void QuadFit::printApproximationErrors(const std::vector<BSSurface> &result) const {
  std::cout << "Curve maximal deviation:\tC0\tG1 (degrees)" << std::endl;
  for (const auto &adj : adjacency) {
    if (adj.size() < 2)
      continue;
    size_t side = adj[0].second;
    const auto &q = quads[adj[0].first];
    auto res = q.resolution;
    double max_err = 0, max_t_err = 0;
    size_t max_p = -1, max_t = -1;
    auto cross = edgeSamples(q.samples, q.normals, side, res);
    for (size_t i = 0; i <= res; ++i) {
      auto u = (double)i / res;
      const auto &[p1, n1] = cross[i];
      closestPoint(baseCurve(result[adj[0].first], side), p1, u, 20, 0, 0);
      auto [p2, n2] = evalNormal(result[adj[0].first], side, u);
      auto err = (p1 - p2).norm();
      auto t_err = std::acos(std::min(std::abs(n1 * n2), 1.0)) * 180 / M_PI;
      if (err > max_err) {
        max_err = err;
        max_p = i;
      }
      if (t_err > max_t_err) {
        max_t_err = t_err;
        max_t = i;
      }
    }
    std::cout << "Quads #" << adj[0].first << " - #" << adj[1].first
              << " (" << max_p << " / " << max_t << "):\t"
              << max_err << "\t" << max_t_err << std::endl;
  }
  std::cout << std::endl;
}

std::vector<BSSurface> QuadFit::fit(const std::vector<std::string> &switches) {
  // 1. Simple C0 fit
  auto result = initialFit(parseSwitch<bool>(switches, "preliminary-fit-tangents"),
                           parseSwitch<bool>(switches, "preliminary-fit-twists"));

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
  for (size_t i = 0; i < quads.size(); ++i)
    correctFirstDerivatives(result[i], i);

  // 4. Compute the correct twists for corner- and side vertices
  bool heuristic = !parseSwitch<bool>(switches, "preliminary-fit-twists");
  for (size_t i = 0; i < quads.size(); ++i)
    correctCubicTwists(result[i], i, heuristic);

  // 5. Compute better second derivatives for the quad boundary curves

  // First degree elevate to quintic, just for storing the second derivatives
  for (auto &s : result)
    s = elevateBezier(s, 5);

  for (size_t i = 0; i < quads.size(); ++i)
    correctSecondDerivatives(result[i], i);

  // 6. Compute twist vectors
  for (size_t i = 0; i < quads.size(); ++i)
    correctTwists(result[i], i);

  // 7. Compute inner boundary ribbons
  size_t extra = 0;
  parseSwitch<size_t>(switches, "extra-knots", &extra, 2);
  bool prelim_normals = parseSwitch<bool>(switches, "preliminary-fit-normals");
  bool fitC0 = parseSwitch<bool>(switches, "fit-curves");
  bool fitG1 = parseSwitch<bool>(switches, "fit-normals");
  for (size_t i = 0; i < quads.size(); ++i) {
    for (size_t side = 0; side < 4; ++side) {
      auto &b = quads[i].boundaries[side];
      if (b.on_ribbon)
        continue;
      b.sextic = innerBoundaryRibbon(result, i, side, extra, prelim_normals, fitC0, fitG1);
    }
  }

  // 8. Fill sextic patches
  for (size_t i = 0; i < quads.size(); ++i) {
    auto &quad = quads[i];
    result[i] = elevateBezier(result[i], 6);
    unifyKnots(result[i], quad.boundaries[0].sextic, quad.boundaries[2].sextic, false);
    unifyKnots(result[i], quad.boundaries[1].sextic, quad.boundaries[3].sextic, true);
    fillSextic(result[i], quad.boundaries[0].sextic, quad.boundaries[1].sextic,
               quad.boundaries[2].sextic, quad.boundaries[3].sextic);
  }

  // 9a. Use a mask to compute the placement of the inner control points
  if (parseSwitch<bool>(switches, "coons-patch"))
    for (auto &r : result)
      applyMask(r, DiscreteMask::C1_COONS);
  else {
    // 9. Fit sampled points using inner control points
    size_t inner_knots = 5;
    bool cubic = parseSwitch<size_t>(switches, "cubic-fit", &inner_knots);
    bool retain_boundaries = parseSwitch<bool>(switches, "retain-boundaries");
    bool retain_blends = parseSwitch<bool>(switches, "retain-direction-blends");
    bool retain_ribbons = parseSwitch<bool>(switches, "retain-ribbons");

    if (cubic && (retain_blends || retain_boundaries)) {
      std::cerr << "ERROR: Direction blends / boundaries cannot be retained for cubic fit"
                << std::endl;
      retain_blends = false;
    }
    if (cubic && retain_ribbons) {
      std::cerr << "ERROR: Ribbons cannot be retained for cubic fit (yet)" << std::endl;
      retain_ribbons = false;
    }

    for (size_t i = 0; i < quads.size(); ++i) {
      const auto &quad = quads[i];
      if (cubic) {
        result[i] = BSSurface(3, 3, PointVector(16));
        for (size_t j = 1; j <= inner_knots; ++j) {
          result[i] = result[i].insertKnotU((double)j / (inner_knots + 1), 1);
          result[i] = result[i].insertKnotV((double)j / (inner_knots + 1), 1);
        }
      }
      auto ncp = result[i].numControlPoints();
      auto constraint = [&](size_t i, size_t j) -> MoveConstraint {
        if (retain_blends && (i < 2 || j < 2 || i >= ncp[0] - 2 || j >= ncp[1] - 2))
          return MoveType::Fixed();
        if (retain_ribbons &&
            ((quad.boundaries[0].on_ribbon && i < 2) ||
             (quad.boundaries[1].on_ribbon && j < 2) ||
             (quad.boundaries[2].on_ribbon && i >= ncp[0] - 2) ||
             (quad.boundaries[3].on_ribbon && j >= ncp[1] - 2)))
          return MoveType::Fixed();
        if (retain_boundaries && (i < 1 || j < 1 || i >= ncp[0] - 1 || j >= ncp[1] - 1))
          return MoveType::Fixed();
        return MoveType::Free();
      };
      if (retain_blends)
        for (size_t j = 0; j < 5; ++j)
          bsplineFit(result[i], quad.resolution, quad.samples, constraint, 0, true);
      else
        bsplineFit(result[i], quad.resolution, quad.samples, constraint, 0, false);
    }
  }

  // 10. Fix C0 - TODO: knot vectors may need to be unified
  if (parseSwitch<bool>(switches, "fix-c0-inside")) {
    for (size_t si = 0; si < adjacency.size(); ++si) {
      const auto &adj = adjacency[si];
      if (adj.size() == 2) {
        const auto &q1 = quads[adj[0].first];
        const auto &q2 = quads[adj[1].first];
        auto &s1 = result[adj[0].first];
        auto &s2 = result[adj[1].first];
        size_t b1 = adj[0].second;
        size_t b2 = adj[1].second;
        bool r1 = q1.boundaries[b1].reversed;
        bool r2 = q2.boundaries[b2].reversed;
        auto setIndex = [](const BSSurface &s, size_t b, size_t &j, size_t &n) {
          size_t m;
          if (b == 0 || b == 2) {
            n = s.numControlPoints()[1];
            m = s.numControlPoints()[0];
          } else {
            n = s.numControlPoints()[0];
            m = s.numControlPoints()[1];
          }
          if (b < 2)
            j = 0;
          else
            j = m - 1;
        };
        size_t j1, j2, n1, n2;
        setIndex(s1, b1, j1, n1);
        setIndex(s2, b2, j2, n2);
        assert(n1 == n2);
        for (size_t i1 = 0; i1 < n1; ++i1) {
          size_t i2 = r1 == r2 ? i1 : n1 - i1 - 1;
          auto &p1 = b1 % 2 == 0 ? s1.controlPoint(j1, i1) : s1.controlPoint(i1, j1);
          auto &p2 = b2 % 2 == 0 ? s2.controlPoint(j2, i2) : s2.controlPoint(i2, j2);
          if (i1 == 0 || i1 == n1 - 1) {
            const auto &v1 = vertices[endpoints[si].first];
            const auto &v2 = vertices[endpoints[si].second];
            p1 = (i1 == 0 && !r1) || (i1 != 0 && r1) ? v1 : v2;
            p2 = (i2 == 0 && !r2) || (i2 != 0 && r2) ? v1 : v2;
          } else {
            auto m = (p1 + p2) / 2;
            p1 = m;
            p2 = m;
          }
        }
      }
    }
  }

  // TODO: only works for polynomial ribbons
  if (parseSwitch<bool>(switches, "fix-c0-outside")) {
    for (size_t i = 0; i < quads.size(); ++i) {
      const auto &quad = quads[i];
      for (size_t side = 0; side < 4; ++side) {
        const auto &b = quad.boundaries[side];
        if (!b.on_ribbon)
          continue;
        auto r = ribbons[b.ribbon][0];
        r = r.insertKnot(b.s0, 3);
        r = r.insertKnot(b.s1, 3);
        PointVector pv;
        if (b.s0 != 0 && b.s1 != 0)
          pv.insert(pv.end(), r.controlPoints().begin() + 3, r.controlPoints().begin() + 7);
        else
          pv.insert(pv.end(), r.controlPoints().begin(), r.controlPoints().begin() + 4);
        if (b.s0 > b.s1)
          std::reverse(pv.begin(), pv.end());
        r = BSCurve(pv);
        size_t k = 0;
        if (side == 2)
          k = result[i].numControlPoints()[0] - 1;
        if (side == 3)
          k = result[i].numControlPoints()[1] - 1;
        if (side == 0 || side == 2) {
          const auto &basis = result[i].basisV();
          for (auto k : basis.knots())
            r = r.insertKnot(k, 1);
          for (size_t j = 0; j < r.controlPoints().size(); ++j)
            result[i].controlPoint(k, j) = r.controlPoints()[j];
        } else {
          const auto &basis = result[i].basisU();
          for (auto k : basis.knots())
            r = r.insertKnot(k, 1);
          for (size_t j = 0; j < r.controlPoints().size(); ++j)
            result[i].controlPoint(j, k) = r.controlPoints()[j];
        }
      }
    }
  }

  if (parseSwitch<bool>(switches, "print-continuity-errors"))
    printContinuityErrors(result);

  if (parseSwitch<bool>(switches, "print-approximation-errors"))
    printApproximationErrors(result);

  return result;
}
