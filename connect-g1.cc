#include <algorithm>
#include <fstream>
#include <functional>
#include <iterator>

#include <geometry.hh>

#include <Eigen/Dense>

#include "io.hh"

using namespace Geometry;

static const DoubleMatrix newton_cotes =
  { // n : 1 - 12
    { },
    { 0.5, 0.5 },
    { 0.33333333, 1.33333333, 0.33333333 },
    { 0.375, 1.125, 1.125, 0.375 },
    { 0.31111111, 1.42222222, 0.53333333, 1.42222222, 0.31111111 },
    { 0.32986111, 1.30208333, 0.86805556, 0.86805556, 1.30208333,
      0.32986111 },
    { 0.29285714, 1.54285714, 0.19285714, 1.94285714, 0.19285714,
      1.54285714, 0.29285714 },
    { 0.30422454, 1.4490162 , 0.5359375 , 1.21082176, 1.21082176,
      0.5359375 , 1.4490162 , 0.30422454 },
    { 0.27908289, 1.66151675, -0.26186949, 2.96183422, -1.28112875,
      2.96183422, -0.26186949, 1.66151675, 0.27908289 },
    { 0.28697545, 1.58112723, 0.10848214, 1.94303571, 0.58037946,
      0.58037946, 1.94303571, 0.10848214, 1.58112723, 0.28697545 },
    { 0.26834148, 1.77535941, -0.81043571, 4.54946288, -4.35155123,
      7.1376463, -4.35155123, 4.54946288, -0.81043571, 1.77535941,
      0.26834148 },
    { 0.27426554, 1.70340894, -0.40886155, 3.1862408, -1.21195898,
      1.95690524, 1.95690524, -1.21195898, 3.1862408, -0.40886155,
      1.70340894, 0.27426554 },
    { 0.25967385, 1.88443328, -1.44386357, 6.79798678,
      -9.79806765, 16.65311603, -16.70655744, 16.65311603,
      -9.79806765, 6.79798678, -1.44386357, 1.88443328,
      0.25967385 },
  };
double integrate(const std::function<double(double)> &f, const DoubleVector &intervals,
                 size_t degree) {
  double result = 0;
  size_t n = std::min(12, std::max(1, (int)degree - 1));
  for (size_t i = 1; i < intervals.size(); ++i) {
    double sum = 0;
    double h = (intervals[i] - intervals[i-1]) / n;
    for (size_t j = 0; j <= n; ++j)
      sum += f(intervals[i-1] + h * j) * newton_cotes[n][j];
    result += sum * h;
  }
  return result;
}

BSBasis combineBases(const BSBasis &basis1, const BSBasis &basis2) {
  size_t d1 = basis1.degree();
  size_t d2 = basis2.degree();
  size_t d = d1 + d2;
  const auto &knots1 = basis1.knots();
  const auto &knots2 = basis2.knots();
  size_t n1 = knots1.size() - d1 - 1, n2 = knots2.size() - d2 - 1;
  DoubleVector knots;
  std::fill_n(std::back_inserter(knots), d + 1, knots1.front());
  size_t i1 = d1 + 1, i2 = d2 + 1;
  while (i1 < n1 || i2 < n2) {
    if (knots1[i1] < knots2[i2]) {
      double k = knots1[i1];
      size_t m = 0;
      while (knots1[i1] == k) { m++; i1++; }
      std::fill_n(std::back_inserter(knots), m + d2, k);
    } else if (knots2[i2] < knots1[i1]) {
      double k = knots2[i2];
      size_t m = 0;
      while (knots2[i2] == k) { m++; i2++; }
      std::fill_n(std::back_inserter(knots), m + d1, k);
    } else {
      double k = knots1[i1];
      size_t m1 = 0, m2 = 0;
      while (knots1[i1] == k) { m1++; i1++; }
      while (knots2[i2] == k) { m2++; i2++; }
      std::fill_n(std::back_inserter(knots), d1 + d2 - std::min(d1 - m1, d2 - m2), k);
    }
  }
  std::fill_n(std::back_inserter(knots), d + 1, knots1.back());
  return { d, knots };
}

// coordinates of u when decomposed in the (v,w) system
// (Vector3D instead of Vector2D for convenience)
Vector3D inSystem(const Vector3D &u, const Vector3D &v, const Vector3D &w) {
  auto uv = u ^ v, uw = u ^ w, vw = v ^ w;
  return {
    uw.norm() / vw.norm() * (uw * vw > 0 ? 1 : -1),
    uv.norm() / vw.norm() * (uv * vw < 0 ? 1 : -1),
    0
  };
}

BSSurface connectG1(const BSCurve &c, const BSCurve &c1, const BSCurve &c2) {
  Point3D S_0, S_1, S1_0, S1_1, S2_0, S2_1;
  Vector3D Su_0, Su_1, Suu_0, Suu_1;
  Vector3D S1v_0, S2v_0, S1v_1, S2v_1;
  Vector3D S1uv_0, S2uv_0, S1uv_1, S2uv_1;

  double lo = c.basis().low(), hi = c.basis().high();

  VectorVector der;
  S_0 = c.eval(lo, 2, der);
  Su_0 = der[1];
  Suu_0 = der[2];
  S_1 = c.eval(hi, 2, der);
  Su_1 = der[1];
  Suu_1 = der[2];
  S1_0 = c1.eval(lo, 1, der);
  S1v_0 = S1_0 - S_0;
  S1uv_0 = der[1] - Su_0;
  S2_0 = c2.eval(lo, 1, der);
  S2v_0 = S2_0 - S_0;
  S2uv_0 = der[1] - Su_0;
  S1_1 = c1.eval(hi, 1, der);
  S1v_1 = S1_1 - S_1;
  S1uv_1 = der[1] - Su_1;
  S2_1 = c2.eval(hi, 1, der);
  S2v_1 = S2_1 - S_1;
  S2uv_1 = der[1] - Su_1;

  // Sv = D * alpha + Su * beta
  auto D0 = (S1v_0 - S2v_0) / 2;
  auto D1 = (S1v_1 - S2v_1) / 2;
  auto ab_0 = inSystem(S1v_0, D0, Su_0);
  auto ab_1 = inSystem(S1v_1, D1, Su_1);
  auto ab_0d = inSystem(S1uv_0 - (S1uv_0 - S2uv_0) / 2 * ab_0[0] - Suu_0 * ab_0[1], D0, Su_0);
  auto ab_1d = inSystem(S1uv_1 - (S1uv_1 - S2uv_1) / 2 * ab_1[0] - Suu_1 * ab_1[1], D1, Su_1);
  BSCurve scaling(3, { lo, lo, lo, lo, hi, hi, hi, hi },
                  { ab_0, ab_0 + ab_0d / 3, ab_1 - ab_1d / 3, ab_1 });

  // Same code as in fit-ribbon.cc, just using scaling, and a different h1 & h2
  auto basis = combineBases(c.basis(), scaling.basis());
  size_t d = basis.degree();
  auto N = [&](size_t i, double u) {
    size_t span = basis.findSpan(u);
    if (i + d < span || i > span)
      return 0.0;
    DoubleVector coeff;
    basis.basisFunctions(span, u, coeff);
    return coeff[i+d-span];
  };
  auto h1 = [&](double u) {
    return c.eval(u);
  };
  auto h2 = [&](double u) {
    VectorVector der;
    c.eval(u, 1, der);
    auto p1 = c1.eval(u);
    auto p2 = c2.eval(u);
    auto s = scaling.eval(u);
    return (p1 - p2) / 2 * s[0] + der[1] * s[1];
  };
  const auto &knots = basis.knots();
  DoubleVector intervals = { knots.front() };
  for (double k : knots)
    if (k > intervals.back())
      intervals.push_back(k);
  size_t n = knots.size() - d - 2;
  Eigen::MatrixXd A(n + 1, n + 1), b(n + 1, 6);
  for (size_t i = 0; i <= n; ++i) {
    for (size_t j = 0; j <= n; ++j)
      A(i, j) = integrate([&](double u) { return N(i, u) * N(j, u); }, intervals, 2 * d);
    for (size_t j = 0; j < 3; ++j) {
      b(i, j) = integrate([&](double u) { return h1(u)[j] * N(i, u); }, intervals, 2 * d);
      b(i, 3 + j) = integrate([&](double u) { return h2(u)[j] * N(i, u); }, intervals, 2 * d);
    }
  }
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  PointVector cpts((n + 1) * 2);
  BSSurface result(d, 1, knots, { 0, 0, 1, 1 }, cpts);
  for (size_t i = 0; i <= n; ++i) {
    result.controlPoint(i, 0) = Vector3D(x(i, 0), x(i, 1), x(i, 2));
    result.controlPoint(i, 1) = Vector3D(x(i, 0) + x(i, 3), x(i, 1) + x(i, 4), x(i, 2) + x(i, 5));
  }
  return result;
}


// Testing

int main(int argc, char **argv) {
  DoubleVector knots = { 0, 0, 0, 0, 1, 3, 4, 4, 4, 4 };
  BSCurve c(3, knots, {{0,0,0},{4,1,0},{8,1,0},{12,0,1},{16,0,1},{19,1,0}});
  BSCurve c1(3, knots, {{-2,3,0},{3,5,1},{8,5,2},{12,4,1},{16,4,0},{20,4,0}});
  BSCurve c2(3, knots, {{-1,-4,0},{4,-3,1},{8,-3,2},{13,-4,1},{17,-3,1},{21,-2,0}});
  size_t resolution = 50;
  std::vector<BSSurface> surfaces = { connectG1(c, c1, c2), connectG1(c, c2, c1) };
  surfaces[1].reverseV();
  writeSTL(surfaces, "/tmp/surfaces.stl", resolution);
  writeControlNet(surfaces, "/tmp/controls.obj");
}
