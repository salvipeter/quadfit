#include <Eigen/Dense>

#include "multiply.hh"

using namespace Geometry;

static const DoubleMatrix newton_cotes =
  { // n : 1 - 12; odd n -> exact for degree n; even n -> exact for degree n + 1
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
static double integrate(const std::function<double(double)> &f,
                        const DoubleVector &intervals,
                        size_t degree) {
  double result = 0;
  size_t n = std::min(12, std::max(1, (int)degree));
  for (size_t i = 1; i < intervals.size(); ++i) {
    double sum = 0;
    double h = (intervals[i] - intervals[i-1]) / n;
    for (size_t j = 0; j <= n; ++j)
      sum += f(intervals[i-1] + h * j) * newton_cotes[n][j];
    result += sum * h;
  }
  return result;
}

static BSBasis combineBases(const BSBasis &basis1, const BSBasis &basis2) {
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

// The resulting surface is linear in v, with its point- and derivative computed as given.
BSSurface multiplyBSplines(const BSBasis &basis1, const BSBasis &basis2,
                           const std::function<Point3D(double)> &point,
                           const std::function<Vector3D(double)> &derivative) {
  auto basis = combineBases(basis1, basis2);
  size_t d = basis.degree();
  auto N = [&](size_t i, double u) {
    size_t span = basis.findSpan(u);
    if (i + d < span || i > span)
      return 0.0;
    DoubleVector coeff;
    basis.basisFunctions(span, u, coeff);
    return coeff[i+d-span];
  };
  const auto &h1 = point;
  const auto &h2 = derivative;
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
  PointVector cpts(2 * (n + 1));
  BSSurface result(d, 1, knots, { 0, 0, 1, 1 }, cpts);
  for (size_t i = 0; i <= n; ++i) {
    result.controlPoint(i, 0) = Vector3D(x(i, 0), x(i, 1), x(i, 2));
    result.controlPoint(i, 1) = Vector3D(x(i, 0) + x(i, 3), x(i, 1) + x(i, 4), x(i, 2) + x(i, 5));
  }
  return result;
}
