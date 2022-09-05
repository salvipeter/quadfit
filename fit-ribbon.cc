#include <algorithm>
#include <iterator>

#include <Eigen/Dense>

#include "fit-ribbon.hh"
#include "multiply.hh"

using namespace Geometry;

static BSCurve interpolateWithZeroTangents(const Point2DVector &values) {
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

static BSSurface splitU(BSSurface surface, double from, double to) {
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

std::vector<BSSurface> fitSlices(const BSSurface &ribbon, const Point2DVector &sh) {
  auto alpha = interpolateWithZeroTangents(sh);

  // Indirect solution using integration or interpolation
  // ---
  // auto point = [&](double u) {
  //   return ribbon.eval(u, 0.0);
  // };
  // auto derivative = [&](double u) {
  //   VectorMatrix der;
  //   ribbon.eval(u, 0.0, 1, der);
  //   return der[0][1] * alpha.eval(u)[0];
  // };
  // auto sextic = multiplyBSplines(ribbon.basisU(), alpha.basis(), point, derivative);

  // Direct solution on the Bezier segments
  // ---
  PointVector cpts1, cpts2;
  size_t n = ribbon.numControlPoints()[0];
  size_t cross_deg = ribbon.basisV().degree();
  for (size_t i = 0; i < n; ++i) {
    cpts1.push_back(ribbon.controlPoint(i, 0));
    cpts2.push_back((ribbon.controlPoint(i, 1) - ribbon.controlPoint(i, 0)) * cross_deg);
  }
  BSCurve c(ribbon.basisU().degree(), ribbon.basisU().knots(), cpts1);
  BSCurve d(ribbon.basisU().degree(), ribbon.basisU().knots(), cpts2);
  auto sextic = multiplyBSplinesByBezier(c, d, alpha);

  std::vector<BSSurface> result;
  for (size_t i = 1; i < sh.size(); ++i)
    result.push_back(splitU(sextic, sh[i-1][0], sh[i][0]));
  
  return result;
}
