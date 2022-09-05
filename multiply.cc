#include <algorithm>
#include <cassert>
#include <limits>

#include <Eigen/Dense>

#include <bezier-extractions.hh>

#include "multiply.hh"

using namespace Geometry;

// Generated by newton-cotes.jl
// n : 1 - 12; odd n -> exact for degree n; even n -> exact for degree n + 1
static const DoubleMatrix newton_cotes =
  {{},
   {5.0e-1,5.0e-1},
   {3.333333333333333e-1,1.333333333333333e0,3.333333333333333e-1},
   {3.75e-1,1.125e0,1.125e0,3.75e-1},
   {3.111111111111111e-1,1.422222222222222e0,5.333333333333333e-1,
    1.422222222222222e0,3.111111111111111e-1},
   {3.298611111111111e-1,1.302083333333333e0,8.680555555555556e-1,
    8.680555555555556e-1,1.302083333333333e0,3.298611111111111e-1},
   {2.928571428571429e-1,1.542857142857143e0,1.928571428571429e-1,
    1.942857142857143e0,1.928571428571429e-1,1.542857142857143e0,
    2.928571428571429e-1},
   {3.04224537037037e-1,1.449016203703704e0,5.359375e-1,
    1.210821759259259e0,1.210821759259259e0,5.359375e-1,
    1.449016203703704e0,3.04224537037037e-1},
   {2.790828924162258e-1,1.661516754850088e0,-2.618694885361552e-1,
    2.961834215167549e0,-1.281128747795415e0,2.961834215167549e0,
    -2.618694885361552e-1,1.661516754850088e0,2.790828924162258e-1},
   {2.869754464285714e-1,1.581127232142857e0,1.084821428571429e-1,
    1.943035714285714e0,5.803794642857143e-1,5.803794642857143e-1,
    1.943035714285714e0,1.084821428571429e-1,1.581127232142857e0,
    2.869754464285714e-1},
   {2.683414836192614e-1,1.775359414248303e0,-8.104357062690396e-1,
    4.549462882796216e0,-4.351551226551227e0,7.137646304312971e0,
    -4.351551226551227e0,4.549462882796216e0,-8.104357062690396e-1,
    1.775359414248303e0,2.683414836192614e-1},
   {2.742655400315991e-1,1.703408943727954e0,-4.088615497317754e-1,
    3.186240802744709e0,-1.211958980930335e0,1.956905244157848e0,
    1.956905244157848e0,-1.211958980930335e0,3.186240802744709e0,
    -4.088615497317754e-1,1.703408943727954e0,2.742655400315991e-1},
   {2.596738499595643e-1,1.88443328100471e0,-1.443863565006422e0,
    6.797986775129632e0,-9.798067646639075e0,1.665311602683031e1,
    -1.670655744255744e1,1.665311602683031e1,-9.798067646639075e0,
    6.797986775129632e0,-1.443863565006422e0,1.88443328100471e0,
    2.596738499595643e-1}};

[[maybe_unused]]
static double newtonCotes(const std::function<double(double)> &f, double a, double b,
                          size_t n) {
  assert(n <= 12);
  double sum = 0;
  double h = (b - a) / n;
  for (size_t j = 0; j <= n; ++j)
    sum += f(a + h * j) * newton_cotes[n][j];
  return sum * h;
}

static double trapezoid(const std::function<double(double)> &f, double a, double b,
                        size_t n, double s = 0) {
  if (n == 1)
    return (f(a) + f(b)) / 2 * (b - a);
  size_t k = 1 << (n - 2);
  double h = (b - a) / k;
  double x = a + h / 2;
  double sum = 0;
  for (size_t i = 0; i < k; ++i)
    sum += f(x + h * i);
  return s + h * sum / 2;
}

[[maybe_unused]]
static double simpson(const std::function<double(double)> &f, double a, double b,
                      size_t max_iterations, double tolerance) {
  double s;
  double s_prev = -std::numeric_limits<double>::max(), st_prev = s_prev;
  for (size_t i = 1; i <= max_iterations; ++i) {
    double st = trapezoid(f, a, b, i, st_prev);
    s = (st * 4 - st_prev) / 3;
    if (i > 5 && // avoid exiting too early
        (std::abs(s - s_prev) < tolerance * std::abs(s_prev) ||
         (s == 0 && s_prev == 0)))
      break;
    s_prev = s;
    st_prev = st;
  }
  return s;
}

[[maybe_unused]]
static double integrate(const std::function<double(double)> &f,
                        const DoubleVector &intervals, size_t degree) {
  double result = 0;
  for (size_t i = 1; i < intervals.size(); ++i)
    result += newtonCotes(f, intervals[i-1], intervals[i], degree);
    // result += simpson(f, intervals[i-1], intervals[i], 5, 1e-3);
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

[[maybe_unused]]
static PointVector inBasis(const BSBasis &basis, const std::function<Point3D(double)> &f) {
  const auto &knots = basis.knots();
  size_t d = basis.degree();
  size_t n = knots.size() - d - 1;
  DoubleVector params;

  // TODO: this works only in this special case
  assert(d == 6);
  for (size_t i = d; i < n; ++i) {
    double k1 = knots[i];
    double k2 = knots[i+1];
    if (k1 == k2)
      continue;
    size_t c = (i == d || i == n - 1) ? 5 : 4;
    for (size_t j = 0; j < c; ++j) {
      double alpha = (double)j / c;
      params.push_back(k1 * (1 - alpha) + k2 * alpha);
    }
  }
  params.push_back(knots.back());

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n), b(n, 3);
  DoubleVector coeff;
  for (size_t i = 0; i < n; ++i) {
    double u = params[i];
    size_t span = basis.findSpan(u);
    basis.basisFunctions(span, u, coeff);
    for (size_t j = 0; j <= d; ++j)
      A(i, span + j - d) = coeff[j];
    for (size_t j = 0; j < 3; ++j)
      b(i, j) = f(u)[j];
  }
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  PointVector cpts(n);
  for (size_t i = 0; i < n; ++i)
    cpts[i] = Vector3D(x(i, 0), x(i, 1), x(i, 2));
  return cpts;
}

// The resulting surface is linear in v, with its point- and derivative computed as given.
BSSurface multiplyBSplines(const BSBasis &basis1, const BSBasis &basis2,
                           const std::function<Point3D(double)> &point,
                           const std::function<Vector3D(double)> &derivative) {
  auto basis = combineBases(basis1, basis2);
  size_t d = basis.degree();
  const auto &knots = basis.knots();

#ifdef INTERPOLATION_METHOD
  auto pv1 = inBasis(basis, point);
  auto pv2 = inBasis(basis, derivative);
  size_t n_cp = pv1.size();
  PointVector cp(2 * n_cp);
  for (size_t i = 0; i < n_cp; ++i) {
    cp[2*i] = pv1[i];
    cp[2*i+1] = pv1[i] + pv2[i];
  }
  return { d, 1, knots, { 0, 0, 1, 1 }, cp };
#endif

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

static std::vector<size_t> computeSpans(const BSBasis &basis) {
  const auto &knots = basis.knots();
  std::vector<size_t> spans;
  for (size_t i = 0; i < knots.size() - 1; ++i)
    if (knots[i] != knots[i+1])
      spans.push_back(i);
  return spans;
}

static std::vector<PointVector> segments(const BSCurve &curve) {
  const auto &basis = curve.basis();
  size_t d = basis.degree();
  const auto &knots = basis.knots();
  const auto &cpts = curve.controlPoints();
  auto spans = computeSpans(basis);
  size_t L = spans.size();

  auto C = bezierExtractionMatrices(d, knots);

  std::vector<PointVector> result;
  for (size_t i = 0; i < L; ++i) {
    PointVector pv;
    for (size_t j = 0; j <= d; ++j) {
      Point3D p(0, 0, 0);
      for (size_t k = 0; k <= d; ++k)
        p += cpts[spans[i]+k-d] * C[i](k, j);
      pv.push_back(p);
    }
    result.push_back(pv);
  }

  return result;
}

static PointVector elevateBezier(const PointVector &cpts) {
  size_t n = cpts.size();
  PointVector elevated;
  elevated.push_back(cpts.front());
  for (size_t i = 1; i < n; ++i) {
    double alpha = (double)i / n;
    elevated.push_back(cpts[i-1] * alpha + cpts[i] * (1 - alpha));
  }
  elevated.push_back(cpts.back());
  return elevated;
}

static size_t binomial(size_t n, size_t k) {
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

// Computes c(u) * s(u)[i] in Bezier form
static PointVector bezierProduct(const PointVector &c, const PointVector &s, size_t i) {
  // As in L. Piegl, W. Tiller: Symbolic operators for NURBS. CAD 29(5), pp. 361-368, 1997.
  size_t p = c.size() - 1, q = s.size() - 1;
  PointVector result;
  for (size_t k = 0; k <= p + q; ++k) {
    Point3D sum(0, 0, 0);
    size_t lstart = k > q ? k - q : 0;
    for (size_t l = lstart; l <= std::min(p, k); ++l)
      sum += c[l] * s[k-l][i] * binomial(p, l) * binomial(q, k - l) / binomial(p + q, k);
    result.push_back(sum);
  }
  return result;
}

static BSCurve reconstructSpline(const std::vector<PointVector> &segments, const BSBasis &basis) {
  size_t d = basis.degree();
  const auto &knots = basis.knots();
  auto spans = computeSpans(basis);

  auto C = bezierExtractionMatrices(d, knots);

  PointVector cpts;
  for (size_t k = 0; k < knots.size() - d - 1; ++k) {
    Point3D p(0, 0, 0);
    size_t i = 0;
    while (spans[i] < k)
      ++i;
    auto Cinv = C[i].inverse();
    for (size_t j = 0; j <= d; ++j)
      p += segments[i][j] * Cinv(j, k + d - spans[i]);
    cpts.push_back(p);
  }

  return { basis.degree(), basis.knots(), cpts };
}

BSSurface multiplyBSplinesByBezier(BSCurve point, BSCurve crossder, BSCurve scaling) {
  auto basis = combineBases(point.basis(), scaling.basis());
  size_t d = basis.degree();

  // Insert missing knots for correct segmentation
  for (double k : basis.knots()) {
    if (std::find(point.basis().knots().begin(),
                  point.basis().knots().end(),
                  k) ==
        point.basis().knots().end()) {
      point = point.insertKnot(k, 1);
      crossder = crossder.insertKnot(k, 1);
    }
    if (std::find(scaling.basis().knots().begin(),
                  scaling.basis().knots().end(),
                  k) ==
        scaling.basis().knots().end())
      scaling = scaling.insertKnot(k, 1);
  }

  auto point_segs = segments(point);
  auto crossder_segs = segments(crossder);
  auto scaling_segs = segments(scaling);

  // Compute the elevated derivative of point
  std::vector<PointVector> elevated_der_segs;
  auto point_spans = computeSpans(point.basis());
  size_t point_deg = point_segs[0].size() - 1;
  for (size_t j = 0; j < point_segs.size(); ++j) {
    const auto &ps = point_segs[j];
    PointVector cpts;
    double interval =
      point.basis().knots()[point_spans[j]+1] -
      point.basis().knots()[point_spans[j]];
    for (size_t i = 1; i < ps.size(); ++i)
      cpts.push_back((ps[i] - ps[i-1]) * point_deg / interval);
    elevated_der_segs.push_back(elevateBezier(cpts));
  }

  // Compute the product for each Bezier segment
  size_t L = point_segs.size();
  for (size_t i = 0; i < L; ++i) {
    for (size_t j = 0; j < d - point_deg; ++j)
      point_segs[i] = elevateBezier(point_segs[i]);
    crossder_segs[i] = bezierProduct(crossder_segs[i], scaling_segs[i], 0);
    elevated_der_segs[i] = bezierProduct(elevated_der_segs[i], scaling_segs[i], 1);
    for (size_t j = 0; j <= d; ++j)
      crossder_segs[i][j] += point_segs[i][j] + elevated_der_segs[i][j];
  }

  // Reconstruct the B-Spline curves
  auto c0 = reconstructSpline(point_segs, basis);
  auto c1 = reconstructSpline(crossder_segs, basis);

  // Build the final surface
  size_t n = c0.controlPoints().size();
  PointVector cpts(2 * n);
  BSSurface result(d, 1, basis.knots(), { 0, 0, 1, 1 }, cpts);
  for (size_t i = 0; i < n; ++i) {
    result.controlPoint(i, 0) = c0.controlPoints()[i];
    result.controlPoint(i, 1) = c1.controlPoints()[i];
  }
  return result;
}
