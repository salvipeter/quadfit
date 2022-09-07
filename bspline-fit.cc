#include <map>

#include <Eigen/Dense>

#include "bspline-fit.hh"

using namespace Geometry;

using IndexPair = std::pair<size_t, size_t>;
using VecMap = Eigen::Map<const Eigen::Vector3d>;

void bsplineFit(BSCurve &curve, size_t resolution, const PointVector &samples,
                const std::function<bool(size_t)> &isFixed, double smoothness) {
  const auto &basis = curve.basis();
  size_t p = basis.degree();
  auto &cpts = curve.controlPoints();
  auto n = cpts.size();

  std::map<size_t, size_t> index_map;
  size_t index = 0;
  for (size_t i = 0; i < n; ++i)
    if (!isFixed(i)) {
      index_map[i] = index;
      index++;
    }
  size_t nvars = index;

  size_t n_smoothing = n - 1;
  Eigen::MatrixXd A(samples.size() + n_smoothing, nvars);
  Eigen::MatrixXd b(samples.size() + n_smoothing, 3);

  A.setZero(); b.setZero();

  auto addValue = [&](size_t row, size_t i, double x) {
    if (isFixed(i))
      b.row(row) -= VecMap(cpts[i].data()) * x;
    else
      A(row, index_map.at(i)) = x;
  };

  index = 0;
  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    size_t span = basis.findSpan(u);
    DoubleVector coeff;
    basis.basisFunctions(span, u, coeff);
    b.row(index) = VecMap(samples[index].data());
    for (size_t k = 0; k <= p; ++k) {
      size_t i1 = span - p + k;
      addValue(index, i1, coeff[k]);
    }
    index++;
  }

  // Smoothness terms
  for (size_t i = 1; i < n - 1; ++i) {
    addValue(index, i, smoothness);
    addValue(index, i - 1, -0.5 * smoothness);
    addValue(index, i + 1, -0.5 * smoothness);
    index++;
  }

  Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);

  for (size_t i = 0; i < n; ++i)
    if (!isFixed(i)) {
      size_t k = index_map.at(i);
      cpts[i] = Point3D(x(k, 0), x(k, 1), x(k, 2));
    }
}

void bsplineFit(BSSurface &surface, size_t resolution, const PointVector &samples,
                const std::function<bool(size_t,size_t)> &isFixed, double smoothness) {
  size_t pu = surface.basisU().degree(), pv = surface.basisV().degree();
  auto [nu, nv] = surface.numControlPoints();
  auto &cpts = surface.controlPoints();

  std::map<IndexPair, size_t> index_map;
  size_t index = 0;
  for (size_t i = 0; i < nu; ++i)
    for (size_t j = 0; j < nv; ++j)
      if (!isFixed(i, j)) {
        index_map[{i,j}] = index;
        index++;
      }
  size_t nvars = index;

  size_t n_smoothing = (nu - 1) * (nv - 1);
  Eigen::MatrixXd A(samples.size() + n_smoothing, nvars);
  Eigen::MatrixXd b(samples.size() + n_smoothing, 3);

  A.setZero(); b.setZero();

  auto addValue = [&](size_t row, size_t i, size_t j, double x) {
    if (isFixed(i, j))
      b.row(row) -= VecMap(cpts[i*nv+j].data()) * x;
    else
      A(row, index_map.at({i,j})) = x;
  };

  index = 0;
  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    size_t span_u = surface.basisU().findSpan(u);
    DoubleVector coeff_u;
    surface.basisU().basisFunctions(span_u, u, coeff_u);
    for (size_t j = 0; j <= resolution; ++j) {
      double v = (double)j / resolution;
      size_t span_v = surface.basisV().findSpan(v);
      DoubleVector coeff_v;
      surface.basisV().basisFunctions(span_v, v, coeff_v);

      b.row(index) = VecMap(samples[index].data());
      for (size_t k = 0; k <= pu; ++k)
        for (size_t l = 0; l <= pv; ++l) {
          size_t i1 = span_u - pu + k, j1 = span_v - pv + l;
          addValue(index, i1, j1, coeff_u[k] * coeff_v[l]);
        }

      index++;
    }
  }

  // Smoothness terms
  for (size_t i = 1; i < nu - 1; ++i)
    for (size_t j = 1; j < nv - 1; ++j) {
      addValue(index, i, j, smoothness);
      addValue(index, i - 1, j, -0.25 * smoothness);
      addValue(index, i + 1, j, -0.25 * smoothness);
      addValue(index, i, j - 1, -0.25 * smoothness);
      addValue(index, i, j + 1, -0.25 * smoothness);
      index++;
    }

  Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);

  index = 0;
  for (size_t i = 0; i < nu; ++i)
    for (size_t j = 0; j < nv; ++j) {
      if (!isFixed(i, j)) {
        size_t k = index_map.at({i,j});
        cpts[index] = Point3D(x(k, 0), x(k, 1), x(k, 2));
      }
      index++;
    }
}
