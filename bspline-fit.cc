#include <map>

#include <Eigen/Dense>

#include "bspline-fit.hh"

using namespace Geometry;

using IndexPair = std::pair<size_t, size_t>;
using VecMap = Eigen::Map<const Eigen::Vector3d>;

// Least-squares fit with Lagrange multipliers:
//
// ||Ax - b||^2 -> min   satisfying   Cx = d
//
// [ A^t A   C^t ]  [ x ]     [ A^t b ]
// [             ]  [   ]  =  [       ]
// [   C      0  ]  [ L ]     [   d   ]

void bsplineFit(BSCurve &curve, const PointVector &samples,
                const std::function<MoveConstraint(size_t)> &constraint,
                double smoothness) {
  const auto &basis = curve.basis();
  size_t p = basis.degree();
  auto &cpts = curve.controlPoints();
  size_t n = cpts.size();
  size_t resolution = samples.size() - 1;

  std::map<size_t, size_t> index_map, tangent_map;
  size_t index = 0, n_tconstr = 0;
  for (size_t i = 0; i < n; ++i)
    if (!std::holds_alternative<MoveType::Fixed>(constraint(i))) {
      index_map[i] = index;
      index++;
      if (std::holds_alternative<MoveType::Tangent>(constraint(i))) {
        tangent_map[i] = n_tconstr;
        n_tconstr++;
      }
    }
  size_t nvars = index;

  size_t n_smoothing = smoothness == 0 ? 0 : n - 1;
  size_t n_rows = samples.size() + n_smoothing;
  Eigen::MatrixXd A(n_rows * 3, nvars * 3), C(n_tconstr, nvars * 3);
  Eigen::VectorXd b(n_rows * 3), d(n_tconstr);

  A.setZero(); C.setZero();

  auto addValue = [&](size_t row, size_t i, double x) {
    if (std::holds_alternative<MoveType::Fixed>(constraint(i)))
      for (size_t j = 0; j < 3; ++j)
        b(row + j) -= cpts[i][j] * x;
    else
      for (size_t j = 0; j < 3; ++j)
        A(row + j, 3 * index_map.at(i) + j) = x;
  };

  for (size_t i = 0; i <= resolution; ++i, ++index) {
    double u = (double)i / resolution;
    size_t span = basis.findSpan(u);
    DoubleVector coeff;
    basis.basisFunctions(span, u, coeff);
    for (size_t j = 0; j < 3; ++j)
      b(3 * i + j) = samples[i][j];
    for (size_t k = 0; k <= p; ++k) {
      size_t i1 = span - p + k;
      addValue(3 * i, i1, coeff[k]);
    }
  }

  if (smoothness != 0)
    for (size_t i = 1; i < n - 1; ++i) {
      addValue(index, i, smoothness);
      addValue(index, i - 1, -0.5 * smoothness);
      addValue(index, i + 1, -0.5 * smoothness);
      index += 3;
    }

  for (auto [i, row] : tangent_map) {
    auto normal = std::get<MoveType::Tangent>(constraint(i)).normal;
    for (size_t j = 0; j < 3; ++j)
      C(row, 3 * index_map.at(i) + j) = normal[j];
    d(row) = cpts[i] * normal;
  }

  size_t all_rows = nvars * 3 + n_tconstr;
  Eigen::MatrixXd AC(all_rows, all_rows);
  Eigen::VectorXd bd(all_rows);
  AC << A.transpose() * A, C.transpose(), C, Eigen::MatrixXd::Zero(n_tconstr, n_tconstr);
  bd << A.transpose() * b, d;
  Eigen::VectorXd x = AC.colPivHouseholderQr().solve(bd);

  for (auto [i, k] : index_map)
    cpts[i] = Point3D(x(3 * k), x(3 * k + 1), x(3 * k + 2));
}

void bsplineFit(BSSurface &surface, size_t resolution, const PointVector &samples,
                const std::function<MoveConstraint(size_t,size_t)> &constraint,
                double smoothness) {
  size_t pu = surface.basisU().degree(), pv = surface.basisV().degree();
  auto [nu, nv] = surface.numControlPoints();
  auto &cpts = surface.controlPoints();

  std::map<IndexPair, size_t> index_map;
  size_t index = 0;
  for (size_t i = 0; i < nu; ++i)
    for (size_t j = 0; j < nv; ++j)
      if (std::holds_alternative<MoveType::Free>(constraint(i, j))) {
        index_map[{i,j}] = index;
        index++;
      }
  size_t nvars = index;

  size_t n_smoothing = (nu - 1) * (nv - 1);
  Eigen::MatrixXd A(samples.size() + n_smoothing, nvars);
  Eigen::MatrixXd b(samples.size() + n_smoothing, 3);

  A.setZero(); b.setZero();

  auto addValue = [&](size_t row, size_t i, size_t j, double x) {
    if (std::holds_alternative<MoveType::Fixed>(constraint(i, j)))
      b.row(row) -= VecMap(cpts[i*nv+j].data()) * x;
    else
      A(row, index_map.at({i,j})) = x;
  };

  for (size_t i = 0, index = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    size_t span_u = surface.basisU().findSpan(u);
    DoubleVector coeff_u;
    surface.basisU().basisFunctions(span_u, u, coeff_u);
    for (size_t j = 0; j <= resolution; ++j, ++index) {
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
      if (std::holds_alternative<MoveType::Free>(constraint(i, j))) {
        size_t k = index_map.at({i,j});
        cpts[index] = Point3D(x(k, 0), x(k, 1), x(k, 2));
      }
      index++;
    }
}
