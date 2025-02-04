#include <map>

#include <Eigen/SparseQR>
#include <Eigen/QR>

#include "bspline-fit.hh"
#include "closest-point.hh"

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
                double deviation_penalty, double smoothness) {
  const auto &basis = curve.basis();
  size_t p = basis.degree();
  auto &cpts = curve.controlPoints();
  size_t n = cpts.size();
  size_t resolution = samples.size() - 1;

  std::map<size_t, size_t> index_map, tangent_map, normal_map;
  size_t nvars = 0, n_tconstr = 0, n_nconstr = 0;
  for (size_t i = 0; i < n; ++i)
    if (!std::holds_alternative<MoveType::Fixed>(constraint(i))) {
      index_map[i] = nvars++;
      if (std::holds_alternative<MoveType::Tangent>(constraint(i)))
        tangent_map[i] = n_tconstr++;
      if (std::holds_alternative<MoveType::Normal>(constraint(i)))
        normal_map[i] = n_nconstr++;
    }

  size_t n_deviation = deviation_penalty > 0 ? nvars : 0;
  size_t n_smoothing = smoothness == 0 ? 0 : n - 1;
  size_t n_rows = samples.size() + n_smoothing;
  size_t n_constr = n_tconstr + n_nconstr * 3;
  Eigen::MatrixXd A(n_rows * 3 + n_deviation * 3, nvars * 3), C(n_constr, nvars * 3);
  Eigen::VectorXd b(n_rows * 3 + n_deviation * 3), d(n_constr);

  A.setZero(); C.setZero();

  auto addValue = [&](size_t row, size_t i, double x) {
    if (std::holds_alternative<MoveType::Fixed>(constraint(i)))
      for (size_t j = 0; j < 3; ++j)
        b(row + j) -= cpts[i][j] * x;
    else
      for (size_t j = 0; j < 3; ++j)
        A(row + j, 3 * index_map.at(i) + j) = x;
  };

  for (size_t i = 0; i <= resolution; ++i) {
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

  if (smoothness != 0) {
    size_t row = resolution + 1;
    for (size_t i = 1; i < n - 1; ++i) {
      addValue(3 * row, i, smoothness);
      addValue(3 * row, i - 1, -0.5 * smoothness);
      addValue(3 * row, i + 1, -0.5 * smoothness);
      row++;
    }
  }

  if (deviation_penalty > 0)
    for (size_t i = 0; i < n; ++i)
      if (!std::holds_alternative<MoveType::Fixed>(constraint(i))) {
        size_t k = index_map.at(i);
        for (size_t j = 0; j < 3; ++j) {
          A(n_rows * 3 + k * 3 + j, k * 3 + j) = deviation_penalty;
          b(n_rows * 3 + k * 3 + j) = cpts[i][j] * deviation_penalty;
        }
      }

  for (auto [i, row] : tangent_map) {
    auto normal = std::get<MoveType::Tangent>(constraint(i)).normal;
    for (size_t j = 0; j < 3; ++j)
      C(row, 3 * index_map.at(i) + j) = normal[j];
    d(row) = cpts[i] * normal;
  }
  for (auto [i, r] : normal_map) {
    size_t row = n_tconstr + r * 3;
    auto normal = std::get<MoveType::Normal>(constraint(i)).normal;
    C(row + 0, 3 * index_map.at(i) + 1) = normal[2];
    C(row + 0, 3 * index_map.at(i) + 2) = -normal[1];
    d(row + 0) = cpts[i][1] * normal[2] - cpts[i][2] * normal[1];
    C(row + 1, 3 * index_map.at(i) + 0) = normal[2];
    C(row + 1, 3 * index_map.at(i) + 2) = -normal[0];
    d(row + 1) = cpts[i][0] * normal[2] - cpts[i][2] * normal[0];
    C(row + 2, 3 * index_map.at(i) + 0) = normal[1];
    C(row + 2, 3 * index_map.at(i) + 1) = -normal[0];
    d(row + 2) = cpts[i][0] * normal[1] - cpts[i][1] * normal[0];
  }

  size_t all_rows = nvars * 3 + n_constr;
  Eigen::MatrixXd AC(all_rows, all_rows);
  Eigen::VectorXd bd(all_rows);
  AC << A.transpose() * A, C.transpose(), C, Eigen::MatrixXd::Zero(n_constr, n_constr);
  bd << A.transpose() * b, d;
  Eigen::VectorXd x = AC.colPivHouseholderQr().solve(bd);

  for (auto [i, k] : index_map)
    cpts[i] = Point3D(x(3 * k), x(3 * k + 1), x(3 * k + 2));
}

void bsplineFit(BSCurve &curve, const PointVector &points, const VectorVector &normals,
                const std::function<MoveConstraint(size_t)> &constraint,
                double deviation_penalty, double smoothness) {
  const auto &basis = curve.basis();
  size_t p = basis.degree();
  auto &cpts = curve.controlPoints();
  size_t n = cpts.size();
  size_t resolution = points.size() - 1;

  std::map<size_t, size_t> index_map, tangent_map, normal_map;
  size_t nvars = 0, n_tconstr = 0, n_nconstr = 0;
  for (size_t i = 0; i < n; ++i)
    if (!std::holds_alternative<MoveType::Fixed>(constraint(i))) {
      index_map[i] = nvars++;
      if (std::holds_alternative<MoveType::Tangent>(constraint(i)))
        tangent_map[i] = n_tconstr++;
      if (std::holds_alternative<MoveType::Normal>(constraint(i)))
        normal_map[i] = n_nconstr++;
    }

  size_t n_deviation = deviation_penalty > 0 ? nvars : 0;
  size_t n_smoothing = smoothness == 0 ? 0 : n - 1;
  size_t n_rows = points.size();
  size_t n_constr = n_tconstr + n_nconstr * 3;
  Eigen::MatrixXd A(n_rows * 4 + n_smoothing * 3 + n_deviation * 3, nvars * 3),
    C(n_constr, nvars * 3);
  Eigen::VectorXd b(n_rows * 4 + n_smoothing * 3 + n_deviation * 3), d(n_constr);

  A.setZero(); C.setZero(); b.setZero();

  auto addValue = [&](size_t row, size_t i, double x) {
    if (std::holds_alternative<MoveType::Fixed>(constraint(i)))
      for (size_t j = 0; j < 3; ++j)
        b(row + j) -= cpts[i][j] * x;
    else
      for (size_t j = 0; j < 3; ++j)
        A(row + j, 3 * index_map.at(i) + j) = x;
  };

  auto addVector = [&](size_t row, size_t i, const Vector3D &x) {
    if (std::holds_alternative<MoveType::Fixed>(constraint(i)))
      for (size_t j = 0; j < 3; ++j)
        b(row) -= cpts[i][j] * x[j];
    else
      for (size_t j = 0; j < 3; ++j)
        A(row, 3 * index_map.at(i) + j) = x[j];
  };

  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    size_t span = basis.findSpan(u);
    DoubleMatrix coeff;
    basis.basisFunctionDerivatives(span, u, 1, coeff);
    for (size_t j = 0; j < 3; ++j)
      b(4 * i + j) = points[i][j];
    for (size_t k = 0; k <= p; ++k) {
      size_t i1 = span - p + k;
      addValue(4 * i, i1, coeff[0][k]);
      addVector(4 * i + 3, i1, normals[i] * coeff[1][k]);
    }
  }

  if (smoothness != 0) {
    size_t row = resolution + 1;
    for (size_t i = 1; i < n - 1; ++i) {
      addValue(4 * row, i, smoothness);
      addValue(4 * row, i - 1, -0.5 * smoothness);
      addValue(4 * row, i + 1, -0.5 * smoothness);
      row++;
    }
  }

  if (deviation_penalty > 0)
    for (size_t i = 0; i < n; ++i)
      if (!std::holds_alternative<MoveType::Fixed>(constraint(i))) {
        size_t k = index_map.at(i);
        for (size_t j = 0; j < 3; ++j) {
          A(n_rows * 4 + n_smoothing * 3 + k * 3 + j, k * 3 + j) = deviation_penalty;
          b(n_rows * 4 + n_smoothing * 3 + k * 3 + j) = cpts[i][j] * deviation_penalty;
        }
      }

  for (auto [i, row] : tangent_map) {
    auto normal = std::get<MoveType::Tangent>(constraint(i)).normal;
    for (size_t j = 0; j < 3; ++j)
      C(row, 3 * index_map.at(i) + j) = normal[j];
    d(row) = cpts[i] * normal;
  }
  for (auto [i, r] : normal_map) {
    size_t row = n_tconstr + r * 3;
    auto normal = std::get<MoveType::Normal>(constraint(i)).normal;
    C(row + 0, 3 * index_map.at(i) + 1) = normal[2];
    C(row + 0, 3 * index_map.at(i) + 2) = -normal[1];
    d(row + 0) = cpts[i][1] * normal[2] - cpts[i][2] * normal[1];
    C(row + 1, 3 * index_map.at(i) + 0) = normal[2];
    C(row + 1, 3 * index_map.at(i) + 2) = -normal[0];
    d(row + 1) = cpts[i][0] * normal[2] - cpts[i][2] * normal[0];
    C(row + 2, 3 * index_map.at(i) + 0) = normal[1];
    C(row + 2, 3 * index_map.at(i) + 1) = -normal[0];
    d(row + 2) = cpts[i][0] * normal[1] - cpts[i][1] * normal[0];
  }

  size_t all_rows = nvars * 3 + n_constr;
  Eigen::MatrixXd AC(all_rows, all_rows);
  Eigen::VectorXd bd(all_rows);
  AC << A.transpose() * A, C.transpose(), C, Eigen::MatrixXd::Zero(n_constr, n_constr);
  bd << A.transpose() * b, d;
  Eigen::VectorXd x = AC.colPivHouseholderQr().solve(bd);

  for (auto [i, k] : index_map)
    cpts[i] = Point3D(x(3 * k), x(3 * k + 1), x(3 * k + 2));
}

void bsplineFit(BSCurve &curve, const PointVector &points, const PointVector &normals,
                const std::function<MoveConstraint(size_t)> &constraint,
                double deviation_penalty) {
  const auto &basis = curve.basis();
  size_t p = basis.degree();
  auto &cpts = curve.controlPoints();
  size_t n = cpts.size();
  size_t resolution = points.size() - 1;

  std::map<size_t, size_t> index_map, tangent_map, normal_map;
  size_t nvars = 0, n_tconstr = 0, n_nconstr = 0;
  for (size_t i = 0; i < n; ++i)
    if (!std::holds_alternative<MoveType::Fixed>(constraint(i))) {
      index_map[i] = nvars++;
      if (std::holds_alternative<MoveType::Tangent>(constraint(i)))
        tangent_map[i] = n_tconstr++;
      if (std::holds_alternative<MoveType::Normal>(constraint(i)))
        normal_map[i] = n_nconstr++;
    }

  size_t n_deviation = deviation_penalty > 0 ? nvars : 0;
  size_t n_rows = points.size();
  size_t n_constr = n_tconstr + n_nconstr * 3;
  Eigen::MatrixXd A(n_rows + n_deviation * 3, nvars * 3), C(n_constr, nvars * 3);
  Eigen::VectorXd b(n_rows + n_deviation * 3), d(n_constr);

  A.setZero(); C.setZero();

  auto addValue = [&](size_t row, size_t i, const Vector3D &x) {
    if (std::holds_alternative<MoveType::Fixed>(constraint(i)))
      for (size_t j = 0; j < 3; ++j)
        b(row) -= cpts[i][j] * x[j];
    else
      for (size_t j = 0; j < 3; ++j)
        A(row, 3 * index_map.at(i) + j) = x[j];
  };

  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    size_t span = basis.findSpan(u);
    DoubleVector coeff;
    basis.basisFunctions(span, u, coeff);
    b(i) = points[i] * normals[i];
    for (size_t k = 0; k <= p; ++k) {
      size_t i1 = span - p + k;
      addValue(i, i1, normals[i] * coeff[k]);
    }
  }

  if (deviation_penalty > 0)
    for (size_t i = 0; i < n; ++i)
      if (!std::holds_alternative<MoveType::Fixed>(constraint(i))) {
        size_t k = index_map.at(i);
        for (size_t j = 0; j < 3; ++j) {
          A(n_rows + k * 3 + j, k * 3 + j) = deviation_penalty;
          b(n_rows + k * 3 + j) = cpts[i][j] * deviation_penalty;
        }
      }

  for (auto [i, row] : tangent_map) {
    auto normal = std::get<MoveType::Tangent>(constraint(i)).normal;
    for (size_t j = 0; j < 3; ++j)
      C(row, 3 * index_map.at(i) + j) = normal[j];
    d(row) = cpts[i] * normal;
  }
  for (auto [i, r] : normal_map) {
    size_t row = n_tconstr + r * 3;
    auto normal = std::get<MoveType::Normal>(constraint(i)).normal;
    C(row + 0, 3 * index_map.at(i) + 1) = normal[2];
    C(row + 0, 3 * index_map.at(i) + 2) = -normal[1];
    d(row + 0) = cpts[i][1] * normal[2] - cpts[i][2] * normal[1];
    C(row + 1, 3 * index_map.at(i) + 0) = normal[2];
    C(row + 1, 3 * index_map.at(i) + 2) = -normal[0];
    d(row + 1) = cpts[i][0] * normal[2] - cpts[i][2] * normal[0];
    C(row + 2, 3 * index_map.at(i) + 0) = normal[1];
    C(row + 2, 3 * index_map.at(i) + 1) = -normal[0];
    d(row + 2) = cpts[i][0] * normal[1] - cpts[i][1] * normal[0];
  }

  size_t all_rows = nvars * 3 + n_constr;
  Eigen::MatrixXd AC(all_rows, all_rows);
  Eigen::VectorXd bd(all_rows);
  AC << A.transpose() * A, C.transpose(), C, Eigen::MatrixXd::Zero(n_constr, n_constr);
  bd << A.transpose() * b, d;
  Eigen::VectorXd x = AC.colPivHouseholderQr().solve(bd);

  for (auto [i, k] : index_map)
    cpts[i] = Point3D(x(3 * k), x(3 * k + 1), x(3 * k + 2));
}

void bsplineFit(BSSurface &surface, size_t resolution, const PointVector &samples,
                const std::function<MoveConstraint(size_t,size_t)> &constraint,
                double deviation_penalty, bool parameter_correction) {
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

  Eigen::SparseMatrix<double> A(samples.size() + nvars, nvars);
  Eigen::MatrixXd b(samples.size() + nvars, 3);
  std::vector<Eigen::Triplet<double>> triplets;

  triplets.reserve(samples.size() * (pu + 1) * (pv + 1) + nu * nv);
  b.setZero();

  auto addValue = [&](size_t row, size_t i, size_t j, double x) {
    if (std::holds_alternative<MoveType::Fixed>(constraint(i, j)))
      b.row(row) -= VecMap(cpts[i*nv+j].data()) * x;
    else
      triplets.emplace_back(row, index_map.at({i,j}), x);
  };

  index = 0;
  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    for (size_t j = 0; j <= resolution; ++j, ++index) {
      double v = (double)j / resolution;

      if (parameter_correction)
        closestPoint(surface, samples[index], u, v, 20, 0, 0);

      size_t span_u = surface.basisU().findSpan(u);
      size_t span_v = surface.basisV().findSpan(v);
      DoubleVector coeff_u, coeff_v;
      surface.basisU().basisFunctions(span_u, u, coeff_u);
      surface.basisV().basisFunctions(span_v, v, coeff_v);

      b.row(index) = VecMap(samples[index].data());
      for (size_t k = 0; k <= pu; ++k)
        for (size_t l = 0; l <= pv; ++l) {
          size_t i1 = span_u - pu + k, j1 = span_v - pv + l;
          addValue(index, i1, j1, coeff_u[k] * coeff_v[l]);
        }
    }
  }

  // Deviation terms
  if (deviation_penalty > 0)
    for (size_t i = 0; i < nu; ++i)
      for (size_t j = 0; j < nv; ++j)
        if (std::holds_alternative<MoveType::Free>(constraint(i, j))) {
          b.row(index) = VecMap(cpts[i*nv+j].data()) * deviation_penalty;
          triplets.emplace_back(index, index_map.at({i,j}), deviation_penalty);
          index++;
        }

  A.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver(A);
  Eigen::MatrixXd x = solver.solve(b);

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
