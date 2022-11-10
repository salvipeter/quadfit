#include <Eigen/Dense>

#include "discrete-mask.hh"

using namespace Geometry;

using MaskType = std::array<std::array<double, 5>, 5>;
static constexpr std::array<MaskType, 4> masks = {
  { { { {  0,  0,  0,  0,  0 },
        {  0, -1,  2, -1,  0 },
        {  0,  2, -4,  2,  0 },
        {  0, -1,  2, -1,  0 },
        {  0,  0,  0,  0,  0 } } },

    { { { -1,  4, -6,  4, -1 },
        {  4,-16, 24,-16,  4 },
        { -6, 24,-36, 24, -6 },
        {  4,-16, 24,-16,  4 },
        { -1,  4, -6,  4, -1 } } },

    { { {  0,  0,  0,  0,  0 },
        {  0,  0,  1,  0,  0 },
        {  0,  1, -4,  1,  0 },
        {  0,  0,  1,  0,  0 },
        {  0,  0,  0,  0,  0 } } },

    { { {  0,  0,  1,  0,  0 },
        {  0,  2, -8,  2,  0 },
        {  1, -8, 20, -8,  1 },
        {  0,  2, -8,  2,  0 },
        {  0,  0,  1,  0,  0 } } } }
};

// also of the previous and next 2 control points
static DoubleVector grevilleAbscissae(const BSBasis &basis, size_t index) {
  size_t d = basis.degree();
  const auto &knots = basis.knots();
  DoubleVector result;
  for (size_t i = 0; i < 5; ++i) {
    double sum = 0;
    for (size_t j = 1; j <= d; ++j)
      sum += knots[index + i - 2 + j];
    result.push_back(sum / d);
  }
  return result;
}

static Eigen::MatrixXd setupMask(DiscreteMask type, const BSSurface &s, size_t ui, size_t vi) {
  Eigen::MatrixXd result(5, 5);
  if (type == DiscreteMask::C1_COONS) {
    auto gr_u = grevilleAbscissae(s.basisU(), ui);
    auto gr_v = grevilleAbscissae(s.basisV(), vi);
    DoubleVector du, dv;
    for (size_t i = 1; i < 5; ++i) {
      du.push_back(gr_u[i] - gr_u[i-1]);
      dv.push_back(gr_v[i] - gr_v[i-1]);
    }
    auto wu = [&](size_t i, size_t j, size_t k) { return 1 / (du[i] * du[j] * du[k]); };
    auto wv = [&](size_t i, size_t j, size_t k) { return 1 / (dv[i] * dv[j] * dv[k]); };
    DoubleVector u = {
      wu(0,1,2),
      -(wu(0,1,2) + wu(1,1,1) + wu(1,1,2) + wu(1,2,2)),
      wu(1,1,1) + 2 * wu(1,1,2) + 2 * wu(1,2,2) + wu(2,2,2),
      -(wu(3,2,1) + wu(2,2,2) + wu(2,2,1) + wu(2,1,1)),
      wu(3,2,1)
    };
    DoubleVector v = {
      wv(0,1,2),
      -(wv(0,1,2) + wv(1,1,1) + wv(1,1,2) + wv(1,2,2)),
      wv(1,1,1) + 2 * wv(1,1,2) + 2 * wv(1,2,2) + wv(2,2,2),
      -(wv(3,2,1) + wv(2,2,2) + wv(2,2,1) + wv(2,1,1)),
      wv(3,2,1)
    };
    // auto wu = [&](size_t i, size_t j, size_t k) {
    //   return 1 / (du[i] * (du[i] + du[j]) * (du[i] + du[j] + du[k]));
    // };
    // auto wv = [&](size_t i, size_t j, size_t k) {
    //   return 1 / (dv[i] * (dv[i] + dv[j]) * (dv[i] + dv[j] + dv[k]));
    // };
    // DoubleVector u = {
    //   wu(0,1,2),
    //   -(wu(0,1,2) + wu(1,0,2) + wu(1,2,0) + wu(1,2,3)),
    //   wu(1,0,2) + wu(1,2,0) + wu(2,1,0) + wu(1,2,3) + wu(2,1,3) + wu(2,3,1),
    //   -(wu(2,1,0) + wu(2,1,3) + wu(2,3,1) + wu(3,2,1)),
    //   wu(3,2,1)
    // };
    // DoubleVector v = {
    //   wv(0,1,2),
    //   -(wv(0,1,2) + wv(1,0,2) + wv(1,2,0) + wv(1,2,3)),
    //   wv(1,0,2) + wv(1,2,0) + wv(2,1,0) + wv(1,2,3) + wv(2,1,3) + wv(2,3,1),
    //   -(wv(2,1,0) + wv(2,1,3) + wv(2,3,1) + wv(3,2,1)),
    //   wv(3,2,1)
    // };
    for (size_t i = 0; i < 5; ++i)
      for (size_t j = 0; j < 5; ++j)
        result(i, j) = u[i] * v[j];
    return result;
  }
  const auto &mask = masks[static_cast<size_t>(type)];
  double sum = 0;
  for (size_t i = 0; i < 5; ++i)
    for (size_t j = 0; j < 5; ++j) {
      if (i == 2 && j == 2)
        continue;
      double w = mask[i][j];
      result(i, j) = w;
      sum += w;
    }
  result(2, 2) = -sum;
  return result;
}

using VecMap = Eigen::Map<const Eigen::Vector3d>;

void applyMask(Geometry::BSSurface &surface, DiscreteMask type) {
  auto [n, m] = surface.numControlPoints();
  size_t vars = (n - 4) * (m - 4);
  Eigen::MatrixXd A(vars, vars); A.setZero();
  Eigen::MatrixXd b(vars, 3); b.setZero();
  for (size_t i = 2, k = 0; i < n - 2; ++i)
    for (size_t j = 2; j < m - 2; ++j, ++k) {
      auto mask = setupMask(type, surface, i, j);
      for (size_t di = 0; di < 5; ++di)
        for (size_t dj = 0; dj < 5; ++dj) {
          double w = mask(di, dj);
          if (w == 0)
            continue;
          size_t i1 = i + di - 2, j1 = j + dj - 2;
          if (i1 < 2 || i1 >= n - 2 || j1 < 2 || j1 >= m - 2)
            b.row(k) -= VecMap(surface.controlPoint(i1, j1).data()) * w;
          else
            A(k, (i1 - 2) * (m - 4) + j1 - 2) = w;
        }
    }
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for (size_t i = 2, k = 0; i < n - 2; ++i)
    for (size_t j = 2; j < m - 2; ++j, ++k)
      surface.controlPoint(i, j) = Point3D(x(k, 0), x(k, 1), x(k, 2));
}
