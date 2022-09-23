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

using VecMap = Eigen::Map<const Eigen::Vector3d>;

void applyMask(Geometry::BSSurface &surface, DiscreteMask type) {
  const auto &mask = masks[static_cast<size_t>(type)];
  auto [n, m] = surface.numControlPoints();
  size_t vars = (n - 4) * (m - 4);
  Eigen::MatrixXd A(vars, vars); A.setZero();
  Eigen::MatrixXd b(vars, 3); b.setZero();
  for (size_t i = 2, k = 0; i < n - 2; ++i)
    for (size_t j = 2; j < m - 2; ++j, ++k)
      for (size_t di = 0; di < 5; ++di)
        for (size_t dj = 0; dj < 5; ++dj) {
          double w = mask[di][dj];
          if (w == 0)
            continue;
          size_t i1 = i + di - 2, j1 = j + dj - 2;
          if (i1 < 2 || i1 >= n - 2 || j1 < 2 || j1 >= m - 2)
            b.row(k) -= VecMap(surface.controlPoint(i1, j1).data()) * w;
          else
            A(k, (i1 - 2) * (m - 4) + j1 - 2) = w;
        }
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for (size_t i = 2, k = 0; i < n - 2; ++i)
    for (size_t j = 2; j < m - 2; ++j, ++k)
      surface.controlPoint(i, j) = Point3D(x(k, 0), x(k, 1), x(k, 2));
}
