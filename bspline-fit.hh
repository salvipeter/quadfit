#pragma once

#include <functional>
#include <variant>

#include <geometry.hh>

namespace MoveType {
  struct Free {};                                 // [3D] movement in any direction
  struct Tangent { Geometry::Vector3D normal; };  // [2D] movement only perpendicular to the normal
  struct Normal  { Geometry::Vector3D normal; };  // [1D] movement only in the normal direction
  struct Fixed {};                                // [0D] no movement
};
using MoveConstraint = std::variant<MoveType::Free, MoveType::Tangent,
                                    MoveType::Normal, MoveType::Fixed>;

// C(k/n) = Pk
void bsplineFit(Geometry::BSCurve &curve, const Geometry::PointVector &samples,
                const std::function<MoveConstraint(size_t)> &constraint,
                double smoothness);

// C(k/n) * Nk = Pk * Nk
void bsplineFit(Geometry::BSCurve &curve, const Geometry::PointVector &points,
                const Geometry::VectorVector &normals,
                const std::function<MoveConstraint(size_t)> &constraint);

void bsplineFit(Geometry::BSSurface &surface, size_t resolution,
                const Geometry::PointVector &samples,
                const std::function<MoveConstraint(size_t,size_t)> &constraint,
                double smoothness);
