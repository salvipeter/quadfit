#include "closest-point.hh"

using namespace Geometry;

Point3D closestPoint(const BSCurve &curve, const Point3D &point, double &u,
                     size_t max_iteration, double distance_tol, double cosine_tol) {
  VectorVector der;
  auto lo = curve.basis().low();
  auto hi = curve.basis().high();

  for (size_t iteration = 0; iteration < max_iteration; ++iteration) {
    auto deviation = curve.eval(u, 2, der) - point;
    auto distance = deviation.norm();
    if (distance < distance_tol)
      break;

    double scaled_error = der[1] * deviation;
    double cosine_err = std::abs(scaled_error) / (der[1].norm() * distance);
    if (cosine_err < cosine_tol)
      break;

    double old = u;
    u -= scaled_error / (der[2] * deviation + der[1] * der[1]);
    u = std::min(std::max(u, lo), hi);

    if ((der[1] * (u - old)).norm() < distance_tol)
      break;
  }

  return der[0];
}

Point3D closestPoint(const BSSurface &surface, const Point3D &point, double &u, double &v,
                     size_t max_iteration, double distance_tol, double cosine_tol) {
  VectorMatrix der;
  auto u_lo = surface.basisU().low(), u_hi = surface.basisU().high();
  auto v_lo = surface.basisV().low(), v_hi = surface.basisV().high();

  for (size_t iteration = 0; iteration < max_iteration; ++iteration) {
    auto deviation = surface.eval(u, v, 2, der) - point;
    auto distance = deviation.norm();
    if (distance < distance_tol)
      break;

    auto rSu = deviation * der[1][0];
    auto rSv = deviation * der[0][1];
    if (std::abs(rSu) / (distance * der[1][0].norm()) < cosine_tol ||
        std::abs(rSv) / (distance * der[0][1].norm()) < cosine_tol)
      break;

    double Ja = der[1][0].normSqr() + deviation * der[2][0];
    double Jb = der[1][0] * der[0][1] + deviation * der[1][1];
    double Jd = der[0][1].normSqr() + deviation * der[0][2];
    double D = Ja * Jd - Jb * Jb;
    double du = (Jd * rSu - Jb * rSv) / D;
    double dv = (Ja * rSv - Jb * rSu) / D;

    u = std::min(std::max(u - du, u_lo), u_hi);
    v = std::min(std::max(v - dv, v_lo), v_hi);

    if ((der[1][0] * du + der[0][1] * dv).norm() < distance_tol)
      break;
  }

  return der[0][0];
}
