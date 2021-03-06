#include <algorithm>
#include <iterator>

#include "connect-g1.hh"
#include "multiply.hh"

using namespace Geometry;

// coordinates of u when decomposed in the (v,w) system
// (Vector3D instead of Vector2D for convenience)
static Vector3D inSystem(const Vector3D &u, const Vector3D &v, const Vector3D &w) {
  auto uv = u ^ v, uw = u ^ w, vw = v ^ w;
  return {
    uw.norm() / vw.norm() * (uw * vw > 0 ? 1 : -1),
    uv.norm() / vw.norm() * (uv * vw < 0 ? 1 : -1),
    0
  };
}

BSSurface connectG1(const BSCurve &c, const BSCurve &c1, const BSCurve &c2) {
  Point3D S_0, S_1, S1_0, S1_1, S2_0, S2_1;
  Vector3D Su_0, Su_1, Suu_0, Suu_1;
  Vector3D S1v_0, S2v_0, S1v_1, S2v_1;
  Vector3D S1uv_0, S2uv_0, S1uv_1, S2uv_1;

  double lo = c.basis().low(), hi = c.basis().high();

  VectorVector der;
  S_0 = c.eval(lo, 2, der);
  Su_0 = der[1];
  Suu_0 = der[2];
  S_1 = c.eval(hi, 2, der);
  Su_1 = der[1];
  Suu_1 = der[2];
  S1_0 = c1.eval(lo, 1, der);
  S1v_0 = S1_0 - S_0;
  S1uv_0 = der[1] - Su_0;
  S2_0 = c2.eval(lo, 1, der);
  S2v_0 = S2_0 - S_0;
  S2uv_0 = der[1] - Su_0;
  S1_1 = c1.eval(hi, 1, der);
  S1v_1 = S1_1 - S_1;
  S1uv_1 = der[1] - Su_1;
  S2_1 = c2.eval(hi, 1, der);
  S2v_1 = S2_1 - S_1;
  S2uv_1 = der[1] - Su_1;

  // Sv = D * alpha + Su * beta
  auto D0 = (S1v_0 - S2v_0) / 2;
  auto D1 = (S1v_1 - S2v_1) / 2;
  auto ab_0 = inSystem(S1v_0, D0, Su_0);
  auto ab_1 = inSystem(S1v_1, D1, Su_1);
  auto ab_0d = inSystem(S1uv_0 - (S1uv_0 - S2uv_0) / 2 * ab_0[0] - Suu_0 * ab_0[1], D0, Su_0);
  auto ab_1d = inSystem(S1uv_1 - (S1uv_1 - S2uv_1) / 2 * ab_1[0] - Suu_1 * ab_1[1], D1, Su_1);
  BSCurve scaling(3, { lo, lo, lo, lo, hi, hi, hi, hi },
                  { ab_0, ab_0 + ab_0d / 3, ab_1 - ab_1d / 3, ab_1 });

  auto point = [&](double u) {
    return c.eval(u);
  };
  auto derivative = [&](double u) {
    VectorVector der;
    c.eval(u, 1, der);
    auto p1 = c1.eval(u);
    auto p2 = c2.eval(u);
    auto s = scaling.eval(u);
    return (p1 - p2) / 2 * s[0] + der[1] * s[1];
  };
  return multiplyBSplines(c.basis(), scaling.basis(), point, derivative);
}
