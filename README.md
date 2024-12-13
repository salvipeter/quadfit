# quadfit
Fit continuously connected B-spline surfaces on sampled points with boundary constraints.

Uses my [geometry library](https://github.com/salvipeter/libgeom/),
my [jet-fitting wrapper](https://github.com/salvipeter/jet-wrapper/) on CGAL,
and my [Bézier extraction library](https://github.com/salvipeter/bezier-extractions).

Requires:
- A curve network (given as B-spline curves)
- Uniformly sampled points in each quad
- Positional and cross-derivative constraints for the outer boundary (given as pairs of B-spline curves) [optional]

Currently the input can be given in a `.PWGB` file (see `pwgb-spec.txt`).

The test program generates a single mesh for the surface, and another with the control nets.
The boundary curves of all quads are written to an `.OBJ` file as polylines.
The surfaces are also exported in `.QDS` format (see `qds-spec.txt`).

Command line switches:

| Switch                       | Meaning                                                      |
|------------------------------|--------------------------------------------------------------|
| --mesh=*filename.obj*        | Use the specified mesh for jet-fitting curvatures            |
| --resolution=*n*             | Use the specified resolution for mesh output                 |
| --preliminary-fit-tangents   | Extract tangent length from a preliminary C-1 fit            |
| --preliminary-fit-twists     | Extract twist vectors from a preliminary C-1 fit             |
| --preliminary-fit-normals    | Extract normal vectors from a preliminary C-1 fit            |
| --extra-knots=*n*            | Add *2n* extra knots to the subdividing curves               |
| --fit-curves                 | Fit the quartic curves on sampled points                     |
| --fit-normals                | Fit the cross-derivative curves on sampled normals           |
| --coons-patch                | Create Coons patches and do not fit to interior samples      |
| --cubic-fit=*n*              | Use cubic splines with *n* inner knots to fit                |
| --knot-tolerance=*d*         | Minimum interval between knots (only in cubic fit)           |
| --fix-c0-inside              | Fix C0 continuity at the subdivision curves                  |
| --fix-c0-outside             | Fix C0 continuity at the ribbons                             |
| --retain-boundaries          | Retain quad boundaries                                       |
| --retain-direction-blends    | Retain the shared direction blends between the quads         |
| --retain-ribbons             | Retain the original ribbons (scaled unless the fit is cubic) |
| --print-continuity-errors    | Print C0 and G1 errors between the quads and to the ribbons  |
| --print-approximation-errors | Print C0 and G1 errors to the sampled points and normals     |
