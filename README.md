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
