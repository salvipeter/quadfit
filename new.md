# New concept

We have two modes: approximative or "numerical" G1 (NG1),
and a variation where we have exact G1 at the outer boundaries.

## NG1

1. Subdivide the boundary curves and insert them into the respective quad patches.

1. Add extra knots (set by `cubic-fit`) and fit with all the remaining control points.

1. Ensure exact C0 continuity along the inner boundaries by taking the mean 
   of the control points of the adjacent quads.

This gives exact C0 and approximative G1 everywhere with cubic patches.

```
./quadfit-test examples/test.pwgb \
  --cubic-fit=5 \
  --fix-c0-inside --fix-c0-outside \
  --print-continuity-errors --print-approximation-errors
```

## (N)G1

1. Scale and subdivide the ribbons as before (Step 2 in [the notes](notes.md)),
   and insert them in sextic quads.

1. Create interior boundary curves between MAT vertices and side vertices,
   approximating sampled points and also interpolating the outer ribbon's
   cross derivative. For this we use quartic curves of *2k+2* segments,
   where *k* is set by `extra-knots`. These are then degree-elevated and
   also inserted into the respective quads (inserting knots into the quads
   as needed).

1. Fit the sampled points with all unset control points of the quad.

This gives exact G1 interpolation of the ribbons, and approximative G1
at the interior boundaries with sextic patches.

```
./quadfit-test examples/test.pwgb \
  --fit-curves --extra-knots=1 \
  --retain-boundaries --retain-ribbons \
  --print-continuity-errors --print-approximation-errors
```
