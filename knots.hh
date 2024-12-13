#pragma once

#include <vector>

// Generates a clamped knot vector for a degree-`p` basis
// with `n` interior knots of multiplicity 1.
// - Both `a` and `b` are assumed to be clamped knot vectors in [0,1]
// - No knot in the result is closer than `tolerance`
//   to a knot in either `a` or `b`
// - Tries to insert as many knots from `a` and `b` as possible
std::vector<double> generateKnots(const std::vector<double> &a,
                                  const std::vector<double> &b,
                                  size_t p, size_t n, double tolerance);
