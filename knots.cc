#include <algorithm>

#include "knots.hh"

std::vector<double> generateKnots(const std::vector<double> &a,
                                  const std::vector<double> &b,
                                  size_t p, size_t n, double tolerance) {
  std::vector<double> result;

  // 1. Union of `a` and `b`
  {
    std::fill_n(std::back_inserter(result), p + 1, 0);
    size_t i = 0, j = 0;
    while (a[i] < 1 || b[j] < 1) {
      if (a[i] < b[j]) {
        if (a[i] != result.back())
          result.push_back(a[i]);
        i++;
      } else {
        if (b[j] != result.back())
          result.push_back(b[j]);
        j++;
      }
    }
    std::fill_n(std::back_inserter(result), p + 1, 1);
  }

  // 2. Eliminate close knots
  while (true) {
    // Find closest two knots
    double min = 1;
    size_t mini = 0;
    for (size_t i = 1; i < result.size(); ++i)
      if (result[i-1] != result[i] && result[i] - result[i-1] < min) {
        min = result[i] - result[i-1];
        mini = i - 1;
      }
    if (min > tolerance)
      break;
    if (result[mini] == 0 ||
        (result[mini] - result[mini-1] < result[mini+2] - result[mini+1] &&
         result[mini+1] != 1))
      result.erase(result.begin() + mini + 1);
    else
      result.erase(result.begin() + mini);
  }

  // 3. Add more knots when needed
  while (result.size() - 2 * (p + 1) < n) {
    // Find largest interval
    double max = 0;
    size_t maxi = 0;
    for (size_t i = 1; i < result.size(); ++i)
      if (result[i] - result[i-1] > max) {
        max = result[i] - result[i-1];
        maxi = i;
      }
    result.insert(result.begin() + maxi, (result[maxi-1] + result[maxi]) / 2);
  }

  // 4. Delete some knots when needed
  while (result.size() - 2 * (p + 1) > n) {
    // Find shortest interval around a knot
    double min = 1;
    size_t mini = 0;
    for (size_t i = 1; i < result.size() - 1; ++i)
      if (result[i-1] != result[i+1] && result[i+1] - result[i-1] < min) {
        min = result[i+1] - result[i-1];
        mini = i;
      }
    result.erase(result.begin() + mini);
  }

  return result;
}
