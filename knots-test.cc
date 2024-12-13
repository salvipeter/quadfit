#include <algorithm>
#include <iostream>

#include "knots.hh"

std::vector<double> readKnots() {
  std::vector<double> result;
  size_t n;
  std::cin >> n;
  for (size_t i = 0; i < n; ++i) {
    double x;
    std::cin >> x;
    result.push_back(x);
  }
  return result;
}

void printKnots(std::string name, const std::vector<double> &knots) {
  std::cout << name << ":";
  for (auto k : knots)
    std::cout << ' ' << k;
  std::cout << std::endl;
}

int main(int argc, char **argv) {
  std::cin.exceptions(std::ios::failbit | std::ios::badbit);
  while (true) {
    try {
      size_t p, n;
      double tol;
      std::cin >> p >> n >> tol;
      auto a = readKnots();
      auto b = readKnots();
      auto knots = generateKnots(a, b, p, n, tol);
      std::cout << "**********" << std::endl;
      std::cout << "Degree: " << p << std::endl;
      std::cout << "Inner knots: " << n << std::endl;
      std::cout << "Tolerance: " << tol << std::endl;
      printKnots("A", a);
      printKnots("B", b);
      printKnots("Knots", knots);
      double smallest = 1;
      for (size_t i = 1; i < knots.size(); ++i)
        if (knots[i-1] != knots[i] && knots[i] - knots[i-1] < smallest)
          smallest = knots[i] - knots[i-1];
      std::cout << "-- Smallest interval: " << smallest << std::endl;
      size_t A = 0, B = 0;
      for (auto k : knots)
        if (k != 0 && k != 1) {
          if (std::find(a.begin(), a.end(), k) != a.end())
            A++;
          if (std::find(b.begin(), b.end(), k) != b.end())
            B++;
        }
      std::cout << "-- Same knots (A/B): " << A << '/' << B << std::endl;
    } catch(...) {
      break;
    }
  }
}
