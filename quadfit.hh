#pragma once

#include <set>

#include <geometry.hh>

class QuadFit {
public:
  std::string readPWGB(std::string filename);      // returns model description
  std::vector<Geometry::BSSurface> fit();

private:
  using Ribbon = std::array<Geometry::BSCurve, 2>; // outer, inner
  struct Quad {
    struct Boundary {
      size_t segment;
      bool reversed;
      bool on_ribbon;
      size_t ribbon;         // when on_ribbon
      double s0, s1, h0, h1; //   is true
      Geometry::BSSurface sextic;
    };
    std::array<Boundary, 4> boundaries; // u=0, v=0, u=1, v=1
    size_t resolution;
    Geometry::PointVector samples;
  };

  std::vector<Geometry::BSSurface> initialFit() const;

  std::vector<Ribbon> ribbons;
  std::vector<Geometry::BSCurve> segments;
  std::vector<Quad> quads;
};
