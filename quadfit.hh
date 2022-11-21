#pragma once

#include <optional>

#include <geometry.hh>
#include <jet-wrapper.hh>

class QuadFit {
public:
  std::string readPWGB(std::string filename);      // returns model description
  void update(const std::vector<std::string> &switches);
  std::vector<Geometry::BSSurface> fit(const std::vector<std::string> &switches);

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
    std::array<Geometry::Vector3D, 4> twists; // u=v=0; u=1,v=0; u=v=1; u=0,v=1
    size_t resolution;
    Geometry::PointVector samples;
    Geometry::VectorVector normals;
    mutable std::optional<Geometry::BSSurface> preliminary_fit;
  };

  Geometry::BSSurface preliminaryFit(size_t i) const;
  std::vector<Geometry::BSSurface> initialFit(bool fit_tangents) const;
  std::vector<std::pair<size_t, size_t>> ribbonSegments(size_t i) const;
  void correctFirstDerivatives(Geometry::BSSurface &cubic, size_t quad_index) const;
  void correctCubicTwists(Geometry::BSSurface &cubic, size_t quad_index) const;
  void correctSecondDerivatives(Geometry::BSSurface &quintic, size_t quad_index) const;
  void correctTwists(Geometry::BSSurface &quintic, size_t quad_index) const;
  Geometry::BSSurface innerBoundaryRibbon(const std::vector<Geometry::BSSurface> &quintic_patches,
                                          size_t quad_index, size_t side, size_t extra_knots,
                                          bool prelim_normals, bool fitC0, bool fitG1) const;
  void printContinuityErrors(const std::vector<Geometry::BSSurface> &result) const;
  void printApproximationErrors(const std::vector<Geometry::BSSurface> &result) const;

  std::vector<Ribbon> ribbons;
  std::vector<Geometry::BSCurve> segments;
  std::vector<Quad> quads;

  // Additional topology & geometry structures
  // (updated automatically by update())
  using IndexPair = std::pair<size_t, size_t>;
  using IndexPairVector = std::vector<IndexPair>;
  std::vector<Geometry::Point3D> vertices; // all segment endpoints
  IndexPairVector endpoints;               // segment -> from_vertex, to_vertex
  std::vector<IndexPairVector> adjacency;  // segment -> (quad,side)*
  std::vector<JetWrapper::JetData> jet;    // jet-fitted principal curvatures
};
