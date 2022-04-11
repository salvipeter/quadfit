#include <geometry.hh>

class QuadFit {
public:
  void readPWGB(std::string filename);
  std::vector<Geometry::BSSurface> fit();

private:
  using Ribbon = std::array<Geometry::BSCurve, 2>;
  struct Quad {
    struct Boundary {
      size_t segment;
      bool reversed;
      bool on_ribbon;
      size_t ribbon;       // when on_ribbon
      double s_min, s_max; //   is true
    };
    std::array<Boundary, 4> boundaries; // u=0, v=0, u=1, v=1
    size_t resolution;
    std::vector<Geometry::Point3D> samples;
  };

  size_t n;
  std::vector<Ribbon> ribbons;
  std::vector<Geometry::BSCurve> segments;
  std::vector<Quad> quads;
};
