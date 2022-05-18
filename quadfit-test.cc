#include <fstream>

#include "io.hh"
#include "quadfit.hh"

using namespace Geometry;

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input.pwgb> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 50;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  QuadFit qf;
  auto description = qf.readPWGB(argv[1]);
  std::cout << "Input file \"" << argv[1] << "\" read:" << std::endl;
  std::cout << "\t" << description << std::endl;

  auto surfaces = qf.fit();
  std::cout << "Fit of " << surfaces.size() << " surfaces completed" << std::endl;

  std::cout << "Output resolution: " << resolution << std::endl;
  writeSTL(surfaces, "/tmp/surfaces.stl", resolution);
  writeControlNet(surfaces, "/tmp/controls.obj");
  writeBoundaries(surfaces, "/tmp/boundaries.obj", resolution);
  std::cout << "Surfaces written to /tmp/surfaces.stl" << std::endl;
  std::cout << "Control nets written to /tmp/controls.obj" << std::endl;
  std::cout << "Boundaries written to /tmp/boundaries.obj" << std::endl;
}
