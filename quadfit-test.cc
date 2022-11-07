#include <fstream>

#include "io.hh"
#include "quadfit.hh"
#include "switches.hh"

using namespace Geometry;

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input.pwgb> [switches]" << std::endl;
    return 1;
  }

  std::vector<std::string> switches;
  for (int i = 2; i < argc; ++i) {
    if (argv[i][0] != '-' || argv[i][1] != '-') {
      std::cerr << "Invalid switch: " << argv[i] << std::endl;
      return 2;
    }
    switches.push_back(argv[i]);
  }

  size_t resolution = 0;
  parseSwitch<size_t>(switches, "resolution", &resolution, 50);

  QuadFit qf;
  auto description = qf.readPWGB(argv[1]);
  std::cout << "Input file \"" << argv[1] << "\" read:" << std::endl;
  std::cout << "\t" << description << std::endl;

  qf.update(switches);
  std::cout << "Topology & surface curvatures updated" << std::endl;

  auto surfaces = qf.fit(switches);
  std::cout << "Fit of " << surfaces.size() << " surfaces completed" << std::endl;

  std::cout << "Output resolution: " << resolution << std::endl;
  writeQDS(surfaces, "/tmp/surfaces.qds");
  writeSTL(surfaces, "/tmp/surfaces.stl", resolution);
  writeControlNet(surfaces, "/tmp/controls.obj");
  writeBoundaries(surfaces, "/tmp/boundaries.obj", resolution);
  std::cout << "Surface data written to /tmp/surfaces.qds" << std::endl;
  std::cout << "Surface mesh written to /tmp/surfaces.stl" << std::endl;
  std::cout << "Control nets written to /tmp/controls.obj" << std::endl;
  std::cout << "Boundaries written to /tmp/boundaries.obj" << std::endl;
}
