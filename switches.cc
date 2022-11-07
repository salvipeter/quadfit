#include "switches.hh"

template <>
int extractValue(std::string str) {
  return std::atoi(str.c_str());
}

template <>
size_t extractValue(std::string str) {
  return std::atoi(str.c_str());
}

template <>
double extractValue(std::string str) {
  return std::strtod(str.c_str(), nullptr);
}

template <>
std::string extractValue(std::string str) {
  return str;
}
