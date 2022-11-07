#pragma once

#include <stdexcept>
#include <string>
#include <vector>

template <typename T>
T extractValue(std::string str) {
  throw std::runtime_error("extractValue not defined for this type");
}

template <> int extractValue(std::string str);
template <> size_t extractValue(std::string str);
template <> double extractValue(std::string str);
template <> std::string extractValue(std::string str);

template <typename T>
bool parseSwitch(const std::vector<std::string> &switches, std::string s,
                 T *value = nullptr, T default_value = T()) {
  size_t length = s.size();
  for (const auto &sw : switches)
    if (sw.substr(2, length) == s) {
      if (value) {
        if (sw.size() > length && sw[length+2] == '=')
          *value = extractValue<T>(sw.substr(length + 3));
        else
          *value = default_value;
      }
      return true;
    }
  if (value)
    *value = default_value;
  return false;
}
