#pragma once

#include <geometry.hh>

enum class DiscreteMask { COONS = 0, HARMONIC, BIHARMONIC };

void applyMask(Geometry::BSSurface &surface, DiscreteMask type);
