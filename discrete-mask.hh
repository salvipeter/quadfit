#pragma once

#include <geometry.hh>

enum class DiscreteMask { C0_COONS = 0, C1_COONS, HARMONIC, BIHARMONIC };

void applyMask(Geometry::BSSurface &surface, DiscreteMask type);
