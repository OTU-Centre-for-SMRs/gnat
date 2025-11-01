#pragma once

#include <vector>

#include "Moose.h"
#include "MooseTypes.h"

// A namespace to store helper functions involving the real spherical harmonics.
// At the present, this just includes a function to statically evaluate
// spherical harmonics (supporting degrees up to 3).
namespace RealSphericalHarmonics
{
  Real evaluate(unsigned int degree, int order, const Real & mu, const Real & eta, const Real & xi);
}; // namespace RealSphericalHarmonics
