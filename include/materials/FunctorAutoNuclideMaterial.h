#pragma once

#include "FunctorMaterial.h"

// A class to compute the spatially-varying diffusion coefficients required by the radionuclide
// advection-diffusion-reaction equation using functors. Currently implements purely laminar and
// mixing length models.
class FunctorAutoNuclideMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  FunctorAutoNuclideMaterial(const InputParameters & parameters);

protected:
  // The Brownian diffusion coefficient is equal to:
  // D_B = \frac{k_B T}{6 \pi \rho \nu r}, \mu = \rho \nu
  // The eddy diffusivity is equal to:
  // D_t = \nu_{t} / Sc_{t}

  static constexpr Real _k_boltzmann_m = 1.380649e-23;
  static constexpr Real _k_boltzmann_cm = 1.380649e-19;

  const unsigned int _mesh_dims;

  // The nuclide scheme.
  enum class NuclideScheme
  {
    SUPGFE = 0u,
    FV = 1u
  } _scheme;

  // TODO: Have the diffusion coefficient rely on fluid properties.
  enum class CoefficientType
  {
    Laminar = 0u,
    MixingLength = 1u
  } _diff_type;

  const Real _radii;

  const Moose::Functor<ADReal> & _temperature;
  const Moose::Functor<ADReal> & _visc_dynamic;

  // Turbulent Schmidt number for D_{t}.
  const Real _schmidt_number;

  // Required properties for the mixing-length turbulence model.
  // Turbulent mixing length.
  const Moose::Functor<ADReal> * const _mixing_len;
  // X-velocity.
  const Moose::Functor<ADReal> * const _u;
  // Y-velocity.
  const Moose::Functor<ADReal> * const _v;
  // Z-velocity.
  const Moose::Functor<ADReal> * const _w;

  // The diffusion coefficient.
  const Moose::Functor<ADReal> * _nuc_diff;
};
