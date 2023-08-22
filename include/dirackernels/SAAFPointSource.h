#pragma once

#include "SAAFBaseDiracKernel.h"

class SAAFPointSource : public SAAFBaseDiracKernel
{
public:
  static InputParameters validParams();

  SAAFPointSource(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;

  // Number of spectral energy groups (G).
  const unsigned int _num_groups;

  const std::vector<Real> & _source_moments;
  const Point _source_location;
  const unsigned int _anisotropy;

  // Storage for the pre-computed spherical harmonics coefficients (Y_{l,m,n}).
  // They are stored in the following order: l -> m.
  std::vector<Real> _y_l_m;
};
