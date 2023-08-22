#pragma once

#include "DiracKernel.h"

class DiffusionIsoPointSource : public DiracKernel
{
public:
  static InputParameters validParams();

  DiffusionIsoPointSource(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;

  // Number of spectral energy groups (G).
  const unsigned int _num_groups;
  const unsigned int _group_index; // g

  const std::vector<Real> & _source_moments;
  const Point _source_location;
  const unsigned int _anisotropy;
}; // class DiffusionIsoPointSource
