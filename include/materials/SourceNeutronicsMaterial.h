#pragma once

#include "ConstantNeutronicsMaterial.h"

class SourceNeutronicsMaterial : public ConstantNeutronicsMaterial
{
public:
  static InputParameters validParams();

  SourceNeutronicsMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  std::vector<Real> _source_moments;

  const unsigned int _source_anisotropy;
  unsigned int _max_source_moments;
};
