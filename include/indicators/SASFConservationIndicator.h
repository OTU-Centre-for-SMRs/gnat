#pragma once

#include "ElementIntegralIndicator.h"

class SASFConservationIndicator : public ElementIntegralIndicator
{
public:
  static InputParameters validParams();

  SASFConservationIndicator(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  // The source location.
  const Point _source_location;

  // g
  const unsigned int _group_index;
  // Total cross-section.
  const ADMaterialProperty<std::vector<Real>> & _sigma_t_g;
};
