#pragma once

#include "InternalSideIndicator.h"

class SASFTransverseJumpIndicator : public InternalSideIndicator
{
public:
  static InputParameters validParams();

  SASFTransverseJumpIndicator(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  // The source location.
  const Point _source_location;
};
