#pragma once

#include "InternalSideIndicator.h"

class SASFStreamingJumpIndicator : public InternalSideIndicator
{
public:
  static InputParameters validParams();

  SASFStreamingJumpIndicator(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  // The source location.
  const Point _source_location;
};
