#pragma once

#include "SAAFDiracKernelBase.h"

class SAAFIsoPointSource : public SAAFDiracKernelBase
{
public:
  static InputParameters validParams();

  SAAFIsoPointSource(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const std::vector<Real> _source_intensities; // S_{g, 0, 0} for all points.
  const std::vector<Point> _source_locations;
  std::map<Point, unsigned int> _point_intensity_mapping;
}; // class SAAFIsoPointSource
