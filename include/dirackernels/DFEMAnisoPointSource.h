#pragma once

#include "SNBaseDiracKernel.h"

/*
 * This kernel assumes that the external point sources take the form:
 * S_{g}(\hat{\Omega}) = f_{g}(\hat{\Omega})S_{g}. It is assumed that the
 * integral of f_{g}(\hat{\Omega}) over the unit sphere equals 1.
 */
// TODO: Multiple anisotropic point sources?
class DFEMAnisoPointSource : public SNBaseDiracKernel
{
public:
  static InputParameters validParams();

  DFEMAnisoPointSource(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const unsigned int _ordinate_index;

  const Real _source_intensity;           // S_{g}
  const Function & _angular_distribution; // f_{g}(\hat{\Omega})
  const Point _source_location;
}; // class DFEMAnisoPointSource
