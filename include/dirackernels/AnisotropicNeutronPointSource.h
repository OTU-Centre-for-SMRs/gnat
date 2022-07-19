#pragma once

#include "DiracKernel.h"
#include "GnatBase.h"

/*
 * This kernel assumes that the external point sources take the form:
 * S_{g}(\hat{\Omega}) = f_{g}(\hat{\Omega})S_{g}. It is assumed that the
 * integral of f_{g}(\hat{\Omega}) over the unit sphere equals 1.
*/
// TODO: Multiple anisotropic point sources?
class AnisotropicNeutronPointSource : public DiracKernel
{
public:
  static InputParameters validParams();

  AnisotropicNeutronPointSource(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;

  const ProblemType _type;
  Real _symmetry_factor;

  const unsigned int _ordinate_index;

  const Real _source_intensity; // S_{g}
  const Function & _angular_distribution; // f_{g}(\hat{\Omega})
  const Point _source_location;

  const MaterialProperty<std::vector<RealVectorValue>> & _directions;
}; // class AnisotropicNeutronPointSource
