#pragma once

#include "SAAFBaseKernel.h"

// A class which implements the anisotropic particle source term for a field-based source in the
// SAAF particle transport equation. Use cases include the radiation fields off of a radionuclide
// plume.
class SAAFFieldSource : public SAAFBaseKernel
{
public:
  static InputParameters validParams();

  SAAFFieldSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  // Number of spectral energy groups (G).
  const unsigned int _num_groups;

  /*
   * We assume that the vector of source moments is stored in order of group first,
   * then moment indices. As an example for 2 energy groups (G = 2) and a 2nd
   * order real spherical harmonics source expansion with L = 2 is given below.
   *
   * The source moments are indexed as S_{g, l, m}:
   * _source_moments[0] = S_{1, 0, 0}
   * _source_moments[1] = S_{1, 1, -1}
   * _source_moments[2] = S_{1, 1, 0}
   * _source_moments[3] = S_{1, 1, 1}
   * _source_moments[4] = S_{1, 1, 1}
   * _source_moments[5] = S_{1, 2, -2}
   * _source_moments[6] = S_{1, 2, -1}
   * _source_moments[7] = S_{1, 2, 0}
   * _source_moments[8] = S_{1, 2, 1}
   * _source_moments[9] = S_{1, 2, 2}
   * _source_moments[10] = S_{2, 0, 0}
   * _source_moments[11] = S_{2, 1, -1}
   * _source_moments[12] = S_{2, 1, 0}
   * _source_moments[13] = S_{2, 1, 1}
   * _source_moments[14] = S_{2, 1, 1}
   * _source_moments[15] = S_{2, 2, -2}
   * _source_moments[16] = S_{2, 2, -1}
   * _source_moments[17] = S_{2, 2, 0}
   * _source_moments[18] = S_{2, 2, 1}
   * _source_moments[19] = S_{2, 2, 2}
   *
   * This is a flattening of the source moments matrix to preserve coherency in
   * memory. The user is expected to format them according to this arrangement.
   */
  std::vector<const VariableValue *> _source_moments;
  // Degree of anisotropy (Legendre polynomial order L) for the material source.
  const unsigned int _anisotropy;

  // Storage for the pre-computed spherical harmonics coefficients (Y_{l,m,n}).
  // They are stored in the following order: l -> m.
  std::vector<Real> _y_l_m;

  const Real _scale_factor;
}; // class SAAFFieldSource
