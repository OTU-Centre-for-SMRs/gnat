#pragma once

#include "ADSNBaseKernel.h"

class ADDFEMMaterialSource : public ADSNBaseKernel
{
public:
  static InputParameters validParams();

  ADDFEMMaterialSource(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n
  const unsigned int _group_index; // g
  const unsigned int _num_groups; // G

  /*
   * We assume that the vector of source moments is stored in order of group first,
   * then moment indices. As an example for 2 energy groups (G = 2) and a 2nd
   * order real spherical harmonics source expansion with L = 2 is given below.
   *
   * The source moments are indexed as S_{g, l, m}:
   * _source_moments[_qp][0] = S_{1, 0, 0}
   * _source_moments[_qp][1] = S_{1, 1, -1}
   * _source_moments[_qp][2] = S_{1, 1, 0}
   * _source_moments[_qp][3] = S_{1, 1, 1}
   * _source_moments[_qp][4] = S_{1, 1, 1}
   * _source_moments[_qp][5] = S_{1, 2, -2}
   * _source_moments[_qp][6] = S_{1, 2, -1}
   * _source_moments[_qp][7] = S_{1, 2, 0}
   * _source_moments[_qp][8] = S_{1, 2, 1}
   * _source_moments[_qp][9] = S_{1, 2, 2}
   * _source_moments[_qp][10] = S_{2, 0, 0}
   * _source_moments[_qp][11] = S_{2, 1, -1}
   * _source_moments[_qp][12] = S_{2, 1, 0}
   * _source_moments[_qp][13] = S_{2, 1, 1}
   * _source_moments[_qp][14] = S_{2, 1, 1}
   * _source_moments[_qp][15] = S_{2, 2, -2}
   * _source_moments[_qp][16] = S_{2, 2, -1}
   * _source_moments[_qp][17] = S_{2, 2, 0}
   * _source_moments[_qp][18] = S_{2, 2, 1}
   * _source_moments[_qp][19] = S_{2, 2, 2}
   *
   * This is a flattening of the source moments matrix to preserve coherency in
   * memory. The material providing the source moments is expected to format them
   * according to this arrangement.
  */
  const ADMaterialProperty<std::vector<Real>> & _source_moments;
  // Degree of anisotropy (Legendre polynomial order L) for the material source.
  const MaterialProperty<unsigned int> & _anisotropy;
}; // class ADDFEMMaterialSource
