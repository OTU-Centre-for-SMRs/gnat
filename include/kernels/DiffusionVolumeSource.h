#pragma once

#include "Kernel.h"

// A class which computes the external source contribution to the particle diffusion equation using
// moments provided by the user.
class DiffusionVolumeSource : public Kernel
{
public:
  static InputParameters validParams();

  DiffusionVolumeSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  // The current group (g) and the number of spectral energy groups (G).
  const unsigned int _group_index;
  const unsigned int _num_groups;

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
   * memory. The user is expected to format them according to this arrangement. We use this format
   * for the diffusion solver to ensure compatibility with input decks that provide source moments
   * for the transport solvers. This kernel only uses the 0th degree moment of the external source.
   */
  const std::vector<Real> _source_moments;
}; // class DiffusionVolumeSource
