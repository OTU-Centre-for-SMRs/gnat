#pragma once

#include "ADSNBaseKernel.h"

class ADDFEMInternalScattering : public ADSNBaseKernel
{
public:
  static InputParameters validParams();

  ADDFEMInternalScattering(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n
  const unsigned int _group_index;    // g
  const unsigned int _num_groups;     // G

  /*
   * We assume that the vector of all group flux moments is stored in order of
   * moment indices. An example for a spherical harmonics expansion of L = 2 is
   * provided below.
   *
   * The within-group flux moments are indexed as Phi_{l, m}:
   * _within_group_flux_moments[_qp][0] = Phi_{0, 0}
   * _within_group_flux_moments[_qp][1] = Phi_{1, -1}
   * _within_group_flux_moments[_qp][2] = Phi_{1, 0}
   * _within_group_flux_moments[_qp][3] = Phi_{1, 1}
   * _within_group_flux_moments[_qp][4] = Phi_{1, 1}
   * _within_group_flux_moments[_qp][5] = Phi_{2, -2}
   * _within_group_flux_moments[_qp][6] = Phi_{2, -1}
   * _within_group_flux_moments[_qp][7] = Phi_{2, 0}
   * _within_group_flux_moments[_qp][8] = Phi_{2, 1}
   * _within_group_flux_moments[_qp][9] = Phi_{2, 2}
   *
   * This is a flattening of the flux moments matrix to preserve
   * coherency in memory. The action providing the moments is expected
   * to format them according to this arrangement.
   */
  std::vector<const ADVariableValue *> _within_group_flux_moments;
  unsigned int _provided_moment_degree;

  /*
   * We assume that the vector of scattering cross-sections is stored in the
   * order of initial group first, final group second, and Legendre polynomial
   * third. An example for two energy groups (G = 2) and an anisotropy of 2
   * (Legendre polynomial order L = 2) is given below.
   *
   * The scattering cross-sections are indexed as SigmaS_{g', g, l}:
   * _sigma_s_g_prime_g_l[_qp][0] = SigmaS_{1, 1, 0}
   * _sigma_s_g_prime_g_l[_qp][1] = SigmaS_{1, 1, 1}
   * _sigma_s_g_prime_g_l[_qp][2] = SigmaS_{1, 1, 2}
   * _sigma_s_g_prime_g_l[_qp][3] = SigmaS_{1, 2, 0}
   * _sigma_s_g_prime_g_l[_qp][4] = SigmaS_{1, 2, 1}
   * _sigma_s_g_prime_g_l[_qp][5] = SigmaS_{1, 2, 2}
   * _sigma_s_g_prime_g_l[_qp][6] = SigmaS_{2, 1, 0}
   * _sigma_s_g_prime_g_l[_qp][7] = SigmaS_{2, 1, 1}
   * _sigma_s_g_prime_g_l[_qp][8] = SigmaS_{2, 1, 2}
   * _sigma_s_g_prime_g_l[_qp][9] = SigmaS_{2, 2, 0}
   * _sigma_s_g_prime_g_l[_qp][10] = SigmaS_{2, 2, 1}
   * _sigma_s_g_prime_g_l[_qp][11] = SigmaS_{2, 2, 2}
   *
   * This is a flattening of the scattering matrix to preserve coherency in memory.
   * The material providing the cross-sections moments is expected to format them
   * according to this arrangement.
   */
  const ADMaterialProperty<std::vector<Real>> & _sigma_s_g_prime_g_l;
  // Degree of anisotropy (Legendre polynomial order L) for the medium.
  const MaterialProperty<unsigned int> & _anisotropy;
}; // class ADDFEMInternalScattering
