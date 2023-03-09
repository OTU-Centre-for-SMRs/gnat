#pragma once

#include "SAAFBaseKernel.h"

class SAAFScattering : public SAAFBaseKernel
{
public:
  static InputParameters validParams();

  SAAFScattering(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  Real computeFluxMoment(unsigned int g_prime, unsigned int sh_offset);

  const unsigned int _num_groups;     // G
  const unsigned int _max_anisotropy; // L
  unsigned int _num_dir_sh;           // Number of spherical harmonics evaluations per direction.

  /*
   * We assume that the vector of flux ordinates is stored in order of group
   * first, direction second. An example for 2 energy groups (G = 2) and a
   * quadrature set with 2 elements (N = 2) is given below:
   *
   * The flux ordinates are indexed as Psi_{g, n}:
   * _group_flux_ordinates[0] = Psi_{1, 1}
   * _group_flux_ordinates[0] = Psi_{1, 2}
   * _group_flux_ordinates[0] = Psi_{2, 1}
   * _group_flux_ordinates[0] = Psi_{2, 1}
   */
  std::vector<unsigned int> _group_flux_ordinate_ids;
  std::vector<const VariableValue *> _group_flux_ordinates;

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

  // Storage for the pre-computed spherical harmonics coefficients (Y_{l,m,n}).
  // They are stored in the following order: n -> l -> m.
  std::vector<Real> _y_l_m_n;
}; // class SAAFScattering
