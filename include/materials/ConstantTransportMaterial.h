#pragma once

#include "AbsorbingTransportMaterial.h"

class ConstantTransportMaterial : public AbsorbingTransportMaterial
{
public:
  static InputParameters validParams();

  ConstantTransportMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

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
   * The user is expected to format them according to this arrangement in the
   * input parameters.
   */
  // In-scattering cross-section moments.
  std::vector<Real> _sigma_s_g_prime_g_l;
  // Out-scattering scalar cross-section and within group scattering cross-section.
  std::vector<Real> _sigma_s_g;
  std::vector<Real> _sigma_s_g_g;

  std::vector<Real> _nu_sigma_f_g;
  std::vector<Real> _chi_f_g;

  const unsigned int _anisotropy;
  const unsigned int _max_moments;
}; // class ConstantTransportMaterial
