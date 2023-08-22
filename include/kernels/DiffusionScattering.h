#pragma once

#include "Kernel.h"

class DiffusionScattering : public Kernel
{
public:
  static InputParameters validParams();

  DiffusionScattering(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // The current group (g) and the number of spectral energy groups (G).
  const unsigned int _group_index;
  const unsigned int _num_groups;

  /*
   * Vector of the scalar fluxes (0th degree and moment of the angular flux). We assume that the
   * vector of scalar fluxes is ordered by energy group. An example with 4 energy groups (G = 2) can
   * be found below:
   *
   * The scalar fluxes are indexed as Phi_{g}:
   * _group_scalar_fluxes[0] = Phi_{1}
   * _group_scalar_fluxes[1] = Phi_{2}
   * _group_scalar_fluxes[2] = Phi_{3}
   * _group_scalar_fluxes[3] = Phi_{4}
   */
  std::map<unsigned int, unsigned int> _jvar_map;
  std::vector<const VariableValue *> _group_scalar_fluxes;

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
   * according to this arrangement. We use this format for the diffusion solver to ensure
   * compatibility with materials that provide scattering moments for the transport solvers. This
   * kernel only uses the 0th degree moments of the scattering cross-section.
   */
  const ADMaterialProperty<std::vector<Real>> & _sigma_s_g_prime_g_l;
  // Degree of anisotropy (Legendre polynomial order L) for the medium.
  // Only needed here to index the scattering matix correctly.
  const MaterialProperty<unsigned int> & _anisotropy;
}; // class DiffusionScattering
