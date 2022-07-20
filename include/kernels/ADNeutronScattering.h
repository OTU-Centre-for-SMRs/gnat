#pragma once

#include "ADNeutronBaseKernel.h"

// This kernel evaluates the full scattering contribution. This way the complete
// matrix is assembled and source iteration is avoided. This kernel should only
// be used for debugging purposes as a fully assembled matrix solve for the transport
// equation is quite slow.
// TODO: Finish this kernel.
class ADNeutronScattering : public ADNeutronBaseKernel
{
public:
  static InputParameters validParams();

  ADNeutronScattering(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  ADReal computeFluxMoment(unsigned int g_prime, unsigned int l, int m);

  const unsigned int _ordinate_index; // n
  const unsigned int _group_index; // g
  const unsigned int _num_groups; // G
  const unsigned int _max_anisotropy; // L

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
};
