#pragma once

#include "GnatBase.h"
#include "EmptyTransportMaterial.h"

class FileTransportMaterial : public EmptyTransportMaterial
{
public:
  static InputParameters validParams();

  FileTransportMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // A helper function to turn a string vector of cross-sections deliminated with a ' ' into a
  // vector of numerics.
  static void parseToVector(const std::string & string_rep, std::vector<Real> & real_rep);
  void parseXMLMacroXS();

  // Cross-section information.
  enum class XSUnits
  {
    InvCm = 0u
  } _xs_units;
  enum class EnergyUnits
  {
    eV = 0u,
    keV = 1u,
    MeV = 2u
  } _energy_units;
  std::vector<Real> _group_bounds;

  // The sum of all isotopic material properties are the actual properties provided to the transport
  // solver.
  std::vector<Real> _inv_v_g;
  std::vector<Real> _sigma_t_g;
  std::vector<Real> _sigma_s_g_matrix;
  std::vector<Real> _sigma_s_g_g;
  std::vector<Real> _diffusion_g;
  std::vector<Real> _nu_sigma_f_g;
  std::vector<Real> _chi_f_g;
  std::vector<Real> _heating_g;
  /*
   * We assume that the vector of scattering cross-sections is stored in the
   * order of initial group first, final group second, and Legendre polynomial
   * third. An example for two energy groups (G = 2) and an anisotropy of 2
   * (Legendre polynomial order L = 2) is given below.
   *
   * The scattering cross-sections are indexed as SigmaS_{g', g, l}:
   * _sigma_s_g_prime_g_l[0] = SigmaS_{1, 1, 0}
   * _sigma_s_g_prime_g_l[1] = SigmaS_{1, 1, 1}
   * _sigma_s_g_prime_g_l[2] = SigmaS_{1, 1, 2}
   * _sigma_s_g_prime_g_l[3] = SigmaS_{1, 2, 0}
   * _sigma_s_g_prime_g_l[4] = SigmaS_{1, 2, 1}
   * _sigma_s_g_prime_g_l[5] = SigmaS_{1, 2, 2}
   * _sigma_s_g_prime_g_l[6] = SigmaS_{2, 1, 0}
   * _sigma_s_g_prime_g_l[7] = SigmaS_{2, 1, 1}
   * _sigma_s_g_prime_g_l[8] = SigmaS_{2, 1, 2}
   * _sigma_s_g_prime_g_l[9] = SigmaS_{2, 2, 0}
   * _sigma_s_g_prime_g_l[10] = SigmaS_{2, 2, 1}
   * _sigma_s_g_prime_g_l[11] = SigmaS_{2, 2, 2}
   *
   * This is a flattening of the scattering matrix to preserve coherency in memory.
   * The material providing the cross-sections moments is expected to format them
   * according to this arrangement.
   */
  std::vector<Real> _sigma_s_g_prime_g_l;

  unsigned int _anisotropy;
  unsigned int _max_moments;

  const std::vector<Real> & _source_moments;
  const unsigned int _source_anisotropy;
  unsigned int _max_source_moments;
  bool _has_volumetric_source;

  const std::string _file_name;
  const std::string & _source_material_id;
}; // class FileTransportMaterial
