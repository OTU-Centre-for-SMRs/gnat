#pragma once

#include "EmptyNeutronicsMaterial.h"

class FileNeutronicsMaterial : public EmptyNeutronicsMaterial
{
public:
  static InputParameters validParams();

  FileNeutronicsMaterial(const InputParameters & parameters);

protected:
  enum class PropertyType
  {
    InvVelocity = 0u,
    SigmaR = 1u,
    SigmaA = 2u,
    SigmaS = 3u,
    SigmaSMatrix = 4u
  };

  virtual void computeQpProperties() override;
  void parseProperty(const PropertyType & type, const std::string & property_file);
  void parseGnatProperty(const PropertyType & type, const std::string & property_file);
  void parseOpenMCProperty(const PropertyType & type, const std::string & property_file);

  enum class CrossSectionSource
  {
    Detect = 0u,
    Gnat = 1u,
    OpenMC = 2u
  } _xs_source;

  struct IsotopeProperties
  {
    std::vector<Real> _inv_v_g;
    std::vector<Real> _sigma_r_g;
    std::vector<Real> _sigma_a_g;
    std::vector<Real> _sigma_s_g;
    std::vector<Real> _sigma_s_g_prime_g_l;
  };

  // Nuclear properties for individual isotopes.
  std::unordered_map<std::string, IsotopeProperties> _material_properties;

  // The sum of all isotopic material properties are the actual properties provided to the transport
  // solver.
  std::vector<Real> _inv_v_g;
  std::vector<Real> _sigma_r_g;
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

  const std::string & _file_name;
  const std::string & _source_material_id;
}; // class FileNeutronicsMaterial
