#include "EmptyNeutronicsMaterial.h"

registerMooseObject("GnatApp", EmptyNeutronicsMaterial);

InputParameters
EmptyNeutronicsMaterial::validParams()
{
  auto params = Material::validParams();
  params.addClassDescription("Empty material to declare the required material "
                             "properties for Gnat's neutron transport solver. "
                             "All other neutronics materials must inherit from "
                             "this base class.");
  params.addParam<unsigned int>(
      "num_groups", 0u, "The number of neutron energy groups the energy spectrum is divided into.");
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");
  params.addParam<MooseEnum>("particle_type",
                             MooseEnum("neutron photon", "neutron"),
                             "The type of particle to be consuming material property data.");
  params.addParam<bool>(
      "is_saaf", true, "Whether or not the transport system is initializing a SAAF scheme.");
  params.addParam<bool>(
      "is_diffusion",
      false,
      "Whether or not the transport system is initializing a diffusion approximation scheme.");
  params.addParam<bool>(
      "has_fission", false, "Whether or not the transport system accounts for fission reactions.");
  params.addParam<Real>("saaf_eta",
                        0.5,
                        "Stabilization parameter eta for the SAAF-CFEM scheme. "
                        "eta = 0 and c = 1 is equivalent to the SAAF approach "
                        "with no void treatment. Use with caution.");
  params.addParam<Real>("saaf_c",
                        1.0,
                        "Stabilization parameter c for the SAAF-CFEM scheme. "
                        "eta = 0 and c = 1 is equivalent to the SAAF approach "
                        "with no void treatment. Use with caution.");

  return params;
}

EmptyNeutronicsMaterial::EmptyNeutronicsMaterial(const InputParameters & parameters)
  : Material(parameters),
    _scaling(1.0),
    _num_groups(getParam<unsigned int>("num_groups")),
    _particle(getParam<MooseEnum>("particle_type").getEnum<Particletype>()),
    _is_saaf(getParam<bool>("is_saaf")),
    _is_diffusion(getParam<bool>("is_diffusion")),
    _has_fission(getParam<bool>("has_fission")),
    _mat_inv_v_g(declareADProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                      "inv_v_g")),
    _mat_sigma_t_g(declareADProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g")),
    _mat_surface_source(declareADProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "surface_source")),
    _mat_sigma_s_g_prime_g_l(declareADProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "scattering_matrix")),
    _mat_anisotropy(declareProperty<unsigned int>(getParam<std::string>("transport_system") +
                                                  "medium_anisotropy")),
    _mat_source_moments(declareADProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "source_moments")),
    _mat_src_anisotropy(declareProperty<unsigned int>(getParam<std::string>("transport_system") +
                                                      "medium_source_anisotropy")),
    _mat_sigma_r_g(_is_diffusion ? &declareADProperty<std::vector<Real>>(
                                       getParam<std::string>("transport_system") + "removal_xs_g")
                                 : nullptr),
    _mat_diffusion_g(_is_diffusion ? &declareADProperty<std::vector<Real>>(
                                         getParam<std::string>("transport_system") + "diffusion_g")
                                   : nullptr),
    _mat_nu_sigma_f_g(_has_fission
                          ? &declareADProperty<std::vector<Real>>(
                                getParam<std::string>("transport_system") + "production_xs_g")
                          : nullptr),
    _mat_chi_f_g(_has_fission ? &declareADProperty<std::vector<Real>>(
                                    getParam<std::string>("transport_system") + "fission_spectra_g")
                              : nullptr),
    _mat_saaf_tau(_is_saaf ? &declareADProperty<std::vector<Real>>(
                                 getParam<std::string>("transport_system") + "saaf_tau")
                           : nullptr),
    _saaf_eta(getParam<Real>("saaf_eta")),
    _saaf_c(getParam<Real>("saaf_c"))
{
  if (_num_groups == 0u)
    mooseError("The provided number of energy groups is zero.");
}

void
EmptyNeutronicsMaterial::computeQpProperties()
{
  _mat_anisotropy[_qp] = 0u;
  _mat_src_anisotropy[_qp] = 0u;
}
