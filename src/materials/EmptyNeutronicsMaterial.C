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
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups > 0",
                                                    "The number of neutron "
                                                    "energy groups the energy "
                                                    "spectrum is divided into.");
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
    _num_groups(getParam<unsigned int>("num_groups")),
    _saaf_eta(getParam<Real>("saaf_eta")),
    _saaf_c(getParam<Real>("saaf_c")),
    _mat_inv_v_g(declareADProperty<std::vector<Real>>("inv_v_g")),
    _mat_sigma_t_g(declareADProperty<std::vector<Real>>("total_xs_g")),
    _mat_sigma_r_g(declareADProperty<std::vector<Real>>("removal_xs_g")),
    _mat_diffusion_g(declareADProperty<std::vector<Real>>("diffusion_g")),
    _mat_surface_source(declareADProperty<std::vector<Real>>("surface_source")),
    _mat_sigma_s_g_prime_g_l(declareADProperty<std::vector<Real>>("scattering_matrix")),
    _mat_anisotropy(declareProperty<unsigned int>("medium_anisotropy")),
    _mat_source_moments(declareADProperty<std::vector<Real>>("source_moments")),
    _mat_src_anisotropy(declareProperty<unsigned int>("medium_source_anisotropy")),
    _mat_saaf_eta(declareADProperty<Real>("saaf_eta")),
    _mat_saaf_c(declareADProperty<Real>("saaf_c"))
{
}

void
EmptyNeutronicsMaterial::computeQpProperties()
{
  _mat_anisotropy[_qp] = 0u;
  _mat_src_anisotropy[_qp] = 0u;

  // SAAF stabilization properties.
  _mat_saaf_eta[_qp] = _saaf_eta;
  _mat_saaf_c[_qp] = _saaf_c;
}
