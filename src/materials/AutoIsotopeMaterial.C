#include "AutoIsotopeMaterial.h"

InputParameters
AutoIsotopeMaterial::validParams()
{
  auto params = Material::validParams();
  params.addClassDescription(
      "This material provides the diffusion coefficient for isotope mass transport problems. It "
      "should not be exposed to the user, but should be automatically generated using a mass "
      "transport action.");

  params.addRequiredParam<NonlinearVariableName>("isotope_name",
                                                 "The name of the isotope variable.");
  params.addRequiredParam<MooseEnum>("diffusion_coefficient_type",
                                     MooseEnum("constant"),
                                     "The type of diffusion coefficient to use.");

  params.addParam<Real>(
      "diffusion_coefficient_base", 0.0, "A simple constant diffusion coefficient.");

  return params;
}

AutoIsotopeMaterial::AutoIsotopeMaterial(const InputParameters & parameters)
  : Material(parameters),
    _diff_type(getParam<MooseEnum>("diffusion_coefficient_type").getEnum<CoefficientType>()),
    _diffusion_coefficient(getParam<Real>("diffusion_coefficient_base")),
    _mat_diff(declareADProperty<Real>(
        "isotope_diff_" + Moose::stringify(getParam<NonlinearVariableName>("isotope_name")))),
    _grad_mat_diff(declareADProperty<RealVectorValue>(
        "grad_isotope_diff_" + Moose::stringify(getParam<NonlinearVariableName>("isotope_name"))))
{
}

void
AutoIsotopeMaterial::computeQpProperties()
{
  switch (_diff_type)
  {
    case CoefficientType::Constant:
      _mat_diff[_qp] = _diffusion_coefficient;
      _grad_mat_diff[_qp] = RealVectorValue(0.0);
      break;

    default:
      break;
  }
}
