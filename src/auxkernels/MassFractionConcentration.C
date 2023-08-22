#include "MassFractionConcentration.h"

#include "Nuclide.h"

registerMooseObject("GnatApp", MassFractionConcentration);

InputParameters
MassFractionConcentration::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("An auxkernel for the purposes of computing a radionuclide density "
                             "from a mass fraction for post-processing. This "
                             "is either a number density or a mass density.");
  params.addRequiredParam<bool>("is_fe", "Whether the mass fraction is finite element or not.");
  params.addCoupledVar("nuclide_mass_fraction_var", "The radionuclide mass fraction.");
  params.addParam<MooseFunctorName>("nuclide_mass_fraction_fun", "The radionuclide mass fraction.");
  params.addRequiredParam<MooseFunctorName>("density", "The density of the bulk fluid.");
  params.addParam<bool>("compute_number_density",
                        true,
                        "Whether a number density should be computed or not. Setting this flag to "
                        "false results in a mass density instead.");
  params.addParam<std::string>("nuclide_name",
                               "",
                               "The name of the nuclide this mass fraction represents. If not set, "
                               "the name of the 'nuclide_mass_fraction' functor will be used.");

  return params;
}

MassFractionConcentration::MassFractionConcentration(const InputParameters & parameters)
  : AuxKernel(parameters),
    _is_fe(getParam<bool>("is_fe")),
    _output_number_density(getParam<bool>("compute_number_density")),
    _density(getFunctor<ADReal>("density")),
    _mass_fraction_var(nullptr),
    _mass_fraction_fun(nullptr)
{
  if (_is_fe)
  {
    _mass_fraction_var = &adCoupledValue("nuclide_mass_fraction_var");

    if (getParam<std::string>("nuclide_name") == "")
      _weight = NuclearData::Nuclide::getAtomicMass(coupledName("nuclide_mass_fraction_var"));
    else
      _weight = NuclearData::Nuclide::getAtomicMass(getParam<std::string>("nuclide_name"));
  }
  else
  {
    _mass_fraction_fun = &getFunctor<ADReal>("nuclide_mass_fraction_fun");

    if (getParam<std::string>("nuclide_name") == "")
      _weight = NuclearData::Nuclide::getAtomicMass(
          getParam<MooseFunctorName>("nuclide_mass_fraction_fun"));
    else
      _weight = NuclearData::Nuclide::getAtomicMass(getParam<std::string>("nuclide_name"));
  }
}

Real
MassFractionConcentration::computeValue()
{
  const auto elem_arg = makeElemArg(_current_elem);

  Real val = 0.0;
  if (_is_fe)
    val = MetaPhysicL::raw_value((*_mass_fraction_var)[_qp]) *
          MetaPhysicL::raw_value(_density(elem_arg, 0u));
  else
    val = MetaPhysicL::raw_value((*_mass_fraction_fun)(elem_arg, 0u)) *
          MetaPhysicL::raw_value(_density(elem_arg, 0u));

  if (_output_number_density)
    val *= _n_avogadro / _weight;

  return val;
}
