#include "NeutronTransportAction.h"

#include "Factory.h"
#include "Parser.h"
#include "NonlinearSystemBase.h"
#include "FEProblemBase.h"
#include "Conversion.h"
#include "MooseTypes.h"
#include "FEProblem.h"

#include "AddVariableAction.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/fe_type.h"

registerMooseAction("GnatApp", NeutronTransportAction, "add_variable");
registerMooseAction("GnatApp", NeutronTransportAction, "add_kernel");
registerMooseAction("GnatApp", NeutronTransportAction, "add_dg_kernel");
registerMooseAction("GnatApp", NeutronTransportAction, "add_dirac_kernel");
registerMooseAction("GnatApp", NeutronTransportAction, "add_bc");
registerMooseAction("GnatApp", NeutronTransportAction, "add_ic");
registerMooseAction("GnatApp", NeutronTransportAction, "add_material");
registerMooseAction("GnatApp", NeutronTransportAction, "add_aux_variable");
registerMooseAction("GnatApp", NeutronTransportAction, "add_aux_kernel");

InputParameters
NeutronTransportAction::validParams()
{
  auto params = Action::validParams();
  params.addClassDescription("This action adds all of the required variables, "
                             "kernels, boundary conditions, initial conditions "
                             "and auxiliary systems required to solve "
                             "source-driven multi-group neutron transport "
                             "problems with Gnat.");

  //----------------------------------------------------------------------------
  // Parameters for variables.
  // Get MooseEnums for the possible order/family options for this variable.
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  // Set the params required to add new variables. These will be used for both
  // nonlinear variables and auxvariables.
  params.addParam<MooseEnum>("family", families,
                             "Specifies the family of FE shape functions to "
                             "use for this variable.");
  params.addParam<MooseEnum>("order", orders,
                             "Specifies the order of the FE shape function to use "
                             "for this variable (additional orders not listed are "
                             "allowed).");
  params.addParam<Real>("scaling", 1.0,
                        "Specifies a scaling factor to apply to "
                        "this variable.");
  params.addParam<std::vector<SubdomainName>>("block",
                                              "The list of blocks (ids or "
                                              "names) that this variable will "
                                              "be applied.");

  //----------------------------------------------------------------------------
  // Basic parameters for the neutron transport simulation.
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups > 0",
                                                    "The number of spectral "
                                                    "energy groups in the "
                                                    "problem.");
  params.addRangeCheckedParam<unsigned int>("n_polar", 3, "n_polar > 0",
                                            "Number of Legendre polar "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 3.");
  params.addRangeCheckedParam<unsigned int>("n_azimuthal", 3, "n_azimuthal > 0",
                                            "Number of Chebyshev azimuthal "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 3.");
  params.addParam<unsigned int>("max_anisotropy", 0,
                                "The maximum degree of anisotropy to evaluate. "
                                "Defaults to 0 for isotropic scattering.");
  params.addParam<std::string>("angular_flux_names",
                               "angular_flux",
                               "Variable names for the angular flux. The output "
                               "format for the group angular fluxes will be of "
                               "the form {angular_flux_names}_g_n.");
  params.addParam<std::string>("flux_moment_names",
                               "flux_moment",
                               "Variable names for the moments of the angular "
                               "flux. The output format for the group flux "
                               "moments will be of the form "
                               "{flux_moment_names}_g_l_m.");

  //----------------------------------------------------------------------------
  // Boundary conditions.
  params.addParam<std::vector<BoundaryName>>("vacuum_boundaries",
                                             "The boundaries to apply vacuum "
                                             "boundary conditions.");
  params.addParam<std::vector<BoundaryName>>("source_boundaries",
                                             "The boundaries to apply incoming "
                                             "flux boundary conditions.");
  params.addParam<std::vector<BoundaryName>>("reflective_boundaries",
                                             "The boundaries to apply reflective "
                                             "boundary conditions.");

  //----------------------------------------------------------------------------
  // Parameters for debugging.
  params.addParam<bool>("debug_disable_scattering", false,
                        "Debug option to disable scattering evaluation.");

  return params;
}

NeutronTransportAction::NeutronTransportAction(const InputParameters & params)
  : Action(params)
  , _num_groups(getParam<unsigned int>("num_groups"))
  , _n_l(getParam<unsigned int>("n_polar"))
  , _n_c(getParam<unsigned int>("n_azimuthal"))
  , _max_eval_anisotropy(getParam<unsigned int>("max_anisotropy"))
  , _angular_flux_name(getParam<std::string>("angular_flux_names"))
  , _flux_moment_name(getParam<std::string>("flux_moment_names"))
  , _vacuum_side_sets(getParam<std::vector<BoundaryName>>("vacuum_boundaries"))
  , _source_side_sets(getParam<std::vector<BoundaryName>>("source_boundaries"))
  , _reflective_side_sets(getParam<std::vector<BoundaryName>>("reflective_boundaries"))
{
  // Grab all the subdomain IDs that the neutron transport action should be
  // applied to.
  std::vector<SubdomainName> block_names(getParam<std::vector<SubdomainName>>("block"));
  for (const auto & block_name : block_names)
    _subdomain_ids.insert(_problem->mesh().getSubdomainID(block_name));

  // Find out what sort of problem we're working with. Error if it's not a
  // cartesian coordinate system.
  for (const SubdomainID & id : _subdomain_ids)
  {
    if (_problem->getCoordSystem(id) != Moose::COORD_XYZ)
      mooseError("Neutron transport simulations currently do not support "
                 "non-cartesian coordinate systems.");
  }

  // Set the enum so quadrature sets can be determined appropriately.
  switch (_problem->mesh().getMesh().spatial_dimension())
  {
    case 1u:
      _p_type = ProblemType::Cartesian1D;
      break;

    case 2u:
      _p_type = ProblemType::Cartesian2D;
      break;

    case 3u:
      _p_type = ProblemType::Cartesian3D;
      break;

    default:
      mooseError("Unknown mesh dimensionality.");
      break;
  }
}

void
NeutronTransportAction::addVariable(const std::string & var_name)
{
  auto fe_type = AddVariableAction::feType(_pars);
  auto type = AddVariableAction::variableType(fe_type, false, false);
  auto var_params = _factory.getValidParams(type);

  var_params.applySpecificParameters(_pars, {"family", "order"});
  var_params.set<std::vector<Real>>("scaling") = std::vector<Real>(getParam<Real>("scaling"));

  if (!_subdomain_ids.empty())
  {
    for (const SubdomainID & id : _subdomain_ids)
      var_params.set<std::vector<SubdomainName>>("block").push_back(Moose::stringify(id));
  }

  _problem->addVariable(type, var_name, var_params);
}

void
NeutronTransportAction::act()
{

}
