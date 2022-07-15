#pragma once

#include "Action.h"

#include "GnatBase.h"

class NeutronTransportAction : public Action
{
public:
  static InputParameters validParams();

  NeutronTransportAction(const InputParameters & params);

  virtual void act() override;

protected:
  void addVariable(const std::string & var_name);

  std::set<SubdomainID> _subdomain_ids;

  // Number of spectral energy groups and number of discrete ordinates.
  const unsigned int _num_groups;
  const unsigned int _n_l;
  const unsigned int _n_c;

  // Maximum degree of the flux moments requested. By default the transport
  // action prepares the scalar flux and that's it.
  // The Gauss-Legendre-Chebyshev quadrature with half-circle orders
  // n_{g} and n_{c} can fully integrate a spherical harmonics expansion of
  // degree N IF: n_{g} <= N + 1 AND n_{c} <= N + 1.
  const unsigned int _max_eval_anisotropy;

  // Base names of the stored directional angular fluxes and angular flux moments.
  // As an example for both:
  // Psi_{g, n} = _flux_moment_name_g_n.
  // Phi_{g, l, m} = _angular_flux_name_g_l_m.
  const std::string _flux_moment_name;
  const std::string _angular_flux_name;

  // The coordinate system type and dimensionality.
  ProblemType _p_type;

  // Boundary side-sets.
  const std::vector<BoundaryName> _vacuum_side_sets;
  const std::vector<BoundaryName> _source_side_sets;
  const std::vector<BoundaryName> _reflective_side_sets;

}; // class NeutronTransportAction
