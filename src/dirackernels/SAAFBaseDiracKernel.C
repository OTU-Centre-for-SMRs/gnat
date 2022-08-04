#include "SAAFBaseDiracKernel.h"

InputParameters
SAAFBaseDiracKernel::validParams()
{
  auto params = SNBaseDiracKernel::validParams();
  params.addClassDescription("Provides stabalization parameters for SAAF "
                             "Dirac kernels, notably: $h$, $\\tau_{g}$, and "
                             "$\\phi_{j} + \\tau_{g}\\vec{\\nabla}\\phi_{j}"
                             "\\cdot\\hat{\\Omega}$. This kernel does NOT "
                             "implement computeQpResidual() or "
                             "computeQpJacobian().");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "$n$ of the current angular "
                                                    "flux.");

  return params;
}

SAAFBaseDiracKernel::SAAFBaseDiracKernel(const InputParameters & parameters)
  : SNBaseDiracKernel(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _group_index(getParam<unsigned int>("group_index"))
  , _sigma_r_g(getADMaterialProperty<std::vector<Real>>("removal_xs_g"))
  , _saaf_eta(getADMaterialProperty<Real>("saaf_eta"))
  , _saaf_c(getADMaterialProperty<Real>("saaf_c"))
{ }

Real
SAAFBaseDiracKernel::maxVertexSeparation()
{
  const unsigned int n_nodes = _current_elem->n_nodes();

  // Loop over all nodes in the element to find the maximum vertex separation.
  Real separation = std::numeric_limits<Real>::min();
  for (unsigned int i = 0u; i < n_nodes; ++i)
  {
    const auto & point_i = _current_elem->point(i);
    for (unsigned int j = 0u; j < n_nodes; ++j)
    {
      // Ignore the case when the i is equal to the j vertex: they are the same
      // vertex in the element.
      if (i == j)
        continue;

      const auto & point_j = _current_elem->point(j);
      const auto diff = point_j - point_i;
      separation = std::max(separation, diff.norm());
    }
  }

  return separation;
}

Real
SAAFBaseDiracKernel::computeQPTau()
{
  if (_group_index >= _sigma_r_g[_qp].size())
  {
    mooseError("The group index exceeds the number of provided neutron removal "
               "cross-sections.");
  }

  auto h = maxVertexSeparation();
  Real tau = 0.0;
  if (MetaPhysicL::raw_value(_sigma_r_g[_qp][_group_index] * _saaf_c[_qp]) * h
      >= MetaPhysicL::raw_value(_saaf_eta[_qp]))
    tau = 1.0 / MetaPhysicL::raw_value(_sigma_r_g[_qp][_group_index] * _saaf_c[_qp]);
  else
    tau = h / MetaPhysicL::raw_value(_saaf_eta[_qp]);

  return tau;
}

Real
SAAFBaseDiracKernel::computeQPTests()
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  return _test[_i][_qp]
         + computeQPTau() * _grad_test[_i][_qp]
         * _quadrature_set.direction(_ordinate_index);
}
