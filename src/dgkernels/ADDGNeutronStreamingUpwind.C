#include "ADDGNeutronStreamingUpwind.h"

registerMooseObject("GnatApp", ADDGNeutronStreamingUpwind);

InputParameters
ADDGNeutronStreamingUpwind::validParams()
{
  auto params = ADDGKernel::validParams();
  // TODO: Fix the weak form in the class description.
  params.addClassDescription("Computes the discontinuous face term for the "
                             "streaming operator in the discrete ordinates "
                             "neutron transport equation. The weak form is "
                             "given by "
                             "$-(\\nabla \\psi_{j}\\cdot\\vec{\\Omega}, "
                             "\\Psi_{g, n}^{k})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "$n$ of the current angular "
                                                    "flux.");

  return params;
}

ADDGNeutronStreamingUpwind::ADDGNeutronStreamingUpwind(const InputParameters & parameters)
  : ADDGKernel(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _directions(getADMaterialProperty<std::vector<RealVectorValue>>("directions"))
{ }

ADReal
ADDGNeutronStreamingUpwind::computeQpResidual(Moose::DGResidualType type)
{
  if (_ordinate_index >= _directions[_qp].size())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  ADReal res = 0.0;
  ADReal n_dot_omega = _normals[_qp] * _directions[_qp][_ordinate_index];

  if (n_dot_omega >= 0.0)
  {
    // Neighbor is upwind, current element is downwind.
    switch (type)
    {
      // Add the continuous residual contribution from the current element.
      // Addition because n_dot_omega is positive.
      case Moose::Element:
        res += n_dot_omega * _u[_qp] * _test[_i][_qp];
        break;

      // Remove the continuous residual contribution from the neighboring element.
      // Subtraction because n_dot_omega is positive.
      case Moose::Neighbor:
        res -= n_dot_omega * _u[_qp] * _test_neighbor[_i][_qp];
        break;
    }
  }
  else
  {
    // Neighbor is downwind, current element is upwind.
    switch (type)
    {
      // Remove the continuous residual contribution from the current element.
      // Addition because n_dot_omega is negative (net result is a removal).
      case Moose::Element:
        res += n_dot_omega * _u_neighbor[_qp] * _test[_i][_qp];
        break;

      // Add the continuous residual contribution from the neighboring element.
      // Subtraction because n_dot_omega is negative (net result is an addition).
      case Moose::Neighbor:
        res -= n_dot_omega * _u_neighbor[_qp] * _test_neighbor[_i][_qp];
        break;
    }
  }

  return res;
}
