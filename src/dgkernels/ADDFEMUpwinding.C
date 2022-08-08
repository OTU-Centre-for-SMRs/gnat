#include "ADDFEMUpwinding.h"

registerMooseObject("GnatApp", ADDFEMUpwinding);

InputParameters
ADDFEMUpwinding::validParams()
{
  auto params = ADDGKernel::validParams();
  params.addClassDescription("Computes the discontinuous face term for the "
                             "streaming operator in the discrete ordinates "
                             "neutron transport equation using upwinding for "
                             "stabalization. This kernel should not be exposed "
                             "to the user, instead being enabled through a "
                             "transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("n_l",
                                                    "n_l > 0",
                                                    "Order of the polar Gauss-"
                                                    "Legendre quadrature set.");
  params.addRequiredRangeCheckedParam<unsigned int>("n_c",
                                                    "n_c > 0",
                                                    "Order of the azimuthal "
                                                    "Gauss-Chebyshev "
                                                    "quadrature set.");
  params.addParam<MooseEnum>("major_axis", MooseEnum("x y z", "x"),
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis. This parameter is ignored for 1D and 2D "
                             "problems.");
  params.addRequiredParam<MooseEnum>("dimensionality",
                                     MooseEnum("1D_cartesian 2D_cartesian 3D_cartesian"),
                                     "Dimensionality and the coordinate system of the "
                                     "problem.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "$n$ of the current angular "
                                                    "flux.");

  return params;
}

ADDFEMUpwinding::ADDFEMUpwinding(const InputParameters & parameters)
  : ADDGKernel(parameters)
  , _quadrature_set(getParam<unsigned int>("n_c"),
                    getParam<unsigned int>("n_l"),
                    getParam<MooseEnum>("major_axis").getEnum<MajorAxis>(),
                    getParam<MooseEnum>("dimensionality").getEnum<ProblemType>())
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
{ }

ADReal
ADDFEMUpwinding::computeQpResidual(Moose::DGResidualType type)
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  ADReal res = 0;
  ADReal n_dot_omega = _normals[_qp] * _quadrature_set.direction(_ordinate_index);

  switch (type)
  {
    case Moose::Element:
      if (n_dot_omega >= 0)
        res += n_dot_omega * _u[_qp] * _test[_i][_qp];
      else
        res += n_dot_omega * _u_neighbor[_qp] * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      if (n_dot_omega >= 0)
        res -= n_dot_omega * _u[_qp] * _test_neighbor[_i][_qp];
      else
        res -= n_dot_omega * _u_neighbor[_qp] * _test_neighbor[_i][_qp];
      break;
  }
  
  return res;
}
