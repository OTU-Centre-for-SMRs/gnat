#include "ADFVNuclideDiffusion.h"

registerMooseObject("GnatApp", ADFVNuclideDiffusion);

InputParameters
ADFVNuclideDiffusion::validParams()
{
  auto params = FVFluxKernel::validParams();
  params.addClassDescription("Computes the diffusion term for the isotope "
                             "scalar transport equation: "
                             "$( \\vec{\\nabla}\\psi_{j}, "
                             "D_{i}\\vec{\\nabla}N_{i} )_{V}$. This kernel "
                             "expects a diffusion coefficient from the "
                             "material system.");
  params.set<unsigned short>("ghost_layers") = 2;

  return params;
}

ADFVNuclideDiffusion::ADFVNuclideDiffusion(const InputParameters & parameters)
  : FVFluxKernel(parameters),
    _mat_diff(getFunctor<ADReal>("isotope_diff_" +
                                 Moose::stringify(getParam<NonlinearVariableName>("variable"))))
{
}

ADReal
ADFVNuclideDiffusion::computeQpResidual()
{
  // TODO: More sophisticated face interpolation.
  // TODO: skewness correction.
  auto face = makeCDFace(*_face_info);

  return -1 * _mat_diff(face, 0u) * gradUDotNormal(0u);
}
