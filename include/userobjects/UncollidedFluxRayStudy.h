#pragma once

#include "RayTracingStudy.h"

// Angular quadrature sets.
#include "GaussAngularQuadrature.h"

// libMesh includes.
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"

class UncollidedFluxRayStudy : public RayTracingStudy
{
public:
  static InputParameters validParams();

  UncollidedFluxRayStudy(const InputParameters & parameters);

protected:
  virtual void generateRays() override;

  static void cartesianToSpherical(const RealVectorValue & direction, Real & mu, Real & omega);
  Real computeSHSource(unsigned int source,
                       const std::vector<std::vector<Real>> & moments,
                       const std::vector<unsigned int> & anisotropy,
                       const RealVectorValue & direction,
                       unsigned int group_index);

  // The index into the data on the ray which contains the source intensity and spatial quadratures,
  // pre-multiplied.
  std::vector<RayDataIndex> _source_spatial_weights;
  // The index into the ray which contains a sign to indicate if the target point and source point
  // are in the same element.
  const RayDataIndex _target_in_element;

  const unsigned int _dim;

  // Radiation transport parameters.
  const unsigned int _num_groups;
  const Real _symmetry_factor;

  // Spatial quadrature rules.
  std::unique_ptr<FEBase> _volume_fe;
  std::unique_ptr<QBase> _q_volume;
  std::unique_ptr<FEBase> _face_fe;
  std::unique_ptr<QBase> _q_face;

  // Settings for source cell -> target cell tracing.
  std::unique_ptr<LegendrePolynomial> _2D_angular_quadrature;
  std::unique_ptr<GaussAngularQuadrature> _3D_angular_quadrature;
  unsigned int _num_dir;

  // Point source information.
  const std::vector<Point> & _point_source_locations;
  const std::vector<std::vector<Real>> & _point_source_moments;
  const std::vector<unsigned int> & _point_source_anisotropy;

  // Surface source information.
  const std::vector<BoundaryName> & _source_boundary_names;
  const std::vector<std::vector<Real>> & _boundary_source_moments;
  const std::vector<unsigned int> & _boundary_source_anisotropy;

  // Volume source information.
  const std::vector<SubdomainName> & _volume_source_blocks;
  const std::vector<std::vector<Real>> & _volume_source_moments;
  const std::vector<unsigned int> & _volume_source_anisotropy;
};
