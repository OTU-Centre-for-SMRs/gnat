#pragma once

#include "GnatBaseAction.h"

class UncollidedFluxAction : public GnatBaseAction
{
public:
  static InputParameters validParams();

  UncollidedFluxAction(const InputParameters & parameters);

  virtual void act() override;

protected:
  // Different act functions for each uncollided flux scheme.
  void actUncollidedFluxRT();
  void actUncollidedFluxSASF();

  // Member functions to add objects required for the ray traced uncollided flux treatment.
  void addUncollidedRayKernels();
  void addUncollidedRayStudies();
  void addUncollidedRayAuxVars();
  void addUncollidedRayAuxKernels();
  void addUncollidedRayPostProcessors();

  // Member functions required to add objects required for the SASF uncollided flux treatment.
  void addUncollidedSASFVariables();
  void addUncollidedSASFKernels();
  void addUncollidedSASFBCs();
  void addUncollidedSASFAuxVars();
  void addUncollidedSASFAuxKernels();
  void addUncollidedSASFPostProcessors();

  void modifyRTOutputs();

  // Number of spectral energy groups.
  const unsigned int _num_groups;
  // Maximum degree of the flux moments requested.
  const unsigned int _max_eval_anisotropy;
  // The number of flux moments per energy group.
  unsigned int _num_group_moments;

  const std::string & _uncollided_var_base_name;

  // Uncollided flux treatment.
  enum class UncollidedTreatment
  {
    RT = 0u,
    SASF = 1u
  } _uncollided_treatment;

  const bool _conservative_src;

  // Point source moments.
  const std::vector<Point> & _point_source_locations;
  const std::vector<std::vector<Real>> & _point_source_moments;
  const std::vector<unsigned int> & _point_source_anisotropy;

  // Boundary source properties.
  const std::vector<BoundaryName> & _source_side_sets;
  const std::vector<std::vector<Real>> & _boundary_source_moments;
  const std::vector<unsigned int> & _boundary_source_anisotropy;

  // Volumetric sources.
  const std::vector<SubdomainName> & _volumetric_source_blocks;
  const std::vector<std::vector<Real>> & _volumetric_source_moments;
  const std::vector<unsigned int> & _volumetric_source_anisotropy;

  // SASF parameters.
  const BoundaryName & _sasf_near_source_boundary;
  const std::vector<BoundaryName> & _sasf_vacuum_boundaries;
  const std::vector<Real> & _sasf_near_source_cross_sections;
};
