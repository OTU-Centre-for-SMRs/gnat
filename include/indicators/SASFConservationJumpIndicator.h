#pragma once

#include "InternalSideIndicator.h"

class SASFConservationJumpIndicator : public InternalSideIndicator
{
public:
  static InputParameters validParams();

  SASFConservationJumpIndicator(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  MaterialData & _neighbor_material_data;

  // The source location.
  const Point _source_location;

  // g
  const unsigned int _group_index;
  // Total cross section for the current mesh element.
  const ADMaterialProperty<std::vector<Real>> & _sigma_t_c_g;
  // Total cross section for the neighboring mesh element.
  const ADMaterialProperty<std::vector<Real>> & _sigma_t_n_g;
};
