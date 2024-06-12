#pragma once

#include "ElementPostprocessor.h"

// A class to compute the average value of h_min over the entire domain.
class AverageElementMinSize : public ElementPostprocessor
{
public:
  static InputParameters validParams();

  AverageElementMinSize(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;

  virtual void finalize() override;
  virtual Real getValue() const override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  Real _total_size;
  dof_id_type _elems;
}; // class AverageElementMinSize.
