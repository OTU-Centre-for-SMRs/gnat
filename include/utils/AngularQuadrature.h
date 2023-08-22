#pragma once

#include "MooseTypes.h"
#include "libmesh/vector_value.h"

#include "GnatBase.h"

class AngularQuadrature
{
public:
  AngularQuadrature(MajorAxis axis = MajorAxis::X, ProblemType type = ProblemType::Cartesian3D)
    : _axis(axis), _type(type)
  {
  }

  virtual unsigned int totalOrder() const = 0;
  virtual const RealVectorValue & direction(unsigned int n) const = 0;
  virtual const Real & weight(unsigned int n) const = 0;

  virtual const std::vector<RealVectorValue> & getDirections() const = 0;
  virtual const std::vector<Real> & getWeights() const = 0;

  virtual const Real & getPolarRoot(unsigned int n) const = 0;
  virtual const Real & getAzimuthalAngularRoot(unsigned int n) const = 0;

  MajorAxis getAxis() const { return _axis; }
  ProblemType getProblemType() const { return _type; }

protected:
  const MajorAxis _axis;
  const ProblemType _type;
};
