#pragma once

#include "ThreadedGeneralUserObject.h"

#include "AngularQuadrature.h"

// A userobject which provides data regarding the angular quadrature set to discrete ordinate
// transport schemes.
class AQProvider : public ThreadedGeneralUserObject
{
public:
  static InputParameters validParams();

  AQProvider(const InputParameters & parameters);

  virtual void execute() final {}
  virtual void initialize() final {}
  virtual void finalize() final {}

  virtual void threadJoin(const UserObject &) final {}
  virtual void subdomainSetup() final {}

  unsigned int totalOrder() const { return _aq->totalOrder(); }
  const RealVectorValue & direction(unsigned int n) const { return _aq->direction(n); }
  const Real & weight(unsigned int n) const { return _aq->weight(n); }
  const std::vector<RealVectorValue> & getDirections() const { return _aq->getDirections(); }
  const std::vector<Real> & getWeights() const { return _aq->getWeights(); }

  const Real & getPolarRoot(unsigned int n) const { return _aq->getPolarRoot(n); }
  const Real & getAzimuthalAngularRoot(unsigned int n) const
  {
    return _aq->getAzimuthalAngularRoot(n);
  }

  MajorAxis getAxis() const { return _aq->getAxis(); }
  ProblemType getProblemType() const { return _aq->getProblemType(); }

protected:
  enum class AQType
  {
    GaussChebyshev = 0u
  } _aq_type;

  std::unique_ptr<AngularQuadrature> _aq;
}; // class ThreadedGeneralUserObject
