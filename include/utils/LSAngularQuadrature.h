#pragma once

#include <map>

#include "AngularQuadrature.h"

class LSAngularQuadrature : public AngularQuadrature
{
public:
  LSAngularQuadrature(unsigned int q_order, MajorAxis axis = MajorAxis::X, ProblemType type = ProblemType::Cartesian3D);

  unsigned int totalOrder() const override;
  const RealVectorValue & direction(unsigned int n) const override;
  const Real & weight(unsigned int n) const override;
  const std::vector<RealVectorValue> & getDirections() const override;
  const std::vector<Real> & getWeights() const override;

private:
  void addAllOctantPermutations(const Real & mu_1, const Real & mu_2, const Real & mu_3, const Real & w);

  std::vector<RealVectorValue> _quadrature_set_omega;
  std::vector<Real> _quadrature_set_weight;

  static const std::map<unsigned int, std::vector<Real>> _ls_mus;
  static const std::map<unsigned int, std::vector<Real>> _ls_weights;
}; // class LSAngularQuadrature
