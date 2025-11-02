#include "LSAngularQuadrature.h"

LSAngularQuadrature::LSAngularQuadrature(unsigned int q_order, MajorAxis axis, ProblemType type)
  : AngularQuadrature(std::move(axis), std::move(type))
{
  if (type == ProblemType::Cartesian1D)
    mooseError("Level-symmetric quadratures only support 2D / 3D problem.");

  if (q_order >= 20)
    mooseError("Level-symmetric quadratures only support up to order 18 at present.");

  if (q_order % 2 == 1)
    mooseError("Odd quadrature orders are not supported by level symmetric quadrature rules.");

  const auto & mus = _ls_mus.at(q_order);
  const auto & weights = _ls_weights.at(q_order);

  // Build the quadrature set in 3D by uncompacting the LS quadrature rule.
  // Adapted from the LS quadrature rule in Moltres:
  // https://github.com/arfc/moltres/blob/devel/src/utils/MoltresUtils.C#L162
  unsigned int w_idx = 0;
  for (unsigned int i = 0; i < (q_order + 5) / 6; ++i)
  {
    for (unsigned int j = i; j < ((q_order / 2 - i) + 1) / 2; ++j)
    {
      unsigned int k = q_order / 2 - 1 - i - j;
      // For each unique combination of x,y,z component magnitudes, permute through all possible
      // permutations (including for positive and negative components).
      addAllOctantPermutations(mus[i], mus[j], mus[k], weights[w_idx]);
      if (i == j && i == k)
      {
        // No extra permutations. Do nothing
      }
      else if (i == j)
      {
        addAllOctantPermutations(mus[k], mus[i], mus[i], weights[w_idx]);
        addAllOctantPermutations(mus[i], mus[k], mus[i], weights[w_idx]);
      }
      else if (j == k)
      {
        addAllOctantPermutations(mus[j], mus[i], mus[j], weights[w_idx]);
        addAllOctantPermutations(mus[j], mus[j], mus[i], weights[w_idx]);
      }
      else
      {
        addAllOctantPermutations(mus[k], mus[i], mus[j], weights[w_idx]);
        addAllOctantPermutations(mus[j], mus[k], mus[i], weights[w_idx]);
        addAllOctantPermutations(mus[k], mus[j], mus[i], weights[w_idx]);
        addAllOctantPermutations(mus[i], mus[k], mus[j], weights[w_idx]);
        addAllOctantPermutations(mus[j], mus[i], mus[k], weights[w_idx]);
      }
      w_idx += 1;
    }
  }

  if (q_order * (q_order + 2) != _quadrature_set_omega.size() && type == ProblemType::Cartesian3D)
    mooseError("Quadrature does not meet the level symmetric quadrature identity!");
  if (q_order * (q_order + 2) / 2u != _quadrature_set_omega.size() && type == ProblemType::Cartesian2D)
    mooseError("Quadrature does not meet the level symmetric quadrature identity!");
}

void
LSAngularQuadrature::addAllOctantPermutations(const Real & mu_1, const Real & mu_2, const Real & mu_3, const Real & w)
{
  _quadrature_set_omega.emplace_back( mu_1,  mu_2,  mu_3);
  _quadrature_set_omega.emplace_back(-mu_1,  mu_2,  mu_3);
  _quadrature_set_omega.emplace_back( mu_1, -mu_2,  mu_3);
  _quadrature_set_omega.emplace_back(-mu_1, -mu_2,  mu_3);
  if (_type == ProblemType::Cartesian3D)
  {
    _quadrature_set_omega.emplace_back( mu_1,  mu_2, -mu_3);
    _quadrature_set_omega.emplace_back(-mu_1,  mu_2, -mu_3);
    _quadrature_set_omega.emplace_back( mu_1, -mu_2, -mu_3);
    _quadrature_set_omega.emplace_back(-mu_1, -mu_2, -mu_3);

    for (unsigned int i = 0u; i < 4u; ++i)
      _quadrature_set_weight.emplace_back(w);
  }

  for (unsigned int i = 0u; i < 4u; ++i)
    _quadrature_set_weight.emplace_back(w);
}

unsigned int
LSAngularQuadrature::totalOrder() const
{
  return _quadrature_set_omega.size();
}

const RealVectorValue &
LSAngularQuadrature::direction(unsigned int n) const
{
  return _quadrature_set_omega[n];
}

const Real &
LSAngularQuadrature::weight(unsigned int n) const
{
  return _quadrature_set_weight[n];
}

const std::vector<RealVectorValue> &
LSAngularQuadrature::getDirections() const
{
  return _quadrature_set_omega;
}

const std::vector<Real> &
LSAngularQuadrature::getWeights() const
{
  return _quadrature_set_weight;
}

// Adapted from the LS quadrature rule in Moltres:
// https://github.com/arfc/moltres/blob/devel/src/utils/MoltresUtils.C#L11-L136
const std::map<unsigned int, std::vector<Real>> LSAngularQuadrature::_ls_mus
{
  { 2u,
    { 0.5773502691896257 }
  },
  { 4u,
    { 0.3500211745815407,
     0.8688903007222011 }
  },
  { 6u,
    { 0.2666354015167047,
      0.6815077265365469,
      0.9261809355174891 }
  },
  { 8u,
    { 0.2182178902359924,
      0.5773502691896257,
      0.7867957924694432,
      0.9511897312113419 }
  },
  { 10u,
    { 0.1893213264780105,
      0.5088817555826189,
      0.6943188875943843,
      0.8397599622366847,
      0.9634909811104685 }
  },
  { 12u,
    { 0.1672324971414912,
      0.4595505230549429,
      0.6280180398523312,
      0.7600175218505538,
      0.8722652169264515,
      0.9716308886607312 }
  },
  { 14u,
    { 0.1519858614610319,
      0.4221569823047970,
      0.5773502691896257,
      0.6988920867759013,
      0.8022262552314120,
      0.8936910988743567,
      0.9766271529257704 }
  },
  { 16u,
    { 0.1389756946126266,
      0.3922930709799660,
      0.5370972569141674,
      0.6504253018068668,
      0.7467480720272295,
      0.8319929644667765,
      0.9092809811978063,
      0.9804955444130666 }
  },
  { 18u,
    { 0.1293481120858299,
      0.3680446084547432,
      0.5041653831086007,
      0.6106624544193616,
      0.7011665515053579,
      0.7812556768832803,
      0.8538655235895108,
      0.9207671975517081,
      0.9831267119754519 }
  }
};

const std::map<unsigned int, std::vector<Real>> LSAngularQuadrature::_ls_weights
{
  { 2u,
    { 1.0 }
  },
  { 4u,
    { 0.33333333333333333 }
  },
  { 6u,
    { 0.17612613086338347,
      0.15720720246994987 }
  },
  { 8u,
    { 0.12098765432098764,
      0.09074074074074068,
      0.09259259259259299 }
  },
  { 10u,
    { 0.08930314798435690,
      0.07252915171236500,
      0.04504376743641028,
      0.05392811448783617 }
  },
  { 12u,
    { 0.07077307572802660,
      0.05587899520720631,
      0.03734279426952707,
      0.05025193141404591,
      0.02586474723779406 }
  },
  { 14u,
    { 0.05799704089700514,
      0.04890079763678386,
      0.02279353424122489,
      0.03941320059498668,
      0.03809908614408365,
      0.02583940764180201,
      0.00826957997290902 }
  },
  { 16u,
    { 0.04899634978321998,
      0.04132620269749813,
      0.02031234525671927,
      0.02654915237867853,
      0.03787631708434387,
      0.01355542870643232,
      0.03259039327572704,
      0.01038401511138594 }
  },
  { 18u,
    { 0.04226679857587908,
      0.03760984374258203,
      0.01227772553689008,
      0.03240535837316794,
      0.00666131638823868,
      0.03120767364508285,
      0.01601092912564359,
      0.02004736607851901,
      0.00012094114211235,
      0.01637415786841503 }
  }
};
