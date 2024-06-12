#pragma once

#include "MooseTypes.h"
#include "GnatBase.h"

#include <vector>

namespace NuclearData
{
class ZAI
{
public:
  ZAI(unsigned int z, unsigned int a, unsigned int i) : _z(z), _a(a), _i(i) {}

  const unsigned int _z; // Atomic number of the nuclide.
  const unsigned int _a; // Mass number of the nuclide.
  const unsigned int _i; // The isomeric state of the nuclide.

  unsigned int index() const { return 10000u * _z + 10u * _a + _i; }
};
}

namespace NuclearData::Decay
{
enum class Mode
{
  IT = 0u,              // Isomeric transition.
  BetaMinus = 1u,       // Beta- (electron emission).
  BetaPlus = 2u,        // Beta+ (positron emission).
  ECBetaPlus = 3u,      // Electron capture resulting in the emission of a positron.
  Alpha = 4u,           // Alpha (He4 emission).
  SF = 5u,              // Spontaneous fission.
  BetaMinusAlpha = 6u,  // Beta- (electron emission) and an alpha emission (He4).
  BetaMinusNeutron = 7u // Beta- (electron emission) and a neutron emission.
};

struct Info
{
  const Mode _mode;             // The decay reaction type.
  const Real _branching_factor; // The decay branching factor.
  const std::string _target;    // The resulting radionuclide.

  Info(Mode mode, const Real & factor, const std::string & target)
    : _mode(mode), _branching_factor(factor), _target(target)
  {
  }
};
}

namespace NuclearData::Source
{
struct Info
{
  const Particletype _particle;         // The decay reaction type.
  std::vector<Real> _p_energies;        // The particle energies
  std::vector<Real> _p_decay_constants; // The decay energy decay constants.
  Real _sum; // Total particle decay constant (regardless of the particle energies).

  Info(Particletype p, const std::vector<Real> & energies_factors) : _particle(p), _sum(0.0)
  {
    for (unsigned int i = 0u; i < energies_factors.size() / 2u; ++i)
    {
      _p_energies.emplace_back(energies_factors[i]);
      _p_decay_constants.emplace_back(energies_factors[(energies_factors.size() / 2u) + i]);

      _sum += energies_factors[(energies_factors.size() / 2u) + i];
    }
  }
};
}

namespace NuclearData::Reaction
{
// TODO: support more reaction types:
// https://github.com/openmc-dev/openmc/blob/develop/openmc/deplete/chain.py
enum class Mode
{
  NGamma = 0u,  // (n, \gamma)
  NProton = 1u, // (n, p)
  NAlpha = 2u,  // (n, \alpha)
  N2N = 3u,     // (n, 2n)
  N3N = 4u,     // (n, 3n)
  N4N = 5u      // (n, 4n)
  // N2Alpha = 6u,   // (n, 2\alpha)
  // NNProton = 7u,  // (n, n + p)
  // N2NProton = 8u, // (n, 2n + p)
  // NDeuteron = 9u, // (n, d)
  // NTriton = 10u,  // (n, t)
  // Fission = 11u   // (n, f)
};

struct Info
{
  const Mode _mode; // The neutron reaction type.

  const Real _branching_factor; // The reaction branching factor.

  const int _dz; // The change in atomic number due to this reaction.
  const int _da; // The change in mass number due to this reaction.

  const Real _q; // The Q-value for this reaction.

  std::vector<Real> _cross_sections; // The group-wise cross-sections for this reaction.
  std::string _target; // The resulting radionuclide (if provided). Required for certain neutron
                       // reactions which result in an nuclear isomer.

  Info(Mode mode,
       const Real & branching_factor,
       const std::string & target,
       int dz,
       int da,
       const Real & q)
    : _mode(mode), _branching_factor(branching_factor), _dz(dz), _da(da), _q(q), _target(target)
  {
  }
};
}
