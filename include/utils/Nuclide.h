#pragma once

#include "NuclearData.h"

namespace NuclearData
{
class Nuclide
{
public:
  static NuclearData::ZAI getZAI(const std::string & nuclide_name);
  static std::string getElement(const std::string & nuclide_name);
  static bool isElement(const std::string & name);
  static std::vector<std::pair<std::string, Real>> getAbundances(const std::string & element_name);
  static Real getAtomicMass(const std::string & nuclide_name);

  Nuclide(const std::string & name, const Real & half_life);

  const std::string & name() const { return _name; }
  const Real & halfLife() const { return _half_life; }
  const Real & decayConst() const { return _decay_const; }
  unsigned int z() const { return _z_a_i._z; }
  unsigned int a() const { return _z_a_i._a; }
  unsigned int i() const { return _z_a_i._i; }

  void addDecay(Decay::Mode mode, const Real & factor, const std::string & target);
  std::string addReaction(Reaction::Mode mode,
                          const Real & factor,
                          const std::string & target,
                          int dz,
                          int da,
                          const Real & q);
  bool addReactionCrossSections(Reaction::Mode mode, const std::vector<Real> & cross_sections);

  const std::vector<Decay::Info> & getDecays() const { return _decays; };
  const std::vector<Reaction::Info> & getReactions() const { return _reactions; };

private:
  const std::string _name;
  const std::string _element_name;
  // A negative half-life indicates that the nuclide is stable. Half-life in
  // seconds.
  const Real _half_life;
  // Decay constant in s^{-1}.
  const Real _decay_const;

  const NuclearData::ZAI _z_a_i; // The Z-A-I index.

  std::vector<Decay::Info> _decays;
  std::vector<Reaction::Info> _reactions;
};
}
