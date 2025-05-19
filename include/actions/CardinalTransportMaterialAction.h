#pragma once

#include "GnatBaseAction.h"

class CardinalTransportMaterialAction : public GnatBaseAction
{
public:
  static InputParameters validParams();

  CardinalTransportMaterialAction(const InputParameters & parameters);

  virtual void act() override;

protected:
  void addAuxVariables();
  void addMaterialAuxVar(const std::string & name);

  void addMaterials();

  void addTransfers();
  void addMaterialTransfer(const std::string & name);

  void copyOnRestart();

  const MooseEnum _xs_source;

  const MultiAppName _xs_multi_app;

  std::string _parent_transport_system;
  unsigned int _num_groups;
  MooseEnum _particle;
  MooseEnum _scheme;
  bool _disable_fission;
  const unsigned int _anisotropy;

  bool _is_init;

  std::vector<std::string> _total_var_names;

  std::vector<std::string> _scattering_var_names;

  std::vector<std::string> _nu_fission_var_names;
  std::vector<std::string> _chi_var_names;

  const bool _add_kappa_fission;
  std::vector<std::string> _kf_var_names;

  std::vector<std::string> _inv_v_var_names;

  std::vector<std::string> _diff_var_names;

  std::vector<std::string> _abs_var_names;
};
