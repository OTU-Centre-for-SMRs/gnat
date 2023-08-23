#include "FileNeutronicsMaterial.h"

#include <fstream>
#include <filesystem>
#include "pugixml.h"

registerMooseObject("GnatApp", FileNeutronicsMaterial);

// #define DEBUG_CROSS_SECTIONS

InputParameters
FileNeutronicsMaterial::validParams()
{
  auto params = EmptyNeutronicsMaterial::validParams();
  params.addClassDescription(
      "Provides the neutron group velocity ($v_{g}$), "
      "neutron group total cross-section "
      "($\\Sigma_{t,g}$), and the scattering cross-"
      "section moments "
      "($\\Sigma_{s, g', g, l}$) for "
      "transport problems. The material properties "
      "are imported from files provided by the user. The user many also provide volumetric neutron "
      "source moments through the input parameter system.");
  params.addRequiredParam<std::string>("file_name",
                                       "The file to extract cross-sections from. The path of the "
                                       "file should be relative to the input deck (.i) file.");
  params.addRequiredParam<std::string>("source_material_id",
                                       "The material ID used by the cross-section source to "
                                       "identify cross-sections for different geometric regions.");
  params.addParam<MooseEnum>("cross_section_source",
                             MooseEnum("detect gnat openmc", "detect"),
                             "The source which generated the cross-sections. Used for "
                             "parsing different cross-section formats. If set to 'detect' the "
                             "material will attempt to determine the source and parse "
                             "accordingly.");

  // Optional parameters for a volumetric particle source.
  params.addParam<std::vector<Real>>("group_source",
                                     "The external source moments for all energy groups.");
  params.addParam<unsigned int>(
      "source_anisotropy", 0u, "The external source anisotropy of the medium.");

  return params;
}

FileNeutronicsMaterial::FileNeutronicsMaterial(const InputParameters & parameters)
  : EmptyNeutronicsMaterial(parameters),
    _anisotropy(0u),
    _max_moments(0u),
    _source_moments(getParam<std::vector<Real>>("group_source")),
    _source_anisotropy(getParam<unsigned int>("source_anisotropy")),
    _max_source_moments(0u),
    _has_volumetric_source(false),
    _file_name((std::filesystem::path(_app.getInputFileName()).parent_path() /
                std::filesystem::path(getParam<std::string>("file_name")))
                   .string()),
    _source_material_id(getParam<std::string>("source_material_id"))
{
  // Parse the cross-section XML file.
  parseXMLMacroXS();

  // Validate the resulting cross-sections.
  // TODO: validate Chi and NuSigmaF.
  if (_inv_v_g.size() != _num_groups)
    mooseError("The inverse velocity data failed to parse properly.");
  if (_sigma_t_g.size() != _num_groups)
    mooseError("The total cross-section data failed to parse properly.");
  if (_sigma_a_g.size() != _num_groups)
    mooseError("The absorption cross-section data failed to parse properly.");
  if (_sigma_s_g_prime_g_l.size() != _max_moments)
    mooseError("The scattering matrix cross-section data failed to parse properly.");

  // Resize the properties and initialize with 0.0.
  _sigma_s_g_matrix.resize(_num_groups, 0.0);
  _sigma_s_g_g.resize(_num_groups, 0.0);

  _has_volumetric_source = _source_moments.size() > 0u;
  if (_has_volumetric_source)
  {
    switch (_mesh.dimension())
    {
      case 1u:
        _max_source_moments = (_source_anisotropy + 1u);
        _max_source_moments *= _num_groups;
        break;

      case 2u:
        _max_source_moments = (_source_anisotropy + 1u) * (_source_anisotropy + 2u) / 2u;
        _max_source_moments *= _num_groups;
        break;

      case 3u:
        _max_source_moments = (_source_anisotropy + 1u) * (_source_anisotropy + 1u);
        _max_source_moments *= _num_groups;
        break;

      default:
        mooseError("Unknown mesh dimensionality.");
        break;
    }

    // Warn the user if more parameters have been provided than required.
    if (_source_moments.size() > _max_source_moments)
    {
      mooseWarning("More source moments have been provided than possibly "
                   "supported with the given maximum source anisotropy and "
                   "number of groups. The vector will be truncated.");
    }

    // Error if the user did not provide enough parameters.
    if (_source_moments.size() < _max_source_moments)
      mooseError("Not enough source moments have been provided.");
  }

  // Required for neutron diffusion coefficients and removal cross-sections.
  // Compute the out-scattering cross-section. This is the sum of the 0'th
  // moments of the scattering cross-sections from the current group into all
  // other groups (excluding the current group).
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    _sigma_s_g_g[g] =
        _sigma_s_g_prime_g_l[g * _num_groups * (_anisotropy + 1u) + g * (_anisotropy + 1u)];

    for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
      _sigma_s_g_matrix[g] +=
          _sigma_s_g_prime_g_l[g * _num_groups * (_anisotropy + 1u) + g_prime * (_anisotropy + 1u)];
  }

  // Recompute the neutron diffusion coefficients to account for scattering (transport
  // approximation).
  // Sum the first order Legendre out-scattering cross-sections for all groups.
  if (_diffusion_g.size() == 0u && _is_diffusion)
  {
    std::vector<Real> _sigma_s_g_prime_g_1;
    _sigma_s_g_prime_g_1.resize(_num_groups, 0.0);

    _diffusion_g.resize(_num_groups, 0.0);
    bool warning = false;
    for (unsigned int g = 0u; g < _diffusion_g.size(); ++g)
    {
      if (3.0 * (_sigma_t_g[g] - _sigma_s_g_prime_g_1[g]) < 10.0 * libMesh::TOLERANCE)
      {
        _diffusion_g[g] = 1.0 / (10.0 * libMesh::TOLERANCE);
        warning = true;
      }
      else
        _diffusion_g[g] = 1.0 / (3.0 * (_sigma_t_g[g] - _sigma_s_g_prime_g_1[g]));
    }
    if (warning)
      mooseWarning(
          "3.0 * (_sigma_t_g[g] - _sigma_s_g_prime_g_1[g]) < 10.0 * libMesh::TOLERANCE for the "
          "provided cross-section(s). Using a diffusion coefficient of 1 / "
          "10.0 * libMesh::TOLERANCE for those values.");
  }

#ifdef DEBUG_CROSS_SECTIONS
  // Output material properties for verification.
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    _console << COLOR_GREEN << "Group " << g << " cross-sections for " << _name << ":\n"
             << COLOR_DEFAULT;
    _console << "_sigma_t_g              = " << _sigma_t_g[g] << "\n";
    _console << "_sigma_a_g              = " << _sigma_a_g[g] << "\n";
    _console << "_sigma_s_g              = " << _sigma_s_g[g] << "\n";
    _console << "_sigma_s_g + _sigma_a_g = " << _sigma_s_g[g] + _sigma_a_g[g] << "\n";
    _console << "_sigma_s_g_g            = " << _sigma_s_g_g[g] << "\n";
    _console << "IsConservative:           " << std::boolalpha
             << ((_sigma_s_g[g] + _sigma_a_g[g]) <= _sigma_t_g[g]);

    if (_has_fission)
    {
      _console << "\n_nu_sigma_f_g       = " << _nu_sigma_f_g[g] << "\n";
      _console << "_chi_f_g              = " << _chi_f_g[g];
    }
    _console << std::endl;
  }
#endif
}

void
FileNeutronicsMaterial::parseToVector(const std::string & string_rep, std::vector<Real> & real_rep)
{
  auto current_delim_pos = string_rep.find(" ");
  std::size_t offset = 0u;
  auto previous_delim_pos = 0u;
  while (current_delim_pos != std::string::npos)
  {
    real_rep.emplace_back(std::stod(string_rep.substr(
        previous_delim_pos + offset, current_delim_pos - (previous_delim_pos + offset))));
    previous_delim_pos = current_delim_pos;
    current_delim_pos = string_rep.find(" ", current_delim_pos + 1);
    offset = 1u;
  }
  real_rep.emplace_back(std::stod(string_rep.substr(previous_delim_pos)));
}

// TODO: parse Chi and NuSigmaF.
void
FileNeutronicsMaterial::parseXMLMacroXS()
{
  pugi::xml_document doc;
  if (!doc.load_file(_file_name.c_str()))
    mooseError("Failed to load the following microscopic cross-section file: " + _file_name);

  auto all_materials = doc.child("macroscopic_cross_sections");

  // Parse the energy units.
  if (std::string(all_materials.attribute("energy_units").as_string()) == "eV")
    _energy_units = EnergyUnits::eV;
  else if (std::string(all_materials.attribute("energy_units").as_string()) == "keV")
    _energy_units = EnergyUnits::keV;
  else if (std::string(all_materials.attribute("energy_units").as_string()) == "MeV")
    _energy_units = EnergyUnits::MeV;
  else
    _console << COLOR_YELLOW << "Unsupported units of energy: "
             << all_materials.attribute("energy_units").as_string() << ". Defauling to eV.\n"
             << COLOR_DEFAULT;

  // Parse the cross-section units. cm^{-1} is currently the only supported unit.
  if (std::string(all_materials.attribute("xs_units").as_string()) == "cm^-1")
    _xs_units = XSUnits::InvCm;
  else
    _console << COLOR_YELLOW << "Cross-section units cannot be detected. Defauling to cm^-1.\n"
             << COLOR_DEFAULT;

  // Parse the number of groups and group boundaries.
  unsigned int num_groups_in_file = all_materials.attribute("num_groups").as_uint(0u);
  if (num_groups_in_file != _num_groups)
    mooseError("The number of groups provided (" + Moose::stringify(num_groups_in_file) +
               ") is not consistent with the number of groups in the simulation (" +
               Moose::stringify(_num_groups) + ").");
  parseToVector(std::string(all_materials.attribute("group_bounds").as_string()), _group_bounds);

  // Attempt to find the material by ID the user provided in the input syntax.
  for (auto & material_node : all_materials)
  {
    // We have a match.
    if (std::string(material_node.attribute("id").as_string()) == _source_material_id &&
        std::string(material_node.name()) == "domain")
    {
      _anisotropy = material_node.attribute("num_legendre").as_uint(0u);
      _max_moments = _num_groups * _num_groups * (_anisotropy + 1u);

      for (auto & xs_node : material_node)
      {
        if (std::string(xs_node.attribute("type").as_string()) == "total")
          parseToVector(std::string(xs_node.attribute("mgxs").as_string()), _sigma_t_g);

        if (std::string(xs_node.attribute("type").as_string()) == "absorption")
          parseToVector(std::string(xs_node.attribute("mgxs").as_string()), _sigma_a_g);

        if (std::string(xs_node.attribute("type").as_string()) == "scatter")
          parseToVector(std::string(xs_node.attribute("mgxs").as_string()), _sigma_s_g_prime_g_l);

        if (std::string(xs_node.attribute("type").as_string()) == "inverse-velocity")
          parseToVector(std::string(xs_node.attribute("mgxs").as_string()), _inv_v_g);
      }

      // Can stop searching after all cross-section data is parsed.
      return;
    }
  }

  mooseError("Failed to find a material with the ID " + _source_material_id +
             " in the provided cross-section file.");
}

void
FileNeutronicsMaterial::computeQpProperties()
{
  EmptyNeutronicsMaterial::computeQpProperties();

  // SAAF tau.
  if (_is_saaf)
  {
    (*_mat_saaf_tau)[_qp].resize(_num_groups, 0.0);

    auto h = _current_elem->hmin();
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      if (_sigma_t_g[g] * _saaf_c * h >= _saaf_eta)
        (*_mat_saaf_tau)[_qp][g] = 1.0 / (_sigma_t_g[g] * _saaf_c);
      else
        (*_mat_saaf_tau)[_qp][g] = h / _saaf_eta;
    }
  }

  _mat_sigma_t_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int i = 0; i < _num_groups; ++i)
    _mat_sigma_t_g[_qp][i] = _sigma_t_g[i];

  if (_is_diffusion)
  {
    (*_mat_sigma_r_g)[_qp].resize(_num_groups, 0.0);
    (*_mat_diffusion_g)[_qp].resize(_num_groups, 0.0);
    for (unsigned int i = 0; i < _num_groups; ++i)
    {
      // Have to sum the absorption and out-scattering cross-section to form the
      // total cross-section. This sums all g -> g_prime scattering cross-sections and then
      // subtracts the g -> g cross-section.
      (*_mat_sigma_r_g)[_qp][i] = (_sigma_t_g[i] - _sigma_s_g_g[i]);
      (*_mat_diffusion_g)[_qp][i] = _diffusion_g[i];
    }
  }

  // Particle speeds.
  _mat_inv_v_g[_qp].resize(_num_groups, 0.0);
  if (_particle == Particletype::Neutron)
  {
    for (unsigned int i = 0; i < _num_groups; ++i)
      _mat_inv_v_g[_qp][i] = _inv_v_g[i];
  }
  if (_particle == Particletype::GammaPhoton)
  {
    for (unsigned int i = 0; i < _num_groups; ++i)
      _mat_inv_v_g[_qp][i] = _inv_c_cm;
  }

  // Scattering moments and anisotropy.
  _mat_anisotropy[_qp] = _anisotropy;
  _mat_sigma_s_g_prime_g_l[_qp].resize(_max_moments, 0.0);
  for (unsigned int i = 0u; i < _max_moments; ++i)
    _mat_sigma_s_g_prime_g_l[_qp][i] = _sigma_s_g_prime_g_l[i];

  // Source moments and anisotropy.
  if (_has_volumetric_source)
  {
    _mat_src_anisotropy[_qp] = _source_anisotropy;
    _mat_source_moments[_qp].resize(_max_source_moments, 0.0);
    for (unsigned int i = 0u; i < _max_source_moments; ++i)
      _mat_source_moments[_qp][i] = _source_moments[i];
  }

  // Fission production cross-sections and spectra.
  if (_has_fission)
  {
    (*_mat_nu_sigma_f_g)[_qp].resize(_num_groups, 0.0);
    (*_mat_chi_f_g)[_qp].resize(_num_groups, 0.0);
    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      (*_mat_nu_sigma_f_g)[_qp][g] = _nu_sigma_f_g[g];
      (*_mat_chi_f_g)[_qp][g] = _chi_f_g[g];
    }
  }
}
