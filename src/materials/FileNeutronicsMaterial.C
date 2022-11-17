#include "FileNeutronicsMaterial.h"

#include <fstream>
#include <filesystem>

#include "SimpleCSVReader.h"

registerMooseObject("GnatApp", FileNeutronicsMaterial);

InputParameters
FileNeutronicsMaterial::validParams()
{
  auto params = EmptyNeutronicsMaterial::validParams();
  params.addClassDescription(
      "Provides the neutron group velocity ($v_{g}$), "
      "neutron group absorption cross-section "
      "($\\Sigma_{a,g}$), and the scattering cross-"
      "section moments "
      "($\\Sigma_{s, g', g, l}$) for "
      "transport problems. The material properties "
      "are imported from files provided by the user. The user many also provide volumetric neutron "
      "source moments through the input parameter system.");
  params.addRequiredParam<std::string>("file_name", "The file to extract cross-sections from.");
  params.addRequiredParam<std::string>("source_material_id",
                                       "The material ID used by the cross-section source to "
                                       "identify cross-sections for different geometric regions. "
                                       "This parameter will be cast to a different type depending "
                                       "on the cross-section source.");
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
    _file_name(getParam<std::string>("file_name")),
    _source_material_id(getParam<std::string>("source_material_id")),
    _xs_source(getParam<MooseEnum>("cross_section_source").getEnum<CrossSectionSource>()),
    _anisotropy(0u),
    _max_moments(0u),
    _source_moments(getParam<std::vector<Real>>("group_source")),
    _source_anisotropy(getParam<unsigned int>("source_anisotropy")),
    _max_source_moments(0u),
    _has_volumetric_source(false)
{
  // Open the descriptor file to figure out where the cross-sections are stored and what the files
  // are named.
  std::ifstream descriptor_file(_file_name);
  std::string line;
  if (descriptor_file.is_open())
  {
    std::getline(descriptor_file, line);
    if (_xs_source == CrossSectionSource::Detect)
    {
      // First line should always be the cross-section source. Use that to determine the
      // properties.
      if (line == "openmc")
        _xs_source = CrossSectionSource::OpenMC;
      else if (line == "gnat")
        _xs_source = CrossSectionSource::Gnat;
      else
        mooseError("Unknown cross-section source.");
    }

    // Read each line. Assume that anything after the keyword is the relative path of the
    // specific cross-section file to the descriptor file.
    while (std::getline(descriptor_file, line))
    {
      auto pos = line.find("InvVelocity: ");
      if (pos != std::string::npos)
      {
        parseProperty(PropertyType::InvVelocity,
                      line.substr(pos + std::string("InvVelocity: ").size()));
        continue;
      }

      pos = line.find("SigmaR: ");
      if (pos != std::string::npos)
      {
        parseProperty(PropertyType::SigmaR, line.substr(pos + std::string("SigmaR: ").size()));
        continue;
      }

      pos = line.find("SigmaA: ");
      if (pos != std::string::npos)
      {
        parseProperty(PropertyType::SigmaA, line.substr(pos + std::string("SigmaA: ").size()));
        continue;
      }

      pos = line.find("SigmaS: ");
      if (pos != std::string::npos)
      {
        parseProperty(PropertyType::SigmaS, line.substr(pos + std::string("SigmaS: ").size()));
        continue;
      }

      pos = line.find("SigmaSMatrix: ");
      if (pos != std::string::npos)
      {
        parseProperty(PropertyType::SigmaSMatrix,
                      line.substr(pos + std::string("SigmaSMatrix: ").size()));
        continue;
      }
    }

    descriptor_file.close();
  }
  else
    mooseError("Failed to open file " + _file_name + " to read material properties.");

  // Resize the properties and initialize with 0.0.
  _sigma_r_g.resize(_num_groups, 0.0);
  _sigma_s_g_prime_g_l.resize(_num_groups * _num_groups * (_anisotropy + 1u), 0.0);

  // Convert per-nuclide cross-sections to bulk material cross-sections.
  bool first_inv_v = true;
  for (auto & [nuclide, data] : _material_properties)
  {
    // InvVelocity is the exception, not the rule. Properties are identical for each isotope.
    if (first_inv_v)
    {
      if (data._inv_v_g.size() < _num_groups)
      {
        mooseError("Number of group-wise velocities is smaller then the number of groups declared "
                   "in the simulation properties.");
      }
      else if (data._inv_v_g.size() > _num_groups)
        mooseWarning(
            "Number of group-wise velocities is greater then then number of groups declared in the "
            "simulation properties. The number of group-wise velocities will be truncated.");

      for (auto & velocity : data._inv_v_g)
        _inv_v_g.emplace_back(velocity);

      first_inv_v = false;
    }

    // SigmaR.
    {
      if (data._sigma_r_g.size() < _num_groups)
      {
        mooseError("Number of group-wise removal cross-sections is smaller then the number of "
                   "groups declared in the simulation properties.");
      }
      else if (data._sigma_r_g.size() > _num_groups)
        mooseWarning(
            "Number of group-wise removal cross-sections is greater then then number of groups "
            "declared in the simulation properties. The number of group-wise removal "
            "cross-sections will be truncated.");

      for (unsigned int i = 0u; i < _sigma_r_g.size(); ++i)
        _sigma_r_g[i] += data._sigma_r_g[i];
    }

    // SigmaSMatrix.
    {
      if (data._sigma_s_g_prime_g_l.size() < _max_moments)
      {
        mooseError("Number of scattering matrix components is smaller then the number allowed by "
                   "the used simulation properties.");
      }
      else if (data._sigma_s_g_prime_g_l.size() > _max_moments)
        mooseWarning("Number of scattering matrix components is greater then the number required "
                     "for the simulation. The matrix may be out of order.");

      for (unsigned int i = 0u; i < _sigma_s_g_prime_g_l.size(); ++i)
        _sigma_s_g_prime_g_l[i] += data._sigma_s_g_prime_g_l[i];
    }
  }

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
}

void
FileNeutronicsMaterial::parseProperty(const PropertyType & type, const std::string & property_file)
{
  switch (_xs_source)
  {
    case CrossSectionSource::Gnat:
      parseGnatProperty(type, property_file);
      break;

    case CrossSectionSource::OpenMC:
      parseOpenMCProperty(type, property_file);
      break;

    default:
      break;
  }
}

void
FileNeutronicsMaterial::parseGnatProperty(const PropertyType & type,
                                          const std::string & property_file)
{
  mooseError("Parsing of cross-sections generated by Gnat is currently not supported.");

  switch (type)
  {
    case PropertyType::InvVelocity:
      break;

    case PropertyType::SigmaR:
      break;

    case PropertyType::SigmaA:
      break;

    case PropertyType::SigmaS:
      break;

    case PropertyType::SigmaSMatrix:
      break;

    default:
      break;
  }
}

void
FileNeutronicsMaterial::parseOpenMCProperty(const PropertyType & type,
                                            const std::string & property_file)
{
  // Compute the actual path of the property file.
  auto dir = std::filesystem::path(_file_name).parent_path();
  dir /= std::filesystem::path(property_file);

  // Parse!
  SimpleCSVReader reader(dir.string());
  if (reader.read())
  {
    // Loop over the cell ids to find ones that match the provided id in the material.
    int min_row = -1;
    int max_row = -1;
    {
      auto & cell_ids = reader.getColumn("cell");
      // Start by finding the minimum row.
      for (unsigned int row = 0u; row < cell_ids.size(); ++row)
      {
        if (cell_ids[row] == _source_material_id)
        {
          min_row = static_cast<int>(row);
          break;
        }
      }

      if (min_row == -1)
        mooseError("OpenMC cell with label " + _source_material_id + " does not exist in " +
                   property_file);

      // Next find the maximum row (taking advantage of how OpenMC groups cells together).
      for (unsigned int row = static_cast<unsigned int>(min_row) + 1; row < cell_ids.size(); ++row)
      {
        if (cell_ids[row - 1u] == _source_material_id && cell_ids[row] != _source_material_id)
        {
          max_row = static_cast<int>(row - 1u);
          break;
        }
      }

      if (min_row != -1 && max_row == -1)
        max_row = static_cast<int>(cell_ids.size() - 1u);
    }

    // Loop over the data to determine the number of isotopes.
    int num_groups = -1;
    unsigned int num_isotopes = 0;
    {
      std::string first_nuclide = reader.getEntry("nuclide", static_cast<unsigned int>(min_row));
      for (unsigned int row = static_cast<unsigned int>(min_row);
           row <= static_cast<unsigned int>(max_row);
           ++row)
      {
        num_groups = std::max(num_groups, std::stoi(reader.getEntry("group in", row)));

        if (first_nuclide == reader.getEntry("nuclide", row) &&
            row != static_cast<unsigned int>(min_row) && first_nuclide != "")
        {
          num_isotopes = row - static_cast<unsigned int>(min_row);
          first_nuclide = "";
        }
      }
    }

    // Loop over the data and add nuclides to the unordered map of data.
    {
      auto & nuclides = reader.getColumn("nuclide");
      for (unsigned int row = static_cast<unsigned int>(min_row);
           row < static_cast<unsigned int>(min_row) + num_isotopes;
           ++row)
        _material_properties.try_emplace(nuclides[row]);
    }

    // Loop over all of the rows of data and parse accordingly.
    // OpenMC orders energy groups from the largest to smallest. As an example: group 1 would be
    // fast, group 2 would be thermal.
    // The case switch-case works with the different nuclear properties.
    auto & values = reader.getColumn("mean");
    bool has_moments = reader.has("legendre");
    switch (type)
    {
      case PropertyType::InvVelocity:
        for (unsigned int row = static_cast<unsigned int>(min_row);
             row <= static_cast<unsigned int>(max_row);
             ++row)
        {
          _material_properties[reader.getEntry("nuclide", row)]._inv_v_g.emplace_back(
              std::stod(values[row]));
        }
        break;

      case PropertyType::SigmaR:
        for (unsigned int row = static_cast<unsigned int>(min_row);
             row <= static_cast<unsigned int>(max_row);
             ++row)
        {
          _material_properties[reader.getEntry("nuclide", row)]._sigma_r_g.emplace_back(
              std::stod(values[row]));
        }
        break;

      case PropertyType::SigmaA:
        for (unsigned int row = static_cast<unsigned int>(min_row);
             row <= static_cast<unsigned int>(max_row);
             ++row)
        {
          _material_properties[reader.getEntry("nuclide", row)]._sigma_a_g.emplace_back(
              std::stod(values[row]));
        }
        break;

      case PropertyType::SigmaS:
        for (unsigned int row = static_cast<unsigned int>(min_row);
             row <= static_cast<unsigned int>(max_row);
             ++row)
        {
          _material_properties[reader.getEntry("nuclide", row)]._sigma_s_g.emplace_back(
              std::stod(values[row]));
        }
        break;

      case PropertyType::SigmaSMatrix:
        // Parse the material properties and detect the maximum anisotropy.
        for (unsigned int row = static_cast<unsigned int>(min_row);
             row <= static_cast<unsigned int>(max_row);
             ++row)
        {
          _material_properties[reader.getEntry("nuclide", row)]._sigma_s_g_prime_g_l.emplace_back(
              std::stod(values[row]));

          if (has_moments)
          {
            _anisotropy = std::max(
                _anisotropy,
                static_cast<unsigned int>(std::stoi(reader.getEntry("legendre", row).substr(1u))));
          }
        }
        _max_moments = (_anisotropy + 1u) * _num_groups * _num_groups;
        break;

      default:
        break;
    }
  }
  else
    mooseError("Failed to open cross-section file: " + dir.string());
}

void
FileNeutronicsMaterial::computeQpProperties()
{
  EmptyNeutronicsMaterial::computeQpProperties();

  _mat_inv_v_g[_qp].resize(_num_groups, 0.0);
  _mat_sigma_t_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int i = 0; i < _num_groups; ++i)
  {
    _mat_inv_v_g[_qp][i] = _inv_v_g[i];
    _mat_sigma_t_g[_qp][i] = _sigma_r_g[i];
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
}
