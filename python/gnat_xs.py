import os
import warnings

import numpy as np
import lxml.etree as ET

import openmc.mgxs as mgxs

#---------------------------------------------------------------------------------------------------------------------------
# Get a depletion list from a collection of domains.
#---------------------------------------------------------------------------------------------------------------------------
def get_depletion_mgxs_list(domain=None, domain_type=None, energy_groups=None, name='', num_polar=1, num_azimuthal=1):
  depletion_list = []
  for reaction_type in ('(n,2n)', '(n,3n)', '(n,4n)', '(n,gamma)', '(n,p)', '(n,a)'):
    depletion_list.append(mgxs.ArbitraryXS(reaction_type, domain, domain_type, energy_groups, True, name, num_polar, num_azimuthal))

  return depletion_list
#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Get a GNAT compatible MGXS list from a collection of domains.
#---------------------------------------------------------------------------------------------------------------------------
def get_mgxs_list(domain=None, domain_type=None, energy_groups=None, legendre_order=0, correction=None, name='', num_polar=1, num_azimuthal=1, nu=False, fission=False, diffusion=False, kappa=False):
  xs_list = []

  xs_list.append(mgxs.TotalXS(domain, domain_type, energy_groups, False, name, num_polar, num_azimuthal))

  xs_list.append(mgxs.AbsorptionXS(domain, domain_type, energy_groups, False, name, num_polar, num_azimuthal))

  scatter = mgxs.ScatterMatrixXS(domain, domain_type, energy_groups, False, name, num_polar, num_azimuthal, nu)
  scatter.formulation = 'consistent'
  scatter.scatter_format = 'legendre'
  scatter.correction = correction
  scatter.legendre_order = legendre_order
  xs_list.append(scatter)

  xs_list.append(mgxs.InverseVelocity(domain, domain_type, energy_groups, False, name, num_polar, num_azimuthal))

  if fission:
    xs_list.append(mgxs.FissionXS(domain, domain_type, energy_groups, nu=True))
    xs_list.append(mgxs.Chi(domain, domain_type, energy_groups))

  if diffusion:
    xs_list.append(mgxs.DiffusionCoefficient(domain, domain_type, energy_groups, nu, False, name, num_polar, num_azimuthal))

  if kappa:
    xs_list.append(mgxs.KappaFissionXS(domain, domain_type, energy_groups, False, name, num_polar, num_azimuthal))

  return xs_list

# Apply tallies to a tallies file.
def apply_tallies(tallies, xs_list):
  for xs in xs_list:
    tallies += xs.tallies.values()

  return tallies

#---------------------------------------------------------------------------------------------------------------------------
# Process results from an OpenMC cross-section library.
#---------------------------------------------------------------------------------------------------------------------------
def process_all_macro_xs_results_from_lib_xml(state_point, lib, file_name, directory, condense = False,
                                              coarse_group_structure = mgxs.EnergyGroups(np.array([0.0, 20.0e6]))):
  lib.load_from_statepoint(state_point)

  xs_string_list = []

  # Check for mandatory cross-sections. Total, absorption, and a scattering matrix.
  if lib.mgxs_types.count('total') == 0:
    raise Exception("This library doesn't contain an total cross-section.")
  else:
    xs_string_list.append('total')

  if lib.mgxs_types.count('absorption') == 0:
    raise Exception("This library doesn't contain an absorption cross-section.")
  else:
    xs_string_list.append('absorption')

  if lib.mgxs_types.count('inverse-velocity') > 0:
    xs_string_list.append('inverse-velocity')

  has_scattering_xs = False
  if lib.mgxs_types.count('nu-scatter matrix') > 0:
    has_scattering_xs = True
    xs_string_list.append('nu-scatter matrix')

  if lib.mgxs_types.count('scatter matrix') > 0 and not has_scattering_xs:
    has_scattering_xs = True
    xs_string_list.append('scatter matrix')

  if lib.mgxs_types.count('consistent scatter matrix') > 0 and not has_scattering_xs:
    has_scattering_xs = True
    xs_string_list.append('consistent scatter matrix')

  if lib.mgxs_types.count('consistent nu-scatter matrix') > 0 and not has_scattering_xs:
    has_scattering_xs = True
    xs_string_list.append('consistent nu-scatter matrix')

  if not has_scattering_xs:
    raise Exception("This library doesn't contain a scattering matrix.")

  # Check for cross-sections required for fission.
  has_fission = False
  if lib.mgxs_types.count('nu-fission') > 0 and lib.mgxs_types.count('chi') > 0:
    has_fission = True
    xs_string_list.append('nu-fission')
    xs_string_list.append('chi')

  # Add fission heating if it's available.
  if lib.mgxs_types.count('kappa-fission') > 0 and has_fission:
    xs_string_list.append('kappa-fission')

  # Check for diffusion coefficients.
  has_diffusion = False
  if lib.mgxs_types.count('diffusion-coefficient') > 0:
    has_diffusion = True
    xs_string_list.append('diffusion-coefficient')

  if lib.mgxs_types.count('nu-diffusion-coefficient') > 0 and not has_diffusion:
    xs_string_list.append('nu-diffusion-coefficient')

  # Condense cross-sections if requested.
  lib_to_use = lib
  if condense:
    lib_to_use = lib.get_condensed_library(coarse_group_structure)

  root_elem = ET.Element('macroscopic_cross_sections')
  root_elem.attrib['generator'] = 'openmc'
  root_elem.attrib['num_groups'] = str(lib_to_use.energy_groups.group_edges.size - 1)
  root_elem.attrib['group_bounds'] = ''
  for bounds in reversed(lib_to_use.energy_groups.group_edges):
    root_elem.attrib['group_bounds'] += str(bounds) + ' '
  root_elem.attrib['group_bounds'] = root_elem.attrib['group_bounds'].rstrip(root_elem.attrib['group_bounds'][-1])
  root_elem.attrib['energy_units'] = 'eV'
  root_elem.attrib['xs_units'] = 'cm^-1'

  for domain in lib_to_use.domains:
    material_element = ET.Element('domain')
    material_element.attrib['type'] = str(lib_to_use.domain_type)
    material_element.attrib['name'] = str(domain.name)
    material_element.attrib['id'] = str(domain.id)

    for xs_name in xs_string_list:
      xs = lib_to_use.get_mgxs(xs_name)
      reaction = xs.rxn_type

      reaction_element = ET.Element('reaction')
      reaction_element.attrib['type'] = str(reaction)

      xs_data = xs.get_xs(subdomains = 'all', nuclides= 'sum', xs_type='macro', order_groups = 'increasing', value = 'mean', squeeze = True)

      if reaction == 'scatter' or reaction == 'nu-scatter matrix' or reaction == 'scatter matrix' or reaction == 'consistent scatter matrix' or reaction == 'consistent nu-scatter matrix':
        material_element.attrib['num_legendre'] = str(xs.legendre_order)
        reaction_element.attrib['mgxs'] = ''
        reaction_element.attrib['type'] = 'scatter'
        for incoming_group in xs_data:
          for outgoing_group in incoming_group:
            for legendre_moment in outgoing_group:
              reaction_element.attrib['mgxs'] += str(legendre_moment) + ' '
        reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])
      else:
        reaction_element.attrib['mgxs'] = ''
        for xs in xs_data:
          reaction_element.attrib['mgxs'] += str(xs) + ' '
        reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])

      material_element.append(reaction_element)

    root_elem.append(material_element)

  if not os.path.exists(os.path.join('./', directory)):
    os.makedirs(os.path.join('./', directory))

  tree = ET.ElementTree(root_elem)
  tree.write(directory + '/' + file_name + '.xml', encoding='utf-8', pretty_print=True)
#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Process results from an OpenMC cross-section library for running OpenMC in multigroup mode.
#---------------------------------------------------------------------------------------------------------------------------
def process_all_macro_xs_results_from_lib_xml(mgxs_lib, file_name, directory, T_index = 0, warn = False):
  root_elem = ET.Element('macroscopic_cross_sections')
  root_elem.attrib['generator'] = 'openmc'
  root_elem.attrib['num_groups'] = str(mgxs_lib.energy_groups.group_edges.size - 1)
  root_elem.attrib['group_bounds'] = ''
  for bounds in reversed(mgxs_lib.energy_groups.group_edges):
    root_elem.attrib['group_bounds'] += str(bounds) + ' '
  root_elem.attrib['group_bounds'] = root_elem.attrib['group_bounds'].rstrip(root_elem.attrib['group_bounds'][-1])
  root_elem.attrib['energy_units'] = 'eV'
  root_elem.attrib['xs_units'] = 'cm^-1'

  id = 0
  for xs_data in mgxs_lib.xsdatas:
    material_element = ET.Element('domain')
    material_element.attrib['type'] = str('Unknown')
    material_element.attrib['name'] = str(xs_data.name)
    material_element.attrib['id'] = str(id)
    id += 1

    # Make sure we're using Legendre moments.
    if not xs_data.scatter_format == 'legendre':
      raise Exception("This library does not use Legendre moments of the scattering cross-section.")

    material_element.attrib['num_legendre'] = str(xs_data.order)

    # Check for mandatory cross-sections. Total, absorption, and a scattering matrix.
    if type(xs_data.total[T_index]) == type(None):
      raise Exception("This library doesn't contain an total cross-section.")
    if type(xs_data.absorption[T_index]) == type(None):
      raise Exception("This library doesn't contain an absorption cross-section.")
    if type(xs_data.scatter_matrix[T_index]) == type(None):
      raise Exception("This library doesn't contain an scattering matrix.")

    # Total cross-section.
    reaction_element = ET.Element('reaction')
    reaction_element.attrib['type'] = str('total')
    reaction_element.attrib['mgxs'] = ''
    for xs in xs_data.total[T_index]:
      reaction_element.attrib['mgxs'] += str(xs) + ' '
    reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])
    material_element.append(reaction_element)

    # Absorption cross-section.
    reaction_element = ET.Element('reaction')
    reaction_element.attrib['type'] = str('absorption')
    reaction_element.attrib['mgxs'] = ''
    for xs in xs_data.absorption[T_index]:
      reaction_element.attrib['mgxs'] += str(xs) + ' '
    reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])
    material_element.append(reaction_element)

    # Scattering matrix.
    reaction_element = ET.Element('reaction')
    reaction_element.attrib['type'] = str('scatter')
    reaction_element.attrib['mgxs'] = ''
    for incoming_group in xs_data.scatter_matrix[T_index]:
      for outgoing_group in incoming_group:
        for legendre_moment in outgoing_group:
          reaction_element.attrib['mgxs'] += str(legendre_moment) + ' '
    reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])
    material_element.append(reaction_element)

    # Inverse velocity.
    if type(xs_data.inverse_velocity[T_index]) != type(None):
      reaction_element = ET.Element('reaction')
      reaction_element.attrib['type'] = str('inverse-velocity')
      reaction_element.attrib['mgxs'] = ''
      for xs in xs_data.inverse_velocity[T_index]:
        reaction_element.attrib['mgxs'] += str(xs) + ' '
      reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])
      material_element.append(reaction_element)

    # Neutron production cross-section and fission spectra.
    has_fission = (type(xs_data.chi[T_index]) != type(None)) and (type(xs_data.nu_fission) != type(None))
    if (type(xs_data.chi[T_index]) == type(None)) and (type(xs_data.nu_fission) != type(None)) and warn:
      warnings.warn('A fission production matrix was provided for {} as opposed to a fission neutron spectra '
                    'and a neutron production cross-section. GNAT does not support fission production '
                    'matrices, so fission will be ignored.'.format(xs_data.name))
    if (type(xs_data.chi[T_index]) != type(None)) and type(xs_data.nu_fission[T_index]) == type(None) and warn:
      warnings.warn('No neutron production cross-sections have been provided.')

    if has_fission:
      reaction_element = ET.Element('reaction')
      reaction_element.attrib['type'] = str('nu-fission')
      reaction_element.attrib['mgxs'] = ''
      for xs in xs_data.nu_fission[T_index]:
        reaction_element.attrib['mgxs'] += str(xs) + ' '
      reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])
      material_element.append(reaction_element)

      reaction_element = ET.Element('reaction')
      reaction_element.attrib['type'] = str('chi')
      reaction_element.attrib['mgxs'] = ''
      for xs in xs_data.chi[T_index]:
        reaction_element.attrib['mgxs'] += str(xs) + ' '
      reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])
      material_element.append(reaction_element)

    root_elem.append(material_element)

  if not os.path.exists(os.path.join('./', directory)):
    os.makedirs(os.path.join('./', directory))

  tree = ET.ElementTree(root_elem)
  tree.write(directory + '/' + file_name + '.xml', encoding='utf-8', pretty_print=True)
#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Process all cross-section results from a domain-first cross-section second list.
#---------------------------------------------------------------------------------------------------------------------------
def process_all_macro_xs_results_xml(state_point, domain_list, file_name, directory,
                                     condense = False, coarse_group_structure = mgxs.EnergyGroups(np.array([0.0, 20.0e6]))):
  new_domain_list = []

  for xs_list in domain_list:
    new_xs_list = []
    for xs in xs_list:
      xs.load_from_statepoint(state_point)

      if (condense):
        new_xs_list.append(xs.get_condensed_xs(coarse_group_structure))

    new_domain_list.append(new_xs_list)

  if (condense):
    domain_list = new_domain_list

  root_elem = ET.Element('macroscopic_cross_sections')
  root_elem.attrib['generator'] = 'openmc'
  root_elem.attrib['num_groups'] = str(domain_list[0][0].energy_groups.group_edges.size - 1)
  root_elem.attrib['group_bounds'] = ''
  for bounds in reversed(domain_list[0][0].energy_groups.group_edges):
    root_elem.attrib['group_bounds'] += str(bounds) + ' '
  root_elem.attrib['group_bounds'] = root_elem.attrib['group_bounds'].rstrip(root_elem.attrib['group_bounds'][-1])
  root_elem.attrib['energy_units'] = 'eV'
  root_elem.attrib['xs_units'] = 'cm^-1'

  for xs_list in domain_list:
    material_element = ET.Element('domain')
    material_element.attrib['type'] = str(xs_list[0].domain_type)
    material_element.attrib['name'] = str(xs_list[0].domain.name)
    material_element.attrib['id'] = str(xs_list[0].domain.id)

    for xs in xs_list:
      reaction = xs.rxn_type

      reaction_element = ET.Element('reaction')
      reaction_element.attrib['type'] = str(reaction)

      xs_data = xs.get_xs(subdomains = 'all', nuclides= 'sum', xs_type='macro', order_groups = 'increasing', value = 'mean', squeeze = True)

      if reaction == 'scatter':
        material_element.attrib['num_legendre'] = str(xs.legendre_order)
        reaction_element.attrib['mgxs'] = ''
        for incoming_group in xs_data:
          for outgoing_group in incoming_group:
            for legendre_moment in outgoing_group:
              reaction_element.attrib['mgxs'] += str(legendre_moment) + ' '
        reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])
      else:
        reaction_element.attrib['mgxs'] = ''
        for xs in xs_data:
          reaction_element.attrib['mgxs'] += str(xs) + ' '
        reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])

      material_element.append(reaction_element)

    root_elem.append(material_element)

  if not os.path.exists(os.path.join('./', directory)):
    os.makedirs(os.path.join('./', directory))

  tree = ET.ElementTree(root_elem)
  tree.write(directory + '/' + file_name + '.xml', encoding='utf-8', pretty_print=True)
#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Process all depletion results from a depletion list.
#---------------------------------------------------------------------------------------------------------------------------
def process_depletion_results(state_point, depletion_list, file_name, directory):
  for dep in depletion_list:
    dep.load_from_statepoint(state_point)

  root_elem = ET.Element('depletion_chain')
  root_elem.attrib['generator'] = 'openmc'
  root_elem.attrib['num_groups'] = str(depletion_list[0].energy_groups.group_edges.size - 1)
  root_elem.attrib['group_bounds'] = ''
  for bounds in reversed(depletion_list[0].energy_groups.group_edges):
    root_elem.attrib['group_bounds'] += str(bounds) + ' '
  root_elem.attrib['group_bounds'] = root_elem.attrib['group_bounds'].rstrip(root_elem.attrib['group_bounds'][-1])
  root_elem.attrib['energy_units'] = 'eV'
  root_elem.attrib['xs_units'] = 'barns'

  nuclides = depletion_list[0].nuclides
  for nuclide in nuclides:
    nuclide_element = ET.Element('nuclide')
    nuclide_element.attrib['name'] = nuclide

    num_xs = 0
    for dep in depletion_list:
      reaction = dep.rxn_type

      reaction_element = ET.Element('reaction')
      reaction_element.attrib['type'] = str(reaction)

      xs_data = dep.get_xs(groups = 'all', subdomains = 'all', nuclides = [nuclide], xs_type='micro', order_groups = 'increasing', value = 'mean', squeeze = True)

      if (np.all(np.isclose(xs_data, 0.0, rtol=1e-16, atol=1e-16))):
        continue

      reaction_element.attrib['mgxs'] = ''
      for xs in xs_data:
        reaction_element.attrib['mgxs'] += str(xs) + ' '
      reaction_element.attrib['mgxs'] = reaction_element.attrib['mgxs'].rstrip(reaction_element.attrib['mgxs'][-1])
      nuclide_element.append(reaction_element)
      num_xs += 1

    nuclide_element.attrib['reactions'] = str(num_xs)
    root_elem.append(nuclide_element)

  if not os.path.exists(os.path.join('./', directory)):
    os.makedirs(os.path.join('./', directory))

  tree = ET.ElementTree(root_elem)
  tree.write(directory + '/' + file_name + '.xml', encoding='utf-8', pretty_print=True)

  return depletion_list
#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# The SCALE 252 neutron group cross-sections.
#---------------------------------------------------------------------------------------------------------------------------
SCALE_252 = np.array([
    0., 1.e-4, 5.e-4, 7.5e-4, 1.e-3, 1.2e-3, 1.5e-3, 2.e-3, 2.5e-3, 3.e-3,
    4.e-3, 5.e-3, 7.5e-3, 1.e-2, 2.53e-2, 3.e-2, 4.e-2, 5.e-2, 6.e-2, 7.e-2,
    8.e-2, 9.e-2, 1.e-1, 1.25e-1, 1.5e-1, 1.75e-1, 2.e-1, 2.25e-1, 2.5e-1,
    2.75e-1, 3.e-1, 3.25e-1, 3.5e-1, 3.75e-1, 4.e-1, 4.5e-1, 5.e-1, 5.5e-1,
    6.e-1, 6.25e-1, 6.5e-1, 7.e-1, 7.5e-1, 8.e-1, 8.5e-1, 9.e-1, 9.25e-1,
    9.5e-1, 9.75e-1, 1., 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09,
    1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.175, 1.2, 1.225, 1.25, 1.3, 1.35, 1.4,
    1.45, 1.5, 1.59, 1.68, 1.77, 1.86, 1.94, 2., 2.12, 2.21, 2.3, 2.38, 2.47,
    2.57, 2.67, 2.77, 2.87, 2.97, 3., 3.1, 3.2, 3.5, 3.73, 4.1, 4.7, 5., 5.4,
    6., 6.25, 6.5, 6.75, 6.875, 7., 7.15, 8.1, 9.1, 1.e+1, 1.15e+1, 1.19e+1,
    1.29e+1, 1.44e+1, 1.6e+1, 1.7e+1, 1.85e+1, 1.94e+1, 2.e+1, 2.05e+1,
    2.12e+1, 2.175e+1, 2.25e+1, 2.5e+1, 2.75e+1, 3.e+1, 3.125e+1, 3.175e+1,
    3.325e+1, 3.375e+1, 3.5e+1, 3.55e+1, 3.6e+1, 3.7e+1, 3.713e+1, 3.727e+1,
    3.763e+1, 3.8e+1, 3.91e+1, 3.96e+1, 4.1e+1, 4.24e+1, 4.4e+1, 4.52e+1,
    4.83e+1, 5.06e+1, 5.34e+1, 5.8e+1, 6.1e+1, 6.3e+1, 6.5e+1, 6.75e+1, 7.2e+1,
    7.6e+1, 8.e+1, 8.17e+1, 9.e+1, 9.7e+1, 1.012e+2, 1.05e+2, 1.08e+2, 1.13e+2,
    1.16e+2, 1.175e+2, 1.19e+2, 1.22e+2, 1.43e+2, 1.7e+2, 1.8e+2, 1.877e+2,
    1.885e+2, 1.915e+2, 1.93e+2, 2.02e+2, 2.074e+2, 2.095e+2, 2.2e+2, 2.4e+2,
    2.85e+2, 3.05e+2, 5.5e+2, 6.7e+2, 6.83e+2, 9.5e+2, 1.15e+3, 1.5e+3,
    1.55e+3, 1.8e+3, 2.2e+3, 2.25e+3, 2.5e+3, 3.e+3, 3.74e+3, 3.9e+3, 5.7e+3,
    8.03e+3, 9.5e+3, 1.3e+4, 1.7e+4, 2.e+4, 3.e+4, 4.5e+4, 5.e+4, 5.2e+4,
    6.e+4, 7.3e+4, 7.5e+4, 8.2e+4, 8.5e+4, 1.e+5, 1.283e+5, 1.49e+5, 2.e+5,
    2.7e+5, 3.3e+5, 4.e+5, 4.2e+5, 4.4e+5, 4.7e+5, 4.92e+5, 5.5e+5, 5.73e+5,
    6.e+5, 6.7e+5, 6.79e+5, 7.5e+5, 8.2e+5, 8.611e+5, 8.75e+5, 9.e+5, 9.2e+5,
    1.01e+6, 1.1e+6, 1.2e+6, 1.25e+6, 1.317e+6, 1.356e+6, 1.4e+6, 1.5e+6,
    1.85e+6, 2.354e+6, 2.479e+6, 3.e+6, 4.304e+6, 4.8e+6, 6.434e+6, 8.187e+6,
    1.e+7, 1.284e+7, 1.384e+7, 1.455e+7, 1.568e+7, 1.733e+7, 2.e+7])
#---------------------------------------------------------------------------------------------------------------------------
