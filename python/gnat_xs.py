import os

import numpy as np
import lxml.etree as ET

import openmc
import openmc.mgxs as mgxs

def get_depletion_mgxs_list(domain=None, domain_type=None, energy_groups=None, name='', num_polar=1, num_azimuthal=1):
  depletion_list = []
  for reaction_type in ('(n,2n)', '(n,3n)', '(n,4n)', '(n,gamma)', '(n,p)', '(n,a)'):
    depletion_list.append(mgxs.ArbitraryXS(reaction_type, domain, domain_type, energy_groups, True, name, num_polar, num_azimuthal))

  return depletion_list

def get_mgxs_list(domain=None, domain_type=None, energy_groups=None, legendre_order=0, correction=None, name='', num_polar=1, num_azimuthal=1, nu=False):
  depletion_list = []

  depletion_list.append(mgxs.TotalXS(domain, domain_type, energy_groups, False, name, num_polar, num_azimuthal))

  depletion_list.append(mgxs.AbsorptionXS(domain, domain_type, energy_groups, False, name, num_polar, num_azimuthal))

  scatter = mgxs.ScatterMatrixXS(domain, domain_type, energy_groups, False, name, num_polar, num_azimuthal, nu)
  scatter.formulation = 'consistent'
  scatter.scatter_format = 'legendre'
  scatter.correction = correction
  scatter.legendre_order = legendre_order
  depletion_list.append(scatter)

  depletion_list.append(mgxs.InverseVelocity(domain, domain_type, energy_groups, False, name, num_polar, num_azimuthal))

  return depletion_list

def apply_tallies(tallies, xs_list):
  for xs in xs_list:
    tallies += xs.tallies.values()

  return tallies

def process_all_macro_xs_results_xml(state_point, domain_list, file_name, directory):
  for xs_list in domain_list:
    for xs in xs_list:
      xs.load_from_statepoint(state_point)

  root_elem = ET.Element('macroscopic_cross_sections')
  root_elem.attrib['generator'] = 'openmc'
  root_elem.attrib['num_groups'] = str(xs_list[0].energy_groups.group_edges.size - 1)
  root_elem.attrib['group_bounds'] = ''
  for bounds in reversed(xs_list[0].energy_groups.group_edges):
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
