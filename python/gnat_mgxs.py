import os

import numpy as np
import lxml.etree as ET

import openmc

GNAT_SUPPORTED_MACRO_MGXS = ['total', 'scatter', 'inverse-vel', 'nu-fission', 'chi', 'kappa-fission', 'diffusion', 'absorption']
GNAT_SUPPORTED_MICRO_MGXS = ['(n,2n)', '(n,3n)', '(n,4n)', '(n,gamma)', '(n,p)', '(n,a)']

# A class which initializes and stores data required to compute a collection of
# multi-group cross sections for a single homogenization domain.
class GnatOpenMCMGXS:
  def __init__(self, name, groups, domain=None, macro=GNAT_SUPPORTED_MACRO_MGXS,\
               micro=[], correction=None, legendre_order=0):
    if not isinstance(name, str):
      raise TypeError('Name must be a string')
    self._name = name

    if not isinstance(groups, openmc.mgxs.EnergyGroups):
      raise TypeError('Groups must be of type openmc.mgxs.EnergyGroups')
    self._groups = groups

    self._domain = domain
    if domain is not None:
      if isinstance(domain, openmc.Material):
        self._domain_type = 'material'
      elif isinstance(domain, openmc.Cell):
        self._domain_type = 'cell'
      elif isinstance(domain, openmc.Universe):
        self._domain_type = 'universe'
      elif isinstance(domain, openmc.RegularMesh):
        self._domain_type = 'mesh'
      else:
        raise TypeError('Invalid MGXS domain ' + str(domain))

    for mgxs in macro:
      if mgxs not in GNAT_SUPPORTED_MACRO_MGXS:
        raise Exception('Invalid macroscopic cross section type ' + str(mgxs))
    self._macro_types = macro
    for mgxs in micro:
      if mgxs not in GNAT_SUPPORTED_MICRO_MGXS:
        raise Exception('Invalid microscopic cross section type ' + str(mgxs))
    self._micro_types = micro

    if correction != None and correction != 'P0':
      raise Exception('Correction must be either None or P0')
    self._correction = correction

    if not isinstance(legendre_order, int):
      raise TypeError('legendre_order must be an integer')
    if legendre_order < 0:
      raise Exception('legendre_order must be >= 0')
    self._legendre_order = legendre_order

    self._macro_xs_objs = {}
    self._micro_xs_objs = {}

    if domain != None:
      # Prepare macroscopic cross section objects.
      for mgxs in macro:
        if mgxs == 'total':
          self._macro_xs_objs[mgxs] = openmc.mgxs.TotalXS(domain=self._domain,\
                                                          domain_type=self._domain_type,\
                                                          energy_groups=self._groups,\
                                                          name=self._name)
        if mgxs == 'inverse-vel':
          self._macro_xs_objs[mgxs] = openmc.mgxs.InverseVelocity(domain=self._domain,\
                                                                  domain_type=self._domain_type,\
                                                                  energy_groups=self._groups,\
                                                                  name=self._name)
        if mgxs == 'nu-fission':
          self._macro_xs_objs[mgxs] = openmc.mgxs.FissionXS(domain=self._domain,\
                                                            domain_type=self._domain_type,\
                                                            energy_groups=self._groups,\
                                                            nu=True,\
                                                            name=self._name)
        if mgxs == 'chi':
          self._macro_xs_objs[mgxs] = openmc.mgxs.Chi(domain=self._domain,\
                                                      domain_type=self._domain_type,\
                                                      energy_groups=self._groups,\
                                                      name=self._name)
        if mgxs == 'kappa-fission':
          self._macro_xs_objs[mgxs] = openmc.mgxs.KappaFissionXS(domain=self._domain,\
                                                                 domain_type=self._domain_type,\
                                                                 energy_groups=self._groups,\
                                                                 name=self._name)
        if mgxs == 'diffusion':
          self._macro_xs_objs[mgxs] = openmc.mgxs.DiffusionCoefficient(domain=self._domain,\
                                                                       domain_type=self._domain_type,\
                                                                       energy_groups=self._groups,\
                                                                       nu=True,\
                                                                       name=self._name)
        if mgxs == 'scatter':
          self._macro_xs_objs[mgxs] = openmc.mgxs.ScatterMatrixXS(domain=self._domain,\
                                                                  domain_type=self._domain_type,\
                                                                  energy_groups=self._groups,\
                                                                  nu=True,\
                                                                  name=self._name)
          self._macro_xs_objs[mgxs].formulation = 'consistent'
          self._macro_xs_objs[mgxs].legendre_order = legendre_order
          self._macro_xs_objs[mgxs].correction = correction

      # Prepare microscopic cross section objects.
      for mgxs in micro:
        self._micro_xs_objs[mgxs] = openmc.mgxs.ArbitraryXS(rxn_type=mgxs,\
                                                            domain=self._domain,\
                                                            domain_type=self._domain_type,\
                                                            energy_groups=self._groups,\
                                                            by_nuclide=True,\
                                                            name=self._name)
      if len(self._micro_xs_objs) > 0:
        self._nuclides = next(iter(self._micro_xs_objs.values())).nuclides
      else:
        self._nuclides = []

  # Get an OpenMC.Tallies instance containing all of the tallies necessary to compute the requested MGXS'.
  def get_tallies(self):
    tallies = openmc.Tallies()
    for mgxs in self._macro_xs_objs.values():
      tallies += mgxs.tallies.values()
    for mgxs in self._micro_xs_objs.values():
      tallies += mgxs.tallies.values()
    return tallies

  # Load the tally data from a given statepoint.
  def load_from_statepoint(self, sp):
    for mgxs in self._macro_xs_objs.values():
      mgxs.load_from_statepoint(sp)
    for mgxs in self._micro_xs_objs.values():
      mgxs.load_from_statepoint(sp)

  # Write the microscopic cross section results for this specific domain and nuclide set to a xml node.
  def to_micro_xml_node(self):
    xml_element = ET.Element('domain')
    xml_element.attrib['type'] = str(self._domain_type)
    xml_element.attrib['name'] = str(self._name)
    xml_element.attrib['id'] = str(self._domain.id)

    for nuclide in self._nuclides:
      nuclide_xml_element = ET.Element('nuclide')
      nuclide_xml_element.attrib['name'] = nuclide

      num_xs = 0
      for rxn_type, mgxs in zip(self._micro_xs_objs.keys(), self._micro_xs_objs.values()):
        reaction_xml_element = ET.Element('reaction')
        reaction_xml_element.attrib['type'] = str(rxn_type)

        xs_data = mgxs.get_xs(groups='all', subdomains='all', nuclides=[nuclide], xs_type='micro', order_groups='increasing', value='mean', squeeze=True)
        if (np.all(np.isclose(xs_data, 0.0, rtol=1e-16, atol=1e-16))):
          continue

        reaction_xml_element.attrib['mgxs'] = ''
        for xs in xs_data:
          reaction_xml_element.attrib['mgxs'] += str(xs) + ' '
        reaction_xml_element.attrib['mgxs'] = reaction_xml_element.attrib['mgxs'].rstrip(reaction_xml_element.attrib['mgxs'][-1])

        nuclide_xml_element.append(reaction_xml_element)
        num_xs += 1

      nuclide_xml_element.attrib['reactions'] = str(num_xs)
      xml_element.append(nuclide_xml_element)

    return xml_element

  # Write the microscopic cross section results to an xml file readable by Gnat.
  def to_micro_xml_file(self, file, directory):
    root_elem = ET.Element('depletion_chain')
    root_elem.attrib['generator'] = 'openmc'
    root_elem.attrib['num_groups'] = str(self._groups.group_edges.size - 1)
    root_elem.attrib['group_bounds'] = ''
    for bounds in reversed(self._groups.group_edges):
      root_elem.attrib['group_bounds'] += str(bounds) + ' '
    root_elem.attrib['group_bounds'] = root_elem.attrib['group_bounds'].rstrip(root_elem.attrib['group_bounds'][-1])
    root_elem.attrib['energy_units'] = 'eV'
    root_elem.attrib['xs_units'] = 'barns'
    root_elem.append(self.to_micro_xml_node())

    if not os.path.exists(os.path.join('./', directory)):
      os.makedirs(os.path.join('./', directory))

    tree = ET.ElementTree(root_elem)
    tree.write(directory + '/' + file + '.xml', encoding='utf-8', pretty_print=True)

  # Write the macroscopic cross section results for this specific domain to a xml node.
  def to_macro_xml_node(self):
    xml_element = ET.Element('domain')
    xml_element.attrib['type'] = str(self._domain_type)
    xml_element.attrib['name'] = str(self._name)
    xml_element.attrib['id'] = str(self._domain.id)
    for rxn_type, mgxs in zip(self._macro_xs_objs.keys(), self._macro_xs_objs.values()):
      reaction_element = ET.Element('reaction')
      reaction_element.attrib['type'] = str(rxn_type)

      xs_data = mgxs.get_xs(subdomains='all', nuclides='sum', xs_type='macro',\
                            order_groups='increasing', value='mean', squeeze=True)

      if rxn_type == 'scatter':
        xml_element.attrib['num_legendre'] = str(mgxs.legendre_order)
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

      xml_element.append(reaction_element)

    return xml_element

  # Write the macroscopic cross section results to an xml file readable by Gnat.
  def to_macro_xml_file(self, file, directory):
    root_elem = ET.Element('macroscopic_cross_sections')
    root_elem.attrib['generator'] = 'openmc'
    root_elem.attrib['num_groups'] = str(self._groups.group_edges.size - 1)
    root_elem.attrib['group_bounds'] = ''
    for bounds in reversed(self._groups.group_edges):
      root_elem.attrib['group_bounds'] += str(bounds) + ' '
    root_elem.attrib['group_bounds'] = root_elem.attrib['group_bounds'].rstrip(root_elem.attrib['group_bounds'][-1])
    root_elem.attrib['energy_units'] = 'eV'
    root_elem.attrib['xs_units'] = 'cm^-1'
    root_elem.append(self.to_macro_xml_node())

    if not os.path.exists(os.path.join('./', directory)):
      os.makedirs(os.path.join('./', directory))

    tree = ET.ElementTree(root_elem)
    tree.write(directory + '/' + file + '.xml', encoding='utf-8', pretty_print=True)

# A class which manages a collection of multi-group cross sections for a list of
# homogenization domains.
class OpenMCMGXSCollection:
  def __init__(self, names, groups, domains, macro=GNAT_SUPPORTED_MACRO_MGXS,\
               micro=[], correction=None, legendre_order=0):
    if not hasattr(names, '__iter__'):
      raise TypeError('names must be iterable')
    if not hasattr(domains, '__iter__'):
      raise TypeError('domains must be iterable')
    if len(names) != len(domains):
      raise Exception('names and domains must be the same size')

    if not isinstance(groups, openmc.mgxs.EnergyGroups):
      raise TypeError('Groups must be of type openmc.mgxs.EnergyGroups')
    self._groups = groups

    self._mgxs_objs = []
    for name, domain in zip(names, domains):
      self._mgxs_objs.append(GnatOpenMCMGXS(name, groups, domain, macro,\
                                            micro, correction, legendre_order))

  # Get an OpenMC.Tallies instance containing all of the tallies necessary to compute the requested
  # cross sections over all homogenization domains.
  def get_tallies(self):
    tallies = openmc.Tallies()
    for mgxs in self._mgxs_objs:
      tallies += mgxs.get_tallies()
    return tallies

  # Load the tally data from a given statepoint.
  def load_from_statepoint(self, sp):
    for mgxs in self._mgxs_objs:
      mgxs.load_from_statepoint(sp)

  # Write the microscopic cross section results to an xml file readable by Gnat.
  def to_micro_xml_file(self, file, directory):
    root_elem = ET.Element('depletion_chain')
    root_elem.attrib['generator'] = 'openmc'
    root_elem.attrib['num_groups'] = str(self._groups.group_edges.size - 1)
    root_elem.attrib['group_bounds'] = ''
    for bounds in reversed(self._groups.group_edges):
      root_elem.attrib['group_bounds'] += str(bounds) + ' '
    root_elem.attrib['group_bounds'] = root_elem.attrib['group_bounds'].rstrip(root_elem.attrib['group_bounds'][-1])
    root_elem.attrib['energy_units'] = 'eV'
    root_elem.attrib['xs_units'] = 'barns'

    for mgxs in self._mgxs_objs:
      root_elem.append(mgxs.to_micro_xml_node())

    if not os.path.exists(os.path.join('./', directory)):
      os.makedirs(os.path.join('./', directory))

    tree = ET.ElementTree(root_elem)
    tree.write(directory + '/' + file + '.xml', encoding='utf-8', pretty_print=True)

  # Write the macroscopic cross section results to an xml file readable by Gnat.
  def to_macro_xml_file(self, file, directory):
    root_elem = ET.Element('macroscopic_cross_sections')
    root_elem.attrib['generator'] = 'openmc'
    root_elem.attrib['num_groups'] = str(self._groups.group_edges.size - 1)
    root_elem.attrib['group_bounds'] = ''
    for bounds in reversed(self._groups.group_edges):
      root_elem.attrib['group_bounds'] += str(bounds) + ' '
    root_elem.attrib['group_bounds'] = root_elem.attrib['group_bounds'].rstrip(root_elem.attrib['group_bounds'][-1])
    root_elem.attrib['energy_units'] = 'eV'
    root_elem.attrib['xs_units'] = 'cm^-1'

    for mgxs in self._mgxs_objs:
      root_elem.append(mgxs.to_macro_xml_node())

    if not os.path.exists(os.path.join('./', directory)):
      os.makedirs(os.path.join('./', directory))

    tree = ET.ElementTree(root_elem)
    tree.write(directory + '/' + file + '.xml', encoding='utf-8', pretty_print=True)
