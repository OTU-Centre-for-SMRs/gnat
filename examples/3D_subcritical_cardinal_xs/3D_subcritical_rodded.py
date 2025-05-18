import openmc

# Subcritical geometric properties.
## Fuel and sheath.
r_fuel = 3.2766 / 2.0
r_cladding_outer = 3.683 / 2.0

## Graphite block.
r_hole = 3.96748 / 2.0
l_w_block = 10.16

## Graphite plug.
r_plug = 3.81 / 2.0

## Control rod.
r_control = r_plug

# Declare the materials.
## Air
n2 = openmc.Material(name = "n2")
n2.add_element('N', 2.0)

o2 = openmc.Material(name = "o2")
o2.add_element('O', 2.0)

ar = openmc.Material(name = "ar")
ar.add_element('Ar', 1.0)

co2 = openmc.Material(name = "co2")
co2.add_element('C', 1.0)
co2.add_element('O', 2.0)

comp = [n2, o2, ar, co2]
comp_frac = [78.084 / 100.0, 20.946 / 100.0, 0.9340 / 100.0]

total = 0.0
for fraction in comp_frac:
  total = total + fraction

comp_frac.append(1.0 - total)

air = openmc.Material.mix_materials(comp, comp_frac, 'wo')
air.set_density('g/cc', 0.0012)
air.name = 'Air'

## Fuel, metallic uranium.
fuel = openmc.Material(name="Fuel")
fuel.add_element('U', 1.0, enrichment=0.72)
fuel.set_density('g/cc', 19.1)

## Sheath, metallic aluminum.
cladding = openmc.Material(name="Cladding")
cladding.add_element('Al', 1.0)
cladding.set_density('g/cc', 2.699)

## Moderator, graphite.
graphite_1 = openmc.Material(name="Graphite 1 (Internal)")
graphite_1.add_element('C', 1.0)
graphite_1.set_density('g/cc', 1.91)
graphite_1.add_s_alpha_beta('c_Graphite')

graphite_2 = openmc.Material(name="Graphite 2 (Reflector)")
graphite_2.add_element('C', 1.0)
graphite_2.set_density('g/cc', 1.91)
graphite_2.add_s_alpha_beta('c_Graphite')

## Control rods, stainless-steel.
steel = openmc.Material(name='Control')
steel.set_density('g/cm3', 8.00)
steel.add_element('C', 0.08, percent_type='wo')
steel.add_element('Si', 1.00, percent_type='wo')
steel.add_element('Mn', 2.00, percent_type='wo')
steel.add_element('P', 0.045, percent_type='wo')
steel.add_element('S', 0.030, percent_type='wo')
steel.add_element('Cr', 20.0, percent_type='wo')
steel.add_element('Ni', 11.0, percent_type='wo')
steel.add_element('Fe', 65.845, percent_type='wo')

# Declare geometry.
fuel_outer_radius = openmc.ZCylinder(r = r_fuel)
clad_outer_radius = openmc.ZCylinder(r = r_cladding_outer)

control_rod_outer_radius = openmc.ZCylinder(r = r_control)

hole_outer_radius = openmc.ZCylinder(r = r_hole)
single_block_box = openmc.model.RectangularPrism(width=l_w_block, height=l_w_block)

plug_outer_radius = openmc.ZCylinder(r = r_plug)

## A fuel element.
fuel_rod = openmc.Cell(name = 'Fuel Element')
fuel_rod.region = -fuel_outer_radius
fuel_rod.fill = fuel

## The fuel cladding.
fuel_cladding = openmc.Cell(name = 'Fuel Cladding')
fuel_cladding.region = +fuel_outer_radius & -clad_outer_radius
fuel_cladding.fill = cladding

## The fuel cladding -> graphite block air gap.
fuel_clad_block_gap = openmc.Cell(name = 'Fuel Air Gap')
fuel_clad_block_gap.region = +clad_outer_radius & -hole_outer_radius
fuel_clad_block_gap.fill = air

## The graphite plug.
graphite_plug_1 = openmc.Cell(name = 'Graphite Plug 1 (Interior)')
graphite_plug_1.region = -plug_outer_radius
graphite_plug_1.fill = graphite_1

graphite_plug_2 = openmc.Cell(name = 'Graphite Plug 2 (Reflector)')
graphite_plug_2.region = -plug_outer_radius
graphite_plug_2.fill = graphite_2

## The graphite plug -> graphite block air gap.
graphite_plug_block_gap_1 = openmc.Cell(name = 'Plug Air Gap 1')
graphite_plug_block_gap_1.region = +plug_outer_radius & -hole_outer_radius
graphite_plug_block_gap_1.fill = air

graphite_plug_block_gap_2 = openmc.Cell(name = 'Plug Air Gap 2')
graphite_plug_block_gap_2.region = +plug_outer_radius & -hole_outer_radius
graphite_plug_block_gap_2.fill = air

## A control rod.
control_rod = openmc.Cell(name = 'Control Rod')
control_rod.region = -control_rod_outer_radius
control_rod.fill = steel

## The control rod -> graphite block air gap.
control_block_gap = openmc.Cell(name = 'Control Rod Gap')
control_block_gap.region = -hole_outer_radius & +control_rod_outer_radius
control_block_gap.fill = air

## An empty air segment.
empty_air = openmc.Cell(name = 'Empty Block Hole')
empty_air.region = -hole_outer_radius
empty_air.fill = air

## A graphite block.
graphite_block_1 = openmc.Cell(name = 'Graphite Block 1')
graphite_block_1.region = -single_block_box & +hole_outer_radius
graphite_block_1.fill = graphite_1

graphite_block_2 = openmc.Cell(name = 'Graphite Block 2')
graphite_block_2.region = -single_block_box & +hole_outer_radius
graphite_block_2.fill = graphite_1

graphite_block_3 = openmc.Cell(name = 'Graphite Block 3')
graphite_block_3.region = -single_block_box & +hole_outer_radius
graphite_block_3.fill = graphite_1

graphite_block_4 = openmc.Cell(name = 'Removed Block')
graphite_block_4.region = -single_block_box & +hole_outer_radius
graphite_block_4.fill = air

graphite_block_5 = openmc.Cell(name = 'Graphite Block 5')
graphite_block_5.region = -single_block_box & +hole_outer_radius
graphite_block_5.fill = graphite_2

## The individual lattice elements.
f = openmc.Universe(cells=[fuel_rod, fuel_cladding, fuel_clad_block_gap, graphite_block_1])
i = openmc.Universe(cells=[graphite_plug_1, graphite_plug_block_gap_1, graphite_block_2])
c = openmc.Universe(cells=[control_rod, control_block_gap, graphite_block_3])
e = openmc.Universe(cells=[empty_air, graphite_block_4])
r = openmc.Universe(cells=[graphite_plug_2, graphite_plug_block_gap_2, graphite_block_5])

# The entire subcritical assembly
assembly_map_2D = [
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r], #1
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r], #2
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r], #3
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r], #4
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r], #5
  [r,  r,  r,  r,  r,  f,  i,  f,  i,  f,  i,  f,  i,  f,  r,  r,  r,  r,  r], #6  |
  [r,  r,  r,  r,  r,  i,  f,  i,  f,  i,  f,  i,  f,  i,  r,  r,  r,  r,  r], #7  |
  [r,  r,  r,  r,  r,  f,  i,  f,  i,  f,  i,  f,  i,  f,  r,  r,  r,  r,  r], #8  |
  [r,  r,  r,  r,  r,  i,  f,  i,  i,  c,  i,  i,  f,  i,  r,  r,  r,  r,  r], #9  |
  [r,  r,  r,  r,  r,  f,  i,  f,  c,  e,  c,  f,  i,  f,  r,  r,  r,  r,  r], #10 |
  [r,  r,  r,  r,  r,  i,  f,  i,  i,  c,  i,  i,  f,  i,  r,  r,  r,  r,  r], #11 |
  [r,  r,  r,  r,  r,  f,  i,  f,  i,  f,  i,  f,  i,  f,  r,  r,  r,  r,  r], #12 |
  [r,  r,  r,  r,  r,  i,  f,  i,  f,  i,  f,  i,  f,  i,  r,  r,  r,  r,  r], #13 |
  [r,  r,  r,  r,  r,  f,  i,  f,  i,  f,  i,  f,  i,  f,  r,  r,  r,  r,  r], #14 |
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r], #15
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r], #16
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r], #17
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r], #18
  [r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r,  r]  #19
#  1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19
#                     ------------------------------------
]

NUM_LAYERS = 1
assembly_cells = [ assembly_map_2D for i in range(NUM_LAYERS) ]

assembly = openmc.RectLattice(name='Subcritical Assembly')
assembly.pitch = (l_w_block, l_w_block, (152.4 / 2.0) / float(NUM_LAYERS))
assembly.lower_left = (-19.0 * l_w_block / 2.0, -19 * l_w_block / 2.0, 0.0)
assembly.universes = assembly_cells

left   = openmc.XPlane(x0 = 0.0, boundary_type = 'reflective')
right  = openmc.XPlane(x0 = 19.0 * l_w_block / 2.0, boundary_type = 'vacuum')
back   = openmc.YPlane(y0 = 0.0, boundary_type = 'reflective')
front  = openmc.YPlane(y0 = 19.0 * l_w_block / 2.0, boundary_type = 'vacuum')
bottom = openmc.ZPlane(z0 = 0.0, boundary_type = 'reflective')
top    = openmc.ZPlane(z0 = 152.4 / 2.0, boundary_type = 'vacuum')

assembly_region = +left & -right & +back & -front & +bottom & -top
full_assembly_cell = openmc.Cell(name='Assembly Cell', fill = assembly, region = assembly_region)
assembly_universe = openmc.Universe(cells=[full_assembly_cell])

# The D-T neutron generator.
#point = openmc.stats.Point((0.0, 0.0, 0.0))
#energy = openmc.stats.Discrete([14.08e6], [1.0])

# A uniform starter source to reduce variance.
uniform = openmc.stats.Box((0.0, 0.0, 0.0), (19.0 * l_w_block / 2.0, 19.0 * l_w_block / 2.0, 152.4 / 2.0))

# The model
subcritical_model = openmc.Model()
subcritical_model.materials = openmc.Materials([fuel, cladding, graphite_1, graphite_2, steel, air])
subcritical_model.geometry = openmc.Geometry(assembly_universe)

# The simulation settings.
#subcritical_model.settings.source = [openmc.IndependentSource(space = point, energy = energy)]
subcritical_model.settings.source = [openmc.IndependentSource(space = uniform)]
subcritical_model.settings.batches = 100
subcritical_model.settings.inactive = 10
subcritical_model.settings.particles = 1000

ADD_TALLIES = False
RUN_OPENMC = False

if ADD_TALLIES:
  groups =  openmc.mgxs.EnergyGroups([0.0, 2.38, 111000.0, 14900000.0])
  scatter_matrix =  openmc.mgxs.ScatterMatrixXS(domain = full_assembly_cell,
                                                domain_type = 'cell',
                                                energy_groups = groups,
                                                nu = True)
  scatter_matrix.formulation = 'simple'
  scatter_matrix.correction = 'P0'
  scatter_matrix.scatter_format = 'legendre'
  scatter_matrix.legendre_order = 0
  subcritical_model.tallies += scatter_matrix.tallies.values()

  total_xs = openmc.mgxs.TotalXS(domain = full_assembly_cell,
                                 domain_type = 'cell',
                                 energy_groups = groups)
  total_xs.estimator = 'analog'
  subcritical_model.tallies += total_xs.tallies.values()

  chi = openmc.mgxs.Chi(domain = full_assembly_cell,
                        domain_type = 'cell',
                        energy_groups = groups)
  chi.estimator = 'analog'
  subcritical_model.tallies += chi.tallies.values()

  nu_fission_xs = openmc.mgxs.FissionXS(domain = full_assembly_cell,
                                        domain_type = 'cell',
                                        energy_groups = groups,
                                        nu = True)
  nu_fission_xs.estimator = 'analog'
  subcritical_model.tallies += nu_fission_xs.tallies.values()

  if RUN_OPENMC:
    subcritical_model.export_to_model_xml()
    statepoint_filename = subcritical_model.run(threads = 14)
    sp = openmc.StatePoint(statepoint_filename)
    scatter_matrix.load_from_statepoint(sp)
    total_xs.load_from_statepoint(sp)
    chi.load_from_statepoint(sp)
    nu_fission_xs.load_from_statepoint(sp)

    sp.close()

    scatter_matrix.print_xs()
    total_xs.print_xs()
    chi.print_xs()
    nu_fission_xs.print_xs()
  else:
    sp = openmc.StatePoint('./statepoint.100.h5')
    scatter_matrix.load_from_statepoint(sp)

    sp.close()

    scatter_matrix.print_xs()

subcritical_model.export_to_model_xml()
