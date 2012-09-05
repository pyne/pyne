"""This example creates a simple MCNP input file for a infinite lattice of
enriched UO2 fuel pins surrounded by water.

"""
import numpy as np

from pyne import mcnp

# This first (currently, only) example creates an infinite lattice of an
# enriched UO2 pin surrounded by a water moderaor.

# Create the mcnp.Inp object. A single 'modification' is provided.
inp1 = mcnp.Inp('infinitelattice', 
    'Example of PyNE mcnp.Inp class for infinite lattice',
    ("This input is generated using PyNE's mcnp.Inp class, a wrapper for "
    "input files. This particular input is of an infinite lattice of "
    "enriched UO2 fuel pins surrounded by water."), 
    'Chris Dembia', 
    [('Aug 2012', 'creation')])

# Add surfaces.
# A pin, infinite in the z direction, centered on the z axis, and with a radius
# of 0.40 cm.
radius = 0.40 # cm
inp1.add_surface_cylinder('pin', 'z', radius)
pitch = 1.2 # cm
# An RPP, 1.2 x 1.2 x infinite, centered at the origin. Reflecting to impose
# the infinite geometry.
inp1.add_surface_rectangularparallelepiped('cellbound',
        -pitch / 2, pitch / 2, -pitch / 2, pitch / 2, 0, 0,
        reflecting=True)

# Add materials.
enrichment = 0.05
# Ideally the user will be able to pass a pyne.Material object.
# The user provides a name as well as a comment. The comment must be a tuple
# so that the comment is split across multiple lines.
# Densities are provided as number ratios.
# Temperature is 600 K.
uox_temp = 600
inp1.add_material('UOX', ('5% enriched UO2',),
        ['8016', '92235', '92238'], 'atoms/b/cm',
        [2, enrichment, 1 - enrichment], uox_temp)
# Water is at 300 K.
water_temp = 300
inp1.add_material('H2O', ('Water',), 
        ['1001', '8016'], 'atoms/b/cm',
        [2, 1],
        water_temp)
# We specify a scattering law for hydrogen bound in water.
inp1.add_scattering_law('H2O', ['lwtr'])

# Add cells.
pin_vol = 3.14159 * radius**2
uox_density = 11
water_density = 1
# Fuel pin.
inp1.add_cell('UOX pin', # name.
              'UOX', # material, given above.
              uox_density,
              'g/cm^3', # units of density.
              ["pin"], # surfaces for negative sense.
              [], # surfaces for positive sense.
              1, # neutron importance.
              temp=uox_temp, # temperature (K).
              vol=pin_vol) # volume(cm^3).
inp1.add_cell('moderator',
              'H2O',
              water_density,
              'g/cm^3',
              ['cellbound'],
              ['pin'],
              1,
              water_temp,
              vol=pitch**2 - pin_vol)
inp1.add_cell_void('Problem boundary', # name.
                   [], # negative sense.
                   ['cellbound'],# positive sense.
                   0) # neutron importance 

# Source.
# Use the default arguments.
inp1.add_criticality_source()
# Use the default argument.
inp1.add_criticality_source_points()
# Flux tally (F4) in the fuel; for neutrons.
inp1.add_tally_cellflux('fuel spectrum', 'N',
        ['UOX pin'])
# Flux tally (F4) in the moderator; for neutrons.
inp1.add_tally_cellflux('moderator spectrum', 'N',
        ['moderator'])
# Energy groups for all tallies (first argument is 0).
# We use useful array operations from numpy to define the group structure to
# have equal-lethargy bins.
group_def = 10**np.arange(-9.9, 1.1, .1)
inp1.add_tally_energy(0, list(group_def))
# Tally multiplier cards for obtaining reaction rates.
# MT 1: total, MT 2: elastic scattering, MT 3: fission, MT 4: capture
inp1.add_tally_multiplier('fuel spectrum', [(1,),
        (-1, 'UOX', [1, 2, 18, 102])])
inp1.add_tally_multiplier('moderator spectrum', [(1,),
        (-1, 'H2O', [1, 2, 18, 102])])
inp1.add_printdump()

# Write the input file!
inp1.write()
