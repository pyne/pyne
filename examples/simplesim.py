import numpy as np

from pyne.simplesim import cards, definition, inputfile
import pyne.simplesim.nestedgeom as ng

class InfLattice(object):
    """Creates an MCNP input file named `ex_simplesim_inflattice` for a 2-cell
    (fuel and moderator) infinite lattice reactor. The user executes this code
    by instantiating an object of this class and calling :py:meth:`write`::

        inflat = InfLattice()
        inflat.write()

    as is done below. The user can then manipulate the input, and observe the
    change in the output::

        inflat.pin.radius = 0.45
        inflat.write()
        inflat.sim.source['criticality'].keff_guess = 1.5
        inflat.write()
    
    """
    def __init__(self):

        ## Define the system: materials, surfaces, regions, cells.
        self.sys = definition.SystemDefinition(verbose=False)

        ## Materials.
        # Must provide a name as a keyword argument for material cards. See the
        # documentation for :py:mod:`pyne.material` for more information.
        self.uo2 = cards.Material(name='UO2')
        self.uo2.from_atom_frac({'U235': 0.05,
                            'U238': 0.95,
                            'O16' : 2.00})
        self.h2o = cards.Material(name='H2O')
        self.h2o.from_atom_frac({'H1' : 2.0,
                            'O16': 1.0})

        ## Surfaces.
        # There are two surfaces: one for the pin and one for the unit cell
        # boundary.
        radius = 0.40
        # This creates an axis-aligned and axis-centered cylinder along the z
        # axis, with radius 0.40 cm.
        self.pin = cards.AxisCylinder('pin', 'Z', radius)
        # The Parallelepiped is a macrobody. The surface is reflecting,
        # creating an infinte geometry. The surface is infinite in the z
        # direction.
        pitch = 1.2
        self.cellbound = cards.Parallelepiped('bound',
                -pitch / 2, pitch / 2, -pitch / 2, pitch / 2, 0, 0,
                reflecting=True)

        ## Cells.
        # We combine the materials and surfaces above into cells. We use MCNP
        # cells in order to specify particle importances and volumes directly
        # on the cell card. We could alternatively use the
        # :py:class:`Importance` and :py:class:`Volume` cards.

        # fuel cell.
        # The fuel is the region of space inside the pin, pin.neg. 
        self.fuelregion = self.pin.neg
        # The neutron importance is 1, and the user-provided volume is 1 cm^3.
        self.fuel = cards.CellMCNP('fuel', self.fuelregion, self.uo2,
                11.0, 'g/cm^3',
                importance=('neutron', 1),
                volume=1)

        # coolant cell.
        # The region is between the pin and the unit cell boundary.
        self.coolantregion = self.pin.pos | self.cellbound.neg
        self.coolant = cards.CellMCNP('coolant', self.coolantregion, self.h2o,
                1.0, 'g/cm^3',
                importance=('neutron', 1),
                volume=1)

        # graveyard cell: where particles go to die.
        # The region is everything beyond the unit cell boundary.
        self.graveyardregion = self.cellbound.pos
        # This is a void cell, meaning it does not have a material.
        self.graveyard = cards.CellMCNP('graveyard', self.graveyardregion,
                importance=('neutron', 0))

        # We add the cells to the system. The order we add them is the order
        # they are printed in the input file.
        self.sys.add_cell(self.fuel)
        # We can add multiple cells at once.
        self.sys.add_cell(self.coolant, self.graveyard)
        
        ## Define the simulation: sources, tallies, misc. Don't clutter the
        # command window.
        self.sim = definition.MCNPSimulation(self.sys, verbose=False)

        # Specify a thermal scattering law for the H2O material. This is a
        # unique card per material.
        self.sim.add_misc(cards.ScatteringLaw('H2O', {'H1': 'lwtr'}))
        # Add a criticality source, use default values. This is a unique card,
        # so we do not provide a card name.
        self.sim.add_source(cards.Criticality())
        # Add points at which to start neutrons; use default point (0, 0, 0).
        self.sim.add_source(cards.CriticalityPoints())
        # Tally neutron flux in both the fuel and coolant cells.
        self.sim.add_tally(cards.CellFlux('flux', 'neutron', 
                ['fuel', 'coolant']))
        # The energy grid on which to tally neutrons, applied to all tallies.
        self.sim.add_misc(cards.EnergyGrid('egrid0', None,
                10**np.arange(-9.9, 1.1, 0.1)))

    def write(self):
        """Writes the input to 'ex_simplesim_inflattice'."""

        # Create input file, specifying the title of the input.
        self.inp = inputfile.MCNPInput(self.sim, title="Infinite lattice.")
        self.inp.write('ex_simplesim_inflattice')


# Create all relevant objects for the infinite lattice example.
inflat = InfLattice()
# Write to a file.
inflat.write()











"""

# super brief
rxr = simplesim.definition.SystemDefinition()
pinsurf = cards.AxisCylinder('fuelpin', 'x', 0.40)
rxr.add_cell(cards.CellMCNP('fuel', pinsurf.neg,
        material.from_atom_frac({'U235': 0.05, 'U238': 0.95, 'O16': 2.0},
            name='UO2'), 
        neutron_imp=1))
pitch = 1.2
boundsurf = cards.Parallelepiped('bound',
    -pitch / 2, pitch / 2, -pitch / 2, pitch / 2,
    0, 0, reflecting=True),
rxr.add_cell(cards.CellMCNP('coolant', pinsurf.pos & boundsurf.neg,
        material.from_atom_frac({'H1': 2.0, 'O16': 1.0}, name='H2O'),
        neutron_imp=1))
rxr.add_cell(cards.CellVoidMCNP('graveyard', boundsurf.pos, neutron_imp=0))

"""
"""
# Geometry and materials

# Materials.
uo2 = material.from_atom_frac({'U235': 0.05,
                    'U238': 0.95,
                    'O16' : 2.0}, name='UO2')

h2o = material.from_atom_frac({'H1' : 2.0,
                    'O16': 1.0}, name='H2O')


# Surfaces.
radius = 0.40 # cm
pin = cards.AxisCylinder('fuelpin', 'X', radius)

pitch = 1.2 # cm
cellbound = cards.Parallelepiped('bound',
        -pitch / 2, pitch / 2, -pitch / 2, pitch / 2, 0, 0,
        reflecting=True)

# Cells.
fuel = cards.CellMCNP('fuel', pin.neg, uo2, 11.0, 'g/cm^3',
        neutron_imp=1)
coolant = cards.CellMCNP('coolant', pin.pos & cellbound.neg, h2o,
        1.0, 'g/cm^3',
        neutron_imp=1)
graveyard = cards.CellVoidMCNP('graveyard', cellbound.pos, neutron_imp=0)

# Create system definition from the cards above.
rxr = definition.SystemDefinition()
rxr.add_cell(fuel)
rxr.add_cell(coolant)
rxr.add_cell(graveyard)
# The system definition is complete.

# Simulation definition.
sim = definition.MCNPSimulation(rxr)

sim.add_source(cards.Criticality())
sim.add_source(cards.CriticalityPoints())

fueltally = cards.CellFlux('fuel', 'neutron', fuel)
coolanttally = cards.CellFlux('coolant', 'neutron', coolant)
egrid = cards.EnergyGrid('grid0', None, 10**np.arange(-9.9, 1.1, .1))

sim.add_tally(fueltally)
sim.add_tally(coolanttally)
    
sim.add_misc(egrid)

inp = inputfile.MCNPInput("input1", sim)

print fuel.mcnp(sim)
#rxr.save('test')

"""
"""
# 
opts = definition.SimulationDefinition(rxr)

opts.add_card(cards.CriticalitySource(1000, 1, 1, 1))
opts.add_card(cards.CriticalityPoints([[0, 0, 0]]))

rxr.save("system1")

opts.save("options1")

enrichments = [0.01 0.02 0.03]

for this_enrich in enrichments:
    rxr.material['UO2'] = material(blah blah)
    inp = MCNPInput("input1", rxr, opts)
    inp.write()

    inp.system.material['UO2'] = material(blah blah)
    inp.options.card[crit_source_name].cycles = 1000
    inp.write()
    


"""
"""
# Create cards.
channel = cards.AxisCylinder("channel", 'X', 2.54)
leftbound = cards.AxisPlane("leftbound", 'X', -500.0)
rightbound = cards.AxisPlane("rightbound", 'X', 500.0)
polycyl = cards.AxisCylinder("polycyl", 'X', 17.54)
tungstencyl = cards.AxisCylinder("tungstencyl", 'X', 27.54)
coppercyl = cards.AxisCylinder("coppercyl", 'X', 28.04)
shieldleft = cards.AxisPlane("shieldleft", 'X', -25.0)
shieldright = cards.AxisPlane("shieldright", 'X', 25.0)
aperturecyl = cards.AxisCylinder("aperturecyl", 'Z', 0.25)
half = cards.AxisPlane("half", 'Z', 0.0)
gravecyl = cards.AxisCylinder("gravecyl", 'X', 33.04)

# Make regions.
pipemid = leftbound.pos & rightbound.neg # & channel.neg
polyshield = (shieldleft.pos & shieldright.neg & 
        channel.pos & polycyl.neg & 
        aperturecyl.pos )
tungstenshield = (shieldleft.pos & shieldright.neg &
        polycyl.pos & tungstencyl.neg &
        aperturecyl.pos)
coppershield = (shieldleft.pos & shieldright.neg &
        tungstencyl.pos & coppercyl.neg &
        aperturecyl.pos)
aperture = aperturecyl.neg & channel.pos & coppercyl.neg
usefulvoid = ((channel.pos & gravecyl.neg & 
              leftbound.pos & shieldleft.neg)
             |
             (channel.pos & gravecyl.neg &
              shieldright.pos & rightbound.neg)
             |
             (tungstencyl.pos & gravecyl.neg &
              shieldleft.pos & shieldright.neg))
uselessvoid = (leftbound.neg | rightbound.pos |
              (leftbound.pos & leftbound.neg &
               gravecyl.pos))

print pipemid.comment()
print polyshield.comment()
"""

"""
rxr.add_lattice

rxr.add_universe
cards.CellbyUniverse

perhaps add all surfaces and materials at once.





"""















# TODO show a bunch of ways to do a single simulation.
