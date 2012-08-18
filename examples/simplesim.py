import json

import numpy as np

from pyne.simplesim import definition, cards, inputfile
from pyne import material

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
