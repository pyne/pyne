
# Geometry and materials

# Materials.
uo2 = material.Material(name='UOX')
uo2.from_atom_frac({'U235': 0.05,
                    'U238': 0.95,
                    'O16' : 2.0})

h2o = material.Material(name='H2O')
h2o.from_atom_frac({'H1' : 2.0,
                    'O16': 1.0})


# Surfaces.
radius = 0.40 # cm
pin = cards.AxisCylinder('fuelpin', 'X', radius)

pitch = 1.2 # cm
cellbound = cards.Parallelepiped('bound',
        -pitch / 2, pitch / 2, -pitch / 2, pitch / 2, 0, 0,
        reflecting=True)


# Cells.
fuel = cards.CellSimpleMCNP('fuel', ['pin'], [], 'UOX',
        neutron_imp=1)
coolant = cards.CellSimpleMCNP('coolant', ['bound'], ['pin'], 'H2O',
        neutron_imp=1)

# Create system definition from the cards above.
rxr = simulation.definition.SystemDefinition()
rxr.add_material(uo2)
rxr.add_material(h2o)
rxr.add_surface(pin)
rxr.add_surface(cellbound)
rxr.add_cell(fuel)
rxr.add_cell(coolant)

# 
opts = simulationdef.SimulationDefinition(rxr)

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
leftbound = cards.Plane("leftbound", 'X', -500.0)
rightbound = cards.Plane("rightbound", 'X', 500.0)
polycyl = cards.AxisCylinder("polycyl", 'X', 17.54)
tungstencyl = cards.AxisCylinder("tungstencyl", 'X', 27.54)
coppercyl = cards.AxisCylinder("coppercyl", 'X', 28.04)
shieldleft = cards.Plane("shieldleft", 'X', -25.0)
shieldright = cards.Plane("shieldright", 'X', 25.0)
aperturecyl = cards.AxisCylinder("aperturecyl", 'Z', 0.25)
half = cards.Plane("half", 'Z', 0.0)
gravecyl = cards.AxisCylinder("gravecyl", 'X', 33.04)

# Make regions.
pipemid = leftbound.pos & rightbound.neg & channel.neg
polyshield = shieldleft.pos & shieldright.neg & 
        channel.pos & polycyl.neg & 
        aperturecyl.pos 
tungstenshield = shieldleft.pos & shieldright.neg &
        polycyl.pos & tungstencyl.neg &
        aperturecyl.pos
coppershield = shieldleft.pos & shieldright.neg &
        tungstencyl.pos & coppercyl.neg &
        aperturecyl.pos
aperture = aperturecyl.neg & channel.pos & coppercyl.neg
usefulvoid = (channel.pos & gravecyl.neg & 
              leftbound.pos & shieldleft.neg)
             |
             (channel.pos & gravecyl.neg &
              shieldright.pos & rightbound.neg)
             |
             (tungstencyl.pos & gravecyl.neg &
              shieldleft.pos & shieldright.neg)
uselessvoid = leftbound.neg | rightbound.pos |
              (leftbound.pos & leftbound.neg &
               gravecyl.pos)


rxr.add_lattice

rxr.add_universe
cards.CellbyUniverse

























