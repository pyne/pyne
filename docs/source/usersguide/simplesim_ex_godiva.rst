.. _simplesim_ex_godiva:

==============
Godiva example
==============

Not created yet because the sphere surface has not been implemented!


.. There are some errors with the following:
   # card container classes
   system = SystemDefinition()
   sim = SimulationDefinition(system)
   # uranium material
   u = Material({'U235': 1}, name='U')
   # surface with radius 1.0 cm
   sph = Sphere(1.0)
   # region of space on negative sense side of sphere.
   fuelreg = sph.neg
   # region of space on outside of sphere.
   gravereg = sph.pos
   # cell for fuel, with density
   fuel = Cell('fuel', fuelreg, u, density=19, density_units='g/cm^3', importance=('neutron', 1)
   # graveyrad
   grave = Cell('grave', gravereg)
   # register cells with system
   system.add_cell(fuel, grave)
   # source: default params
   source = Criticality()
   # add source to simulation
   sim.add_source(source)
   # write input file
   MCNPInput(sim).write('filename')
