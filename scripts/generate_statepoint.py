import openmc
import os

box = openmc.model.RectangularParallelepiped(-40, 40, -12.5, 12.5, -2.5, 2.5,
                                              boundary_type='vacuum')
box_cell = openmc.Cell(region=-box)
geometry = openmc.Geometry([box_cell])

mesh = openmc.RegularMesh()
mesh.lower_left = (-40, -12.5, -2.5)
mesh.upper_right = (40, 12.5, 2.5)
mesh.dimension = (3,2,1)

energy_bins = [0.0, 1e6, 2e7]

energy_filter = openmc.EnergyFilter(energy_bins)
mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name='flux tally', tally_id=1)
tally.filters = [energy_filter, mesh_filter]
tally.scores = ['flux']
tallies = openmc.Tallies([tally])

settings = openmc.Settings()
settings.particles = 1000
settings.batches = 10
settings.output = {'path':'../tests/files_test_openmc/',
                   'summary': False,
                   'tallies': False}
settings.run_mode = 'fixed source'

model = openmc.Model(geometry=geometry, settings=settings, tallies=tallies)

model.run()

os.remove('model.xml')
