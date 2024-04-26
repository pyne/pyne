# creates a statepoint file for the CI test suite

import openmc

mesh = openmc.RegularMesh()
mesh.lower_left = (-40.0, -12.5, -2.5)
mesh.upper_right = (80.0, 25.0, 5.0)
mesh.dimension = [3, 2, 1]

my_material = openmc.Material()
my_material.add_element("H", 2, percent_type="ao")
my_material.add_element("O", 1, percent_type="ao")
my_material.set_density("g/cm3", 1)
my_materials = openmc.Materials([my_material])

outer_surface = openmc.Sphere(r=100, boundary_type="vacuum")
cell_1 = openmc.Cell(region=-outer_surface, fill=my_material)
my_geometry = openmc.Geometry([cell_1])
my_settings = openmc.Settings()
my_settings.batches = 10
my_settings.particles = 10000
my_settings.run_mode = "fixed source"

my_source = openmc.IndependentSource()
my_source.space = openmc.stats.Point((0, 0, 0))
my_source.angle = openmc.stats.Isotropic()
my_source.energy = openmc.stats.Discrete([14e6], [1])
my_settings.source = my_source

energy_filter = openmc.EnergyFilter([0.0, 1.0e6, 2.0e6])
mesh_filter = openmc.MeshFilter(mesh)
mesh_spectra_tally = openmc.Tally(name="mesh_spectra_tally")
mesh_spectra_tally.scores = ["flux"]
mesh_spectra_tally.filters = [mesh_filter, energy_filter]
my_tallies = openmc.Tallies([mesh_spectra_tally])

model = openmc.model.Model(my_geometry, my_materials, my_settings, my_tallies)
results_filename = model.run()
