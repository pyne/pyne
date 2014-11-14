import numpy as np
from itaps import iMesh

from pyne.mcnp import Meshtal, Wwinp
from pyne.mesh import Mesh, IMeshTag


print("reading_meshtal")
meshtal = Meshtal("seed_meshtal")
m = meshtal.tally[4]

m.n_total = n_total = IMeshTag(1, float, mesh = m, name = "neutron_result_total")
m.ww_n = IMeshTag(1, float)

print("finding max")
max_flux = np.max(m.n_total[:])

print("divding my max")
fluxes = m.n_total[:]

ww = []
for flux in fluxes:
    ww.append(flux/(2*max_flux))

m.ww_n[:] = ww

print("rootSet creation")
root_tag = m.mesh.createTag("n_e_upper_bounds",1, float)
root_tag[m.mesh.rootSet] = 100 # 100 MeV neutrons, not likely.

print("creating wwinp from mesh and writing")
wwinp = Wwinp()
wwinp.read_mesh(m.mesh)
wwinp.write_wwinp("seed_wwinp")

print("reading wwinp file and and making vtk")
new_wwinp = Wwinp()
new_wwinp.read_wwinp("seed_wwinp")
new_wwinp.mesh.save("seed_wwinp.vtk")

