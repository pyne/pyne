""" 
This script demostrates the ray tracing capabilities of
dagmc.discretize_geom(), using the Utah Teapot as an example. A .step file of
the model is borrowed from GradCad user Michael Christensen, found here:

https://grabcad.com/library/ceramic-tea-set-teapot-teacup-1

This .step file was then opened in CubIt. A bounding box volume was added, then
the geometry was imprinted and merged and saved as a .sat file with attributes.
The file was then faceted using dagmc_preproc, an executable found in the /bin
directory of a standard MOAB distribution. The faceted file, teapot.h5m, is
found in this directory (pyne/examples).

This script reads in teapot.vtk, calls dagmc.discretize_geom() and tags the
volume fraction output a MOAB mesh. Output is saved as "teapot.vtk", which can
be visualized in VisIt.

This script takes ~2 mins to run. The runtime can be reduced to near
instanteous by changing the num_divs and num_rays varible to decrease the
resolution or the rays fired per mesh row.
"""
from pyne import dagmc
from pyne.mesh import Mesh
import numpy as np

# Number of division in the x, y, and z directions.
num_divs = 50
# Number of rays fired per mesh row.
num_rays = 30

# Geneate superimposed Mesh.
coords1 = np.linspace(-6, 6, num_divs)
coords2 = np.linspace(-0, 7, num_divs)
coords3 = np.linspace(-4, 4, num_divs)
mesh = Mesh(structured=True, structured_coords=[coords1, coords2, coords3])

# Load file and fire rays. Note that DAGMC instances are singleton.
dagmc.load("teapot.h5m")
results, uncs = dagmc.discretize_geom(mesh, num_rays)

# Tag results onto Mesh.
ves = mesh.structured_iterate_hex('xyz')
tag_1 = mesh.mesh.createTag("vol_1", 1, float)
tag_2 = mesh.mesh.createTag("vol_2", 1, float)

for i, ve in enumerate(ves):
    if 2 in results[i].keys():
       res = results[i][2]
    else:
       res = 0
    tag_2[ve] = res

    if 1 in results[i].keys():
       res = results[i][1]
    else:
       res = 0
    tag_1[ve] = res

mesh.mesh.save("teapot.vtk")
