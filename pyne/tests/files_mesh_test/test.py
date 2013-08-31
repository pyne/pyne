from itaps import iMesh, iBase, iMeshExtensions
from pyne.mesh import Mesh

a = iMesh.Mesh()
a.load("grid543.h5m")
for ent_set in a.getEntSets():
    print a.getTagHandle("BOX_DIMS")[ent_set]

