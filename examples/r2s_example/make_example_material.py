#!/usr/bin/python
#
from pyne.material import Material,MaterialLibrary
print "Welcome!"
mat_lib=MaterialLibrary()
#
mat2 = Material({'Fe':0.655,'Cr':0.170,'Ni':0.120,'Mo': 0.025,'Mn': 0.02, 'Si':.01},density=7.92)
mat2=mat2.expand_elements()
#
# define a simple water since O-18 not in mcnp xs libs
watervec={10010:2,80160:1} # simple water
water = Material()
water.density = 1.0
water.from_atom_frac(watervec)
#
mat_lib["Steel"]=mat2
mat_lib["Water"]=water
#
mat_lib.write_hdf5("example_material_lib.h5")
#
print "All done!"
