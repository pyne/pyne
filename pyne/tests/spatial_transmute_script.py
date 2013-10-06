from pyne import transmute as tm, nucname as nn
import numpy as np
import tables as tb
phi = np.zeros((175,1))
for i in np.arange(phi.shape[0]):
    phi[i] = 1.0E+12
print('Flux initialized.')
nuc = nn.zzaaam('FE56')
inp_1 = {nuc : 1.}
inp_2 = {nuc : 2.}
space = {1 : (phi, inp_1), 2 : (phi, inp_2)}
t_sim = 3153600
tol = 1e-7
print('Begin transmutation.')
space_out = tm.transmute_spatial(space,t_sim,None,tol)
print('space_out keys')
print(space_out.keys())
print('space_out 1')
print(space_out[1])
print('space_out 2')
print(space_out[2])
h5file = tb.openFile('spatial_h5test.h5','w')
group = h5file.createGroup('/','transmute_parentGroup','Transmute Output')
print('Writing to hdf5.')
tm.write_space_hdf5(h5file, group, space_out)
h5file.close()
print('Closed hdf5.')
"""
print('Total densities -')
print('Volume 1: ' + str(sum(space_out[0].values())) + '.')
print('Volume 2: ' + str(sum(space_out[1].values())) + '.')
"""
