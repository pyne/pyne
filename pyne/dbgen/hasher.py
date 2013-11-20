"""
Low level tools to check hashes of datasets in pyne.

Author: crbates
"""
import tables
import numpy as np
import hashlib


def calc_atomic_decay_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        return(hashlib.md5(f.root.atomic_decay[:].tostring()).hexdigest())    
def set_atomic_decay_hash(nuc_data):
    the_hash = calc_atomic_decay_hash(nuc_data)
    with tables.openFile(nuc_data,mode='a') as f:
        f.root.atomic_decay.attrs.hash = the_hash
def check_atomic_decay_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        return calc_atomic_decay_hash(nuc_data) == f.root.atomic_decay.attrs.hash
        
def calc_atomic_weight_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        return(hashlib.md5(f.root.atomic_weight[:].tostring()).hexdigest())
def set_atomic_weight_hash(nuc_data):
    the_hash = calc_atomic_weight_hash(nuc_data)
    with tables.openFile(nuc_data,mode='a') as f:
        f.root.atomic_weight.attrs.hash = the_hash
def check_atomic_weight_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        return calc_atomic_weight_hash(nuc_data) == f.root.atomic_weight.attrs.hash    

    
def calc_neutron_scattering_lengths_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        return(hashlib.md5(f.root.neutron.scattering_lengths[:].tostring()).hexdigest())
def set_neutron_scattering_lengths_hash(nuc_data):
    the_hash = calc_neutron_scattering_lengths_hash(nuc_data)
    with tables.openFile(nuc_data,mode='a') as f:
        f.root.neutron.scattering_lengths.attrs.hash = the_hash
def check_neutron_scattering_lengths_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        return calc_neutron_scattering_lengths_hash(nuc_data) == f.root.neutron.scattering_lengths.attrs.hash   
    
def calc_material_library_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        mhash = hashlib.md5()
        for item in f.root.material_library:
            mhash.update(str(item[:]))
        return mhash.hexdigest()
def set_material_library_hash(nuc_data):
    the_hash = calc_material_library_hash(nuc_data)
    with tables.openFile(nuc_data,mode='a') as f:
        f.root.material_library._v_attrs.hash = the_hash
def check_material_library_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        return calc_material_library_hash(nuc_data) == f.root.material_library._v_attrs.hash
    
    
def calc_neutron_eaf_xs_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        mhash = hashlib.md5()
        for item in f.root.neutron.eaf_xs:
            mhash.update(np.array(item[:]).tostring())
        return mhash.hexdigest()    
def set_neutron_eaf_xs_hash(nuc_data):
    the_hash = calc_neutron_eaf_xs_hash(nuc_data)
    with tables.openFile(nuc_data,mode='a') as f:
        f.root.neutron.eaf_xs._v_attrs.hash = the_hash
def check_neutron_eaf_xs_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        return calc_neutron_eaf_xs_hash(nuc_data) == f.root.neutron.eaf_xs._v_attrs.hash
    
def calc_neutron_simple_xs_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        mhash = hashlib.md5()
        for item in f.root.neutron.simple_xs:
            mhash.update(item[:].tostring())
        return mhash.hexdigest()
def set_neutron_simple_xs_hash(nuc_data):
    the_hash = calc_neutron_simple_xs_hash(nuc_data)
    with tables.openFile(nuc_data,mode='a') as f:
        f.root.neutron.simple_xs._v_attrs.hash = the_hash
def check_neutron_simple_xs_hash(nuc_data):
    with tables.openFile(nuc_data) as f:
        return calc_neutron_simple_xs_hash(nuc_data) == f.root.neutron.simple_xs._v_attrs.hash
        