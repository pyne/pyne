"""
Tools to generate, set and check the hashes of datasets in pyne.
"""
import hashlib
from pyne.utils import QA_warn

import numpy
import tables

from .. import data

QA_warn(__name__)

# list of nodes from distinct data sets
nodelist = [
    "/atomic_mass",
    "/material_library/materials",
    "/neutron/eaf_xs",
    "/neutron/scattering_lengths",
    "/neutron/simple_xs",
    "/decay",
    "/dose_factors",
]


def check_hashes(nuc_data):
    """
    This function checks the hash of all the nodes in nodelist against the
    built-in ones

    Parameters
    ----------
    nuc_data : str
        path to the nuc_data.h5 file

    """
    check_list = []
    for item in data.data_checksums:
        res = calc_hash(item, nuc_data) == data.data_checksums[item]
        if res is False:
            print("Expected hash: " + str(data.data_checksums[item]))
            print("I got:" + str(calc_hash(item, nuc_data)))
        check_list.append([item, res])
    return check_list


def set_internal_hashes(nuc_data):
    """
    This function sets internal hashes for all the nodes in nodelist.

    Parameters
    ----------
    nuc_data : str
        path to the nuc_data.h5 file

    """
    for item in nodelist:
        set_hash(item, nuc_data)


def check_internal_hashes(nuc_data):
    """
    This function checks the hashes of the nodes in nodelist against internally
    saved ones.

    Parameters
    ----------
    nuc_data : str
        path to the nuc_data.h5 file
    """
    check_list = []
    for item in nodelist:
        res = check_hash(item, nuc_data)
        check_list.append([item, res])
    return check_list


def calc_hash(node, nuc_data):
    """
    This function calculates the hash of a dataset or group of datasets in a
    hdf5 file.

    Parameters
    ----------
    node : str
        String with the hdf5 node name
    nuc_data : str
        path to the nuc_data.h5 file

    """
    with tables.open_file(nuc_data) as f:
        node = f.get_node(node)
        if type(node) == tables.group.Group:
            mhash = hashlib.md5()
            for item in node:
                if type(item[:]) == numpy.ndarray:
                    mhash.update(item[:].data)
                else:
                    if type(item[0]) == numpy.ndarray:
                        for tiny_item in item:
                            mhash.update(tiny_item.data)
                    else:
                        for tiny_item in item:
                            mhash.update(str(tiny_item).encode())
            ret = mhash.hexdigest()
        else:
            ret = hashlib.md5(node[:].data).hexdigest()
    return ret


def set_hash(node, nuc_data):
    """
    This function sets the hash of a dataset or group of datasets in an hdf5
    file as an attribute of that node.

    Parameters
    ----------
    node : str
        String with the hdf5 node name
    nuc_data : str
        path to the nuc_data.h5 file

    """
    the_hash = calc_hash(node, nuc_data)
    with tables.open_file(nuc_data, mode="a") as f:
        f.set_node_attr(node, "hash", the_hash)


def check_hash(node, nuc_data):
    """
    This function checks the hash of a dataset or group of datasets and checks
    it against the stored hash attribute.

    Parameters
    ----------
    node : str
        String with the hdf5 node name
    nuc_data : str
        path to the nuc_data.h5 file

    """
    with tables.open_file(nuc_data) as f:
        hash_val = f.get_node_attr(node, "hash")
    return calc_hash(node, nuc_data) == hash_val
