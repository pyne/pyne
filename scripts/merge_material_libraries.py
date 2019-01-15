from pyne.material import MaterialLibrary

"""

This script allows to merge multiple PyNE material hdf5 file into a single one.
"""

def merge_material_library(merged_libname, matlibs, datapaths=[], nucpaths=[]):
    """ Merge the different hdf5 PyNE material libraries into a single one and
    write the merged library in hdf5 format

    Parameters:
    merged_libname (str): name of the new library
    matlibs ([str]): list of the material library name to be merged
    datapaths ([str] -- optional ): list of the datapath for each library to
    merge, using \"/materials\" as default
    nucpaths ([str] -- optional ): list of the nucpath for each library to
    merge, using \"/nucid\" as default

    """


    if len(datapaths) == 0:
        print("No datapaths provided, using \"/Materials\" as default.")
        for i in range(0, len(matlibs)):
            datapaths.append("/Materials")
    elif len(datapaths) != len(matlibs):
        print("!Error! Number of matlibs does not match number of datapaths")

    if len(nucpaths) == 0:
        print("No nucpaths provided, using \"/nucid\" as default.")
        for i in range(0, len(matlibs)):
            nucpaths.append("/nucpaths")
    elif len(nucpaths) != len(matlibs):
        print("!Error! Number of matlibs does not match number of nucpaths")

    merged_mat_library = MaterialLibrary()

    for index, (filename, datapath, nucpath) in enumerate(zip(matlibs, datapaths, nucpaths)):
        mat_lib = MaterialLibrary()
        mat_lib.from_hdf5(filename, datapath, nucpath)
        for item in mat_lib.items():
            merged_mat_library.__setitem__(item[0], item[1])

    merged_mat_library.write_hdf5(
        merged_libname, datapath="/materials", nucpath="/nucid")


