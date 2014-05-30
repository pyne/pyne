About the dagmc_get_materials and tests scripts
====
# General:
'dagmc_get_materials.py' is a Python script to handle the groups names from a CAD model by creating a list of all groups names and then extract the materials groups names which then are used to get the specifications/metadata of each material (for exapmle the composition, density, atoms per molecule and the like). It then creates an output h5m file containing a directory that contains all the materials from the group names on the original CAD model along with all its metadata. In order to be able to copy the material metadata we need a materials library to copy from so some dependencies are required.
- mthe script can be run as a python script and all the flags needed can be found using >>python dagmc_get_materials.py --help

- Generally it is run as > python dagmc_get_materials.py -d ['path to PyNE materials library nuc_data.h5'] -f ['path of the CAD model .h5m file'] -o ['name of the output file as "NAME.h5m"']


# Dependencies:
- PyNE, installation instruction can be found on http://pyne.io/install.html
- PyTAPS, installation instruction can be found on https://pythonhosted.org/PyTAPS/install.html
- nose tools (https://nose.readthedocs.org/en/latest/) for testing

# Sample output:
Provided is a sample CAD model ('sample_output.h5m') created using CUBIT with arbitrary material groups.

# Tests:
- test #1 ('test_functions.py') is a test for the different functions on the 'dagmc_get_materials.py' script. it can by running >> nosetests test_functions.py

- test #2 ('test_output.py') tests the output h5m file created by running the script ans is run as a python script or using nosetests  
