About the dagmc_get_materials and tests scripts
====
# General:
'dagmc_get_materials.py' is a Python script to handle the group names from a CAD model by creating a list of all group names and then extract the materials group names which then are used to get the specifications/metadata of each material (for exapmle the composition, density, atoms per molecule and the like). It then creates an output h5m file containing a directory that contains all the materials from the group names on the original CAD model along with all its metadata. In order to be able to copy the material metadata we need a materials library to copy from so some dependencies are required.
- The script can be run as a python script and all the flags needed can be found using:
  ```python dagmc_get_materials.py --help```

- Generally it is run as:
  ```python dagmc_get_materials.py -d ['path to PyNE materials library nuc_data.h5'] -f ['path of the CAD model .h5m file'] -o ['name of the output file as "NAME.h5m"']```


# Dependencies:
- PyNE, installation instructions can be found here; http://pyne.io/install.html
- PyTAPS, installation instructions can be found here; https://pythonhosted.org/PyTAPS/install.html
- nose tools (https://nose.readthedocs.org/en/latest/), for testing

# Sample output:
-  An output ('sample_output.h5m') obtained by running the script with an arbitrary input h5m file. 
  

# Tests:
- test #1: ('test_functions.py') is a test of the different functions of the 'dagmc_get_materials.py' script. it can be run as ```nosetests test_functions.py```

- test #2:
('test_output.py') tests the output h5m file created by running the script ans is run as a python script or using nosetests  
