Building Conda Binaries
=======================

The basic strategy for building conda binaries is to use the recipe at
https://github.com/conda/conda-recipes to build pyne and change the yaml
file to the desired branch/tag.

Get the latest conda-recipes

.. code-block:: bash

  $ git clone https://github.com/conda/conda-recipes.git
  $ cd conda-recipes

Now edit the "url:" line and the version line to match the desired tag and
version. After this build pyne, install pyne, run nuc_data_make and run the
tests. In order to use the steps shown below you need jinja2, conda-build,
binstar and nose installed.

.. code-block:: bash

  $ conda build pyne
  $ pyneout=$(conda build --output pyne)
  $ conda install $pyneout
  $ nuc_data_make
  $ cd ..
  $ git clone https://github.com/pyne/pyne.git
  $ cd pyne/tests
  $ nosetests

If all expected tests pass the binary should now be ready to upload to binstar.
If you are uploading to your personal repo (a good idea so you can test on a
clean VM before bricking the primary pyne binary) don't specify a user.

.. code-block:: bash

  $ binstar upload $pyneout

In order to upload to the primary pyne binary repo specify the user pyne (you
need the proper priveleges to do this)

.. code-block:: bash

  $ binstar upload -u pyne $pyneout

You may need to add the --force option depending on if binstar has added support
for multiple simultaneous versions (currently not available).

The strategy is similar on windows but you need to:

1. Install your own copy of cmake
2. Use mingw-get to get mingw32-make
3. Add just mingw32-make to your path (Do not install gcc with mingw-get!)
4. Use the pynewin branch from @crbates https://github.com/crbates/pyne/tree/pynewin
5. Download and install the hdf5 1.8.11 shared library version (VS9 or VS10
   32-bit)
6. Add the hdf5 dll's to your path
7. Write a bld.bat script to build pyne and copy the necessary hdf5 libraries
8. Edit the pyne meta.yaml file to include libpython and mingw and remove hdf5
   from both the build and run steps.
9. build, test, and upload pyne
