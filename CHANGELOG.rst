===============
pyne Change Log
===============

.. current developments

Next Version
====================


v0.7.0 RC2
====================

**New Capabilities**

* Subvoxel R2S allows activation and photon source sampling by material within a mesh voxel
   * SubMode parameter in source_sampling
   * Add code for submode: SUBVOXEL
   * Test functions to test SUBVOXEL codes in test_source_sampling.py
   * Add subvoxel r2s workflow support in pyne/scripts/r2s.py
   * Use new sampler for both voxel and sub-voxel R2S.
   * Define cell_mats for subvoxel equals False
   * Function Sampler::mesh_tag_data to be compatible for both DEFAULT and SUBVOXEL mode

* OpenMC support
   * Add CI test for OpenMC as an optional dependency
   * Add options to install OpenMC Python API
   * Functions to read OpenMC state point file to get meshtally data

* File scripts/activation_responses.py to get activation responses data for ALARA output file for better visualization. Responses include:
   * decay heat
   * specific activity
   * alpha heat
   * beta heat
   * gamma heat
   * wdr
   * photon_source
   * multiple response in one output.txt file

* first introduction of GT-CADIS workflows
   * added gtcadis.py script
   * first step for the GT-CADIS workflow, further steps to follow

**Enhancements to Previous Capabilities**

* Enhancements of Material and MaterialLibrary capabilities
   * A C++ API for the MaterialLibrary class has been created for direct
     use in compiled software tools
   * Change the MaterialLibrary Python API to bind to the new c++ one.
   * Overhaul of HDF5 format of Material and MaterialLibrary
      * Capability to read nuclides from a specific path when loading material
        from an hdf5 file (not necessarily relying on autodetect the nucpath).
      * change the structure of material when written in a hdf5 file:
        ::

          /material/
          /material/my_mat/
          /material/my_mat/composition
          /material/my_mat/nuclidelist
          /material/my_mat/composition_metadata
        
        where, ``/material`` and ``/material/my_mat`` are now HDF5 groups
   * Material I/O formats
      * pyne::Material: capability to form PHITS input material card
      * capability to form a uwuw material name
      * A new argument ``mult_den`` (bool) that makes *Material* class's
        **get_density_frac** method (hence also **mcnp** and **phits** methods)
        write atom/mass density if true, otherwise write atom/mass fraction.
      * increased precision of the material density reported in the comment card
        for an MCNP material to 5 decimal places
      * Adds tests to tests/test:materials.py:test_decay_heat(self) to check more isotopes

* Bug fixes in serpent support
   * in serpent.py, function VOL_serp2_fix() to correct for
     _VOL variable being an array. as seen in serpent 2
   * fix error in serpent parse_dep
   * serpent.py function parse_dep.  Catches ValueError that
     occurs when parse_dep attempts to make_mats with a serpent 2 \*_dep.m file
     due to the \*_VOL variable not being a float

* Changes in source sampling for mesh-based Monte Carlo sources
   * Add statistics summary output of find_cell failure in source sampling.
   * Add the ability to allow user turn off the void rejection in source sampling.
   * Add cell_fracs and cell_number tags for both default and subvoxel r2s modes
   * Check for the existence of the e_bounds file. Print error message when it's missing.
   * Check for bias_tag data. Report error when bias tag data are all zero
   * Check 'cell_fracs' tag in source_sampling.cpp when sub_mode is DEFAULT. Prevent wrong use of source.h5m.
   * Fix the problem of reading cell_number_tag with size of 1
   * Change mode range of cell rejection from >3 to >2
   * Sort cell_fracs according to the order of 'idx' and 'vol_frac'. For faster source sampling.
   * Pass cell_list back to Fortran, to speed up source sampling.
   * function to write total photon intensities for subvoxel r2s
   * Removed variables ```icl_tmp``` and ```find_cell``` which are not longer needed.
   * MCNP6 version of source.F90
   * Changed source.F90 to use "implicit none" instead of "implicit real"
   * Addition & updates of unit tests for above improvements

* Improvements in Rigorous-2-Step shutdown dose rate analysis workflow
   * Documentation improvements
   * Provide example files for variety of problems/problem modes
   * Improvements in testing of R2S
      * Use example files for automated testing
   * Streamline code related to addition of subvoxel mode
      * Combine the subvoxel/voxel R2S loops to calculate the total photon source intensities.
      * Keep cell_number, cell_fracs, cell_largest_frac_number and cell_largest_frac tag in r2s step1
      * Use subvoxel and normal r2s compatible workflow parameters
      * Input check of cell_fracs tag under voxel mode. As the cell_fracs tag is there for voxel/sub-voxel mode.
   * Load geom and calculate cell_mats in r2s step2
   * Read decay times from r2s config.ini, and then write them into alara_inp.
   * In R2S step2, add option to write only 'total' to h5 file, reduce the CPU time
   * Error in voxel R2S.
   * Changes in processing of ALARA input/output
      * Change some default names of alara_inp.
      * Decay times in the alara_params.txt.
      * Add input units check to the function utils.py/to_sec
      * Use function utils.py/to_sec to replace alara.py/_TO_SEC
      * Simplify the method to get the list of decay/cooling times

* Nuclear Data Handling and Reporting
   * Fixed issue where some gamma x-rays where throwing ``NotANuclide`` errors
     because the underlying nuclides were being read & recorded with negative ids.
     All nuclide ids are now ensured to be positive.
   * Misidentification of descriptive text in (MF,MT)=(1,451) as contents lines.
   * decay_heat() in material.cpp now calls metastable_id to convert zas_id to state_id
   * Fix ENDF parsing of TSL files with short collision time approximation for non-principal atoms.
   * endf.Library._read_headers() and regular expressions in endf.pyx
        * Removed regexps: CONTENTS_R, SPACE66_R, NUMERICAL_DATA_R
        * Added regexps:   SPACEINT11_R
        * Added methods:   _isContentLine(parts)
   * ENSDF database link to 2019 Oct 4th database
   * Update the C012-n.ace file link.
   * Missing elements name_to_zz dictionary
   * Updated half_life in data.pyx to return nan if isotope not found (#1257)

* Improvements in Mesh capabilities
   * added mesh tally definitions to tallies
   * store multi particle tally (for Volume and Surface)
   * mcnp can write multi-particle tally
   * Move check of tag_names to mesh.py
   * Fix a problem of creating mesh from reading h5m files in unstructued R2S
   * Default initializer pyne.mesh.Mesh() now raises an exception with info on how
     to properly make a mesh
   * Move class MeshTally from mcnp.py to mesh.py
   * Change the method of creating meshtally from mcnp meshtal
   * pyne.mesh now takes advantage of PyMOAB instead of PyTAPS:
      * IMeshTag changed to NativeMeshTag, with according tagetype name change:
        from 'imesh' to 'nat_mesh'
      * write_hdf5(self, filename) -> write_hdf5(self, filename, write_mats)
      * new save(self, filename, write_mats) (alias for write hdf5)
      * new class MeshSetIterator()
      * new get_tag(self, tag_name) and delete_tag(self, tag_name) methods
      * when tagging the root set of a mesh, a new syntax is available:
         * `mymesh.mytag[mesh.mesh.getRootSet()] = val`  can now be written as `mymesh.mytag[mymesh] = val`
      * direct call to the mesh entities change accordingly for example:
         * getEntSets() -> get_entities_by_type( , )
         * getTagHandle('XXX') -> tag_get_handle(types.XXXXX)
         * iterate() -> mesh_iterate()
         * getAllTags(xx) -> tag_get_tags_on_entity(xx)
         * mesh.destroyTag(self, boolean) -> mesh.delete_tag(self)
         * ... (see PyTAPS and PyMOAB respective documentation)
      * those changes have been propagated in mcnp.py, alara.py, ccc.py, dagmc.pyx,
        r2s.py, variancereduction.py, expand_tags.py, and their respective tests...

**Maintenance**

* Documentation Changes
   * Credit Rochman for allowing redistribute TENDL file
   * Fix various typos
   * automatic deployment of a updated version of the website on tags
   * automatic creation of a new version of the website (not deployed) for
     verification purposes in ``pyne.github.com/website_preview``
   * New developers guide: The update adds information about creating an environment,
     updates formatting for more consistency, details considerations and methods to
     check the version of dependencies, and adds additional links to coding resources.
   * In website index, change C++ API link to "C++ API Documentation"
     instead of "C++ & Fortran API Documentation"
   * Added publications to bibliography (PR #1256)
   * Adding contributing guide and code of conduct (#1258)
   * Changed Doc and Tutorial mentions of iPython notebooks to Jupyter notebooks (PR #1262)
    
* Improvements in building and testing
   * require contributor to change CHANGELOG
   * Expand testing matrix to include:
      * python 2 vs 3
      * with vs without PyMOAB
      * with vs without DAGMC
   * Added FindDAGMC.cmake file
   * turn off BLAS/LAPACK & FORTRAN in MOAB build
   * Dockerfile to build many variations of PyNE docker image, with python script CLI
   * Add hdf5-tools as dependency for docker images used in CircleCI, for better nose test comparing h5 files
   * Add future as dependency for docker images used in CircleCI, for python2 and python3 compatibility
   * "--dagmc" flag added to ``setup.py`` in order to build PyNE against DAGMC
   * new check won't now be triggered after a merge only on PRs
   * utils.py: updated the download timeout time to 120sec (from 30sec)
   * updated CI to use CircleCI 2.1 workflows: now build separately from tests with state saved between runs
   * test_fluka:
      * added test to check the data tag name of the different tally part and
        error.
   * revert internal nuc_data_path to origin value after internal data test
   * added DEFINE variable to allow material.cpp amalgamation without decay.cpp
   * now skips endf test when website is not reachable to allow completeness of the other tests.
      * test file for ENDF was wrong
   * Add functions to do file, file block, line, and string almost the same
     compare functions in pyne/utils.py
   * make data available as replacement for data.pyne.io (#1261, #1265)
   * Removed iPython check from ``setup.py`` and added Jupyter to be an optional dependency in documentation (#1273)

* Code cleanup
   * Formatting improvements
   * Compatibility with language updates
      * update the way that ``collections`` is imported in preparation
        for deprecated changes in future python versions
      * removed some imports of ``collections`` that were not necessary
      * change return type of method to avoid compiler compatibility issue
      * Convert some code and tests to enable python2/3 compatibility
   * Clean up some hard coded strings in test_source_sampling.py
   * ``rxname.child()`` and ``rxname.parent()`` now accept ``str`` for the
     ``z`` argument in Python 3.
   * dagmc_bridge: added a static DagMC instance
   * cleanup throws return from ``return (const char *)`` to simple ``return`` 
     (it was suggested that those complicated return might cause seg fault, on some system -- OsX+conda)


v0.5.11
====================

**Added:**

* Function to convert unit to s in pyne/alara.py
* Function to do float match for decay times
* Add SourceParticle class in pyne/src/source_sampling.
* Codes to read ALARA output file under subvoxel R2S condition
* A function to build up a subvoxel_array from mesh and cell_mats information
* A test function to test the process of reading ALARA output file
* Test function for subvoxel with (N, 1) condition in test_mesh.py
* Reshape the array when max_num_cells == 1

**Changed:**

* shape of IMeshTag when input value is a (N, 1) array
* set tag as array rather than number
* decaygen now gets the include dir based on the compiler path.
* Build system now explitily looks for C++11 standard compatability.
* Unit of e_bounds changed from eV to MeV
* Change loop variables to be v for volume elements and e for energy groups (instead of i & j)
* Use bias_mode instead of mode to allow for additional mode types in future
* A parameter in test_alara.py, to test modified match method
* Correct the wrong mode description comment in pyne/src/source_sampling.h
* Change the particle_birth return value from std::vectot<double> to SourceParticle object
* Some code clean up
* Some clean up of white space

**Removed:**

* Code in mesh.py to reshpe a (N,1) to (N, ) array is no longer needed if PR #971 merged

**Fixed:**

* decaygen now can properly produce Clang assembly.
* Build system would always download cram sources, even if they already existed.
  This has been fixed.
* ENDF error bounds bug that was preventing ``nuc_data_make`` from working.
* NNDC no longer provides the ``mednew.dat`` data set. A fallback has been
  supplied.


v0.5.10
====================

**Fixed:**

* Made SSL context creation Python 2 & 3 Compatible.


v0.5.9
====================

**Changed:**

* Downloading files now uses null SSL context.


v0.5.8
====================

**Changed:**

* Downloading data now uses HTTP, rather than HTTPS.


v0.5.7
====================

**Fixed:**

* Occassional bug with downloading URL fix.


v0.5.6
====================


v0.5.5
====================


v0.5.4
====================
