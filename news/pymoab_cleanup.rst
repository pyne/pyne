**Added:** 

* Testing PyNE build without PyMOAB using Python 2 & 3
* Testing PyNE with PyMOAB (and DAGMC) using Python 2
* Added FindDAGMC.cmake file
* Dockerfile to build many variations of PyNE docker image, with python script CLI

**Changed:** 

* PyNE can be built with PyMOAB whitout DAGMC with limited capabilities

* "--dagmc" flag added to `setup.py` in order to build PyNE against DAGMC

* pyne.mesh now takes advantage of PyMOAB instead of PyTAPS:
   - IMeshTag changed to NativeMeshTag, with according tagetype name change:
     from 'imesh' to 'nat_mesh'
   - write_hdf5(self, filename) -> write_hdf5(self, filename, write_mats)
   - new save(self, filename, write_mats) (alias for write hdf5)
   - new class MeshSetIterator()
   - new get_tag(self, tag_name) and delete_tag(self, tag_name) methods
   - when tagging the root set of a mesh, a new syntax is available:
         - `mymesh.mytag[mesh.mesh.getRootSet()] = val`  can now be written as `mymesh.mytag[mymesh] = val`
   - direct call to the mesh entities change accordingly for example:
      - getEntSets() -> get_entities_by_type( , )
      - getTagHandle('XXX') -> tag_get_handle(types.XXXXX)
      - iterate() -> mesh_iterate()
      - getAllTags(xx) -> tag_get_tags_on_entity(xx)
      - mesh.destroyTag(self, boolean) -> mesh.delete_tag(self)
      - ... (see PyTAPS and PyMOAB respective documentation)
   - those changes have been propagated in mcnp.py, alara.py, ccc.py, dagmc.pyx,
     r2s.py, variancereduction.py, expand_tags.py, and their respective tests... 

* test_fluka:
   - added test to check the data tag name of the different tally part and
     error.

* dagmc_bridge: added a static DagMC instance

* utils.py: updated the download timeout time to 120sec (from 30sec)

* updated CI to use CircleCI 2.1 workflows: now build separately from tests with state saved between runs

**Deprecated:** None

**Removed:** None

**Fixed:** None

**Security:** None
