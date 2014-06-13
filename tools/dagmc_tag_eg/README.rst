DAGMC TAG EXAMPLE
====================

This is a simple example that will add a tag named "density" and assign a
value for that tag to one of the volumes.

In these models, each volume is represented by a meshset that contains no
entities, but does have, as its children, the surfaces that define the volume.

Building
---------

To build this example, you must already have MOAB installed, presumably in
some location `MOAB_INSTALL_PREFIX`.

``
MOAB_DIR=$MOAB_INSTALL_PREFIX
``

Running
----------

To run the example, you need to select an input file, a volume ID (integer) to
tag with a density, a density value (double), and an output filename.

``
./dagmc_tag_eg in_filename vol_id density out_filename
``


