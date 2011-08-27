.. _install:

============
Installation
============

------
Binary
------
A binary distribution of PyNE is hopefully coming soon.  Until then, please
install from source.


.. _install_source:

------
Source
------
Installing PyNE from source is a two-step process::

    cd pyne/
    python setup.py install
    nuc_data_make

The ``setup.py`` command compiles and installs the PyNE source code.
The ``nuc_data_make`` builds and installs a database of nuclear data.
Unfortunately, this must be done as a second step because most nuclear 
data is under some form of license restriction or export control which 
prevents the developers from distributing it with PyNE.  However, the 
``nuc_data_make`` program (which is installed by ``setup.py``) will
do its best to find relevant nuclear data elsewhere on your machine
or from public sources on the internet.  

********************
Linux (Ubuntu) + EPD
********************
Assuming you are on some flavor of Linux and you primarily use Python 
through the Enthought Python Distribution (EPD), you can install PyNE
to be based off of the EPD packages.

First, you'll need to know where on your system EPD is installed.
Call this variable ``EPD_DIR``; for example if you have installed it 
to your home directory then ``EPD_DIR=$HOME/epd``.  You'll then need
to add the following lines to your :file:`~/.bashrc` file *after* 
installing EPD but *prior to* installing PyNE:

.. code-block:: bash

    export PATH=$EPD_DIR/bin:$HOME/.local/bin:$PATH
    export CPATH=$EPD_DIR/include:$CPATH
    export LD_LIBRARY_PATH=$EPD_DIR/lib:$LD_LIBRARY_PATH

Or as in the example:

.. code-block:: bash

    export PATH=$HOME/epd/bin:$HOME/.local/bin:$PATH
    export CPATH=$HOME/epd/include:$CPATH
    export LD_LIBRARY_PATH=$HOME/epd/lib:$LD_LIBRARY_PATH

You may now proceed with the PyNE install :ref:`as above <install_source>`.
