.. _usersguide_data:

============
Atomic & Isotopic Data
============

.. currentmodule:: pyne.data

.. automodule:: pyne.data

The isotopic & elemental abundance data is usually loaded at the time the atomic_mass function is called. This data now
exists in the C++ implementation. Behind the scenes the nuc_data.h5 file is loaded and the appropriate data elements populated, now
when this function is called, if the nuc_data.h5 does not exist, the data is instead loaded from that stored in a C++ class. This
is not particularly useful for Python users of PyNE, however the pure C++ users can now use the atomic and isotopic data
from the amalgamated source and do not need to carry the nuc_data.h5 with you.

***************************
C++ example Use of Data Class
***************************

The main use of this feature is to allow C++ users to be able to call the abundance and nuclear data functions without the use of nuc_data.h5. For example, 

.. code-block:: c

   #include "pyne.h"
   #include <iomanip>
   
   int main(int argc, char* argv[]) {
     pyne::NUC_DATA_PATH = ""; // for atomic data
     double atomic_mass = pyne::atomic_mass("2H");
     std::cout << std::setprecision(8);
     std::cout << "Atomic mass of deuterium is " << atomic_mass << std::endl;
   }

To compile & link against your installed version of PyNE

.. code-block:: bash
		
   g++ test.cpp -I$HOME/.local/include/pyne -I<path to hdf5>/include -L$HOME/.local/lib/ -L<path to hdf5>/lib -lhdf5 -lpyne -o test

Running this example gives.

.. code-block:: bash
		
   ./test
   Atomic mass of deuterium is 2.0141018
   
***************************
Python example Use of Data Class
***************************

A Python example for loading data is shown below.

.. code-block:: ipython

   In [1]: from pyne.data import atomic_mass

   In [2]: print atomic_mass('2H')
   2.01410177812

If for whatever reason the nuc_data.h5 file is not found or doesn't exist, the above command will still work. You can force the nuc_data.h5m file to be not found as shown in the below example.

.. code-block:: ipython

   In [1]: from pyne.pyne_config import pyne_conf
   
   In [2]: from pyne.data import atomic_mass
   # note, never do this, this is just for testing and this example
   In [3]: pyne_conf.NUC_DATA_PATH = b'some silly path that doesnt exist'

   In [4]: print atomic_mass('2H')
   2.01410177812

