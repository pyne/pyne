.. _pyne_origen22:

******************************************
ORIGEN 2.2 Support -- :mod:`pyne.origen22`
******************************************
Pyne contains support for reading, writing, and merging 
certain ORIGEN 2.2 input and output files (tapes). 
More information may be found in the ORIGEN users manual.

.. currentmodule:: pyne.origen22

All functionality may be found in the ``origen22`` package::

 from pyne import origen22

Examples of use may be found in the user's guide (:ref:`usersguide_origen22`).


-----
TAPE4
-----

.. autofunction:: write_tape4(mat, outfile="TAPE4.INP")


-----
TAPE5
-----

.. autofunction:: write_tape5_irradiation(irr_type, irr_time, irr_value, outfile="TAPE5.INP", decay_nlb=(1, 2, 3), xsfpy_nlb=(204, 205, 206), cut_off=1E-10, out_table_nes=(False, False, True), out_table_laf=(True,  True,  True), out_table_num=None)



-----
TAPE6
-----

.. autofunction:: parse_tape6(tape6="TAPE6.OUT")



-----
TAPE9
-----

.. autofunction:: parse_tape9(tape9="TAPE9.INP")

.. autofunction:: merge_tape9(tape9s)

.. autofunction:: write_tape9(tape9, outfile="TAPE9.INP", precision=3)


------------------------
Other origen22 functions
------------------------

.. autofunction:: sec_to_time_unit(s)

