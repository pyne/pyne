.. _theorymanual_nucname:

===============================
nucname
===============================

.. currentmodule:: pyne.nucname

The :py:mod:`pyne.nucname` implements a canonical form for nuclide naming that
is used throughout PyNE.  This is called the ``id`` or ``nucid``. This form is
simple but powerful enough to unambiguously handle all nuclide naming conventions
from other codes and conventions. All other nuclide naming schemas are
derived from transformations of this form.

Concisely, the format is three proton number digits (the Z-number) followed
by three nucleon number digits (the A-number) followed by four state digits
(the S-number). Thus we see that the full id form is ``ZZZAAASSSS``.

This format has two important properties:

1. It preserves a natural sort order for nuclides (hydrogen species come before
   helium species which come before lithium, etc.)
2. It maximizes the amount of information that can be stored in a 32-bit signed 
   integer.

In most applications, the state number is interpreted as a metastable state since
this is the most common usage in nuclear engineering. In certain instances when
specified explicitly, the state number may represent the internal excitation
state of the nucleus. This usage arises most frequently in the context of
radioactive decay.

No physics should be performed on the nuclide ``id`` format directly.

-------------------
Rules & Conventions
-------------------
The following are rules about the semantic meaning of specific classes of
identifiers, in no particular order:

* Negative nuclides have no meaning (yet).
* Nuclides where ``AAA < ZZZ`` are physically impossible and thus not allowed.
* Nuclides where ``SSSS == 0`` are considered to be in the ground state.
* Nuclides where ``AAA == 0`` and ``ZZZ > 0`` represent the chemical element.
  In most applications this takes on the natural isotopic abundances. The
  state number is assumed to be zero for elements
* Humans should use a human readable format, computers should use ``id``.

Well-Defined vs Ambiguous Situations
....................................
In situations where the input naming convention is well-defined, it is *highly*
recommended that you use the direct ``<form>_to_id()`` functions (e.g. 
``mcnp_to_id()``) to convert from a nuclide in the given form to the id form 
representation. When a high level of quality assurance is required, it is 
advisable to require an specific input format to leverage the exactness of the 
direct-to-id functions.

However, in situations where arbitrary nuclide naming conventions are allowed, 
you must use the ``id()`` function. An example of such a situation is when accepting 
human input. This function attempts to resolve the underlying nuclide identifier. 
For most nuclides and most normal spellings, this resolution is straightforward. 
However, some nulcides are ambiguous between the various supported naming conventions.
Other nuclide names are completely indeterminate. In the case of indeterminate 
species the ``id()`` function will throw an error. In the case of ambiguous forms 
where the input is an integer input, the form resolution order is:

- id
- zz (elemental z-num only given)
- zzaaam
- cinder (aaazzzm)
- mcnp
- zzaaa

For string (or ``char *``) input, the form resolution order is as follows:
  
- ZZ-LL-AAAM
- Integer form in a string representation, uses integer resolution
- NIST
- name form
- Serpent
- LL (element symbol)

The ``name()`` function first converts functions to id form using the ``id()`` 
function. Thus the form order resolution for ``id()`` also applies to ``name()``.

Nuclide Forms
.............
The following are the currently implemented nuclide naming conventions in pyne:

.. include:: ../nucnameforms.rst
