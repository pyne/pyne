.. _theorymanual_nucname:

===============================
nucname
===============================

.. currentmodule:: pyne.nucname

The :py:mod:`pyne.nucname` implements a canonical form for nuclide naming that is used
throughout PyNE.  This is called the ``id`` or ``nucid``. This form is
simple but powerful enough to unambiguously handle all nuclide naming conventions
from other codes and conventions. All other nuclide naming schemas are
derived from transformations of this form.

Concisely, the format is three proton number digits (the Z-number) followed
by three nucleon number digits (the A-number) followed by four state digits
(the S-number). Thus we see that the full id form is ``ZZZAAASSSS``.

This format has two important properties:

1. It preserves a natural sort order for nuclides (hydrogen species come before helium
   species which come before lithium, etc.)
2. It maximizes the amount of information that can be stored in a 32-bit signed integer.

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
