.. _source_sampling:

=============================
Source Sampling
=============================

.. currentmodule:: pyne.source_sampling

This page presents a detailed explanation of how :py:mod:`Source Sampling` works, including
underlying theory, any assocaited assumptions, and information about correct and
appropriate use of this physics. 

******
Theory
******


**************
Implimentation
**************

Mathematical Details, some example math syntax:

You can use LaTeX format in rst. To do inline expressions (what $expression$ would
do in LaTeX) :math:`A\vec{x} = \vec{b}`.

To have blocks of mathematical content, do this

.. math::
    
    A\vec{x} = \vec{b}

Support is limited to a subset of LaTeX math by the conversion required for many output formats.

***********
Assumptions
***********

The Sampler class chooses the (x, y, z) position within a mesh volume element
with no regard for what geometry cell it lies in. Cell rejection must be
implimented on within the physic-code-specific wrapper script.

**********************
Sample Calculations
**********************
This section provides the sample calculations to justify the results in the nosetests: test_uniform, test_bias, test_bias_spacial.

Consider a mesh with two mesh volume elements with volumes (250, 750). The source on the mesh has two energy groups. The source density distribution is:

.. math::
     q = ((2, 3), (1, 4))

The unnormalized PDF is found by multipling by the volumes:

.. math::
     p_un = ((500, 750), (750, 3000))

Normalizing the PDF yeilds the analog pdf:

.. math::
     p = ((0.1, 0.15), (0.15, 0.6)

Case 1: Uniform Sampling
------------------------
For uniform sampling the biased source density disribution is:

.. math::
     \hat{q} = ((1, 1), (1, 1))

The unnormalized PDF is found by multipling by the volumes:

.. math::
     p_un = ((250, 250), (750, 750))

Normalizing the PDF yeilds the bias pdf:

.. math::
     p = ((0.125, 0.125), (0.375, 0.375)

The weights of particle born from these phase space bins should then be the
ratio of the unbiased to biased PDF values:

.. math::
     weights = ((0.8, 1.2), (0.4, 1.6))

Case 2: User-Specified Biasing
------------------------------
Now consider some user-speficied bias source density distribution:
.. math::
     \hat{q} = ((8, 3), (5, 8))

The unnormalized PDF is found by multipling by the volumes:

.. math::
     p_un = ((2000, 750), (3750, 6000))

Normalizing the PDF yeilds the bias pdf:

.. math::
     p = ((0.16, 0.06), (0.3, 0.48)

The weights of particle born from these phase space bins should then be the
ratio of the unbiased to biased PDF values:

.. math::
     weights = ((0.65, 2.5), (0.5, 1.25))

**********
References
**********

All of the justification for us listening to you.



