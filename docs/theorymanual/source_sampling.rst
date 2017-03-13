.. _source_sampling:

=============================
Source Sampling
=============================

:Author: Elliott Biondo

.. currentmodule:: pyne.source_sampling

******
Theory
******

Meshes can be used to represent particle source distributions for Monte Carlo
radiation transport. On a mesh, source intensity is discretized spatially (into
mesh volume elements) and by energy (into energy bins). In order to randomly
sample these distributions to select particle birth parameters (position,
energy, statistical weight) a discrete probability density function (PDF) must be
created, which can be sampled with pseudo-random variates. It is convenient to
create a single PDF to describe all of phase space; in other words, each bin
within the PDF represents the probability that a particle is born in a
particular energy group within a particular mesh volume element. 

In pyne, meshes define volumetric source density :math:`q'` with units of
:math:`\frac{particles}{time \cdot volume}`. In order to find the source
intensity of a single phase space bin (of index :math:`n`), the density must be
multiplied by the volume of the mesh volume element:

.. math::
     q(n) = q'(n) \cdot V(n)

The probability :math:`p` that a particle is born into a particular phase space
bin is given by the normalized PDF:

.. math::
     p(n) = \frac{q(n)}{\sum_N{\,q(n)}}

where :math:`N` is the total number of phase space bins (the number of mesh
volume elements and energy groups). Phase-space bins can be selected from this
PDF and all particles will have a birth weight of 1. This is known as analog
sampling.  Alternatively, a biased source density distribution :math:`\hat{q}'`
can be specified yielding a biased PDF :math:`\hat{p}(n)`. Sampling the biased
PDF requires that particles have a statistical weight:

.. math::
     w(n) = \frac{p(n)}{\hat{p}(n)}

Once a phase space bin is selected a position must be sampled uniformly within
the selected mesh volume element to determine the (x, y, z) birth position, and
energy must be uniformly sampled uniformly within the selected energy bin.

**************
Implementation
**************

The Sampler class reads :math:`\hat{q}` and optionally :math:`\hat{q}'` from a
MOAB mesh. PDFs are created using the method described above. In order to efficiently
sample these PDFs an alias table is created [1][2]. This data structure requires an
:math:`O(n^2)` setup step, but then allows for :math:`O(1)` sampling. Monte
Carlo radiation transport typically involves the simulation of :math:`10^{6}`
to :math:`10^{12}` particles, so this expensive setup step is well-justified. 

In the analog sampling mode, an alias table is created from :math:`q`. In the
uniform and user-specified sampling modes, an alias table is created from
:math:`\hat{q}` and birth weights are calculated for each phase space bin. In
the uniform sampling mode, :math:`\hat{q}` is created by assigning a total
source density of 1 to each mesh volume element, so that all space is sampled
equally. Within each mesh volume element, a normalized PDF is created on the
basis of source densities at each energy.  An optional threshold can be provided
such that all mesh volume elements with a total source below that threshold will
remain unbiased under the uniform sampling mode.  The default threshold is 0 so
that all mesh elements are biased.  

The method for uniformly sampling within a mesh volume element of Cartesian mesh
is straightforward. A vertex of the hexahedron (:math:`O`) is chosen and three
vectors are created: :math:`\vec{x}`, :math:`\vec{y}`, and :math:`\vec{z}`.
Each vector points to an adjacent vertex (in the x, y, z, direction
respectively) with a magnitude equal to the length of the edge connecting the
vertex to the adjacent vertex. Three random variates are chosen (:math:`v_1`,
:math:`v_2`, :math:`v_3`) in order to randomly select a position (:math:`P`)
within the hexahedron:

.. math::
     P = O + v_1 \cdot \vec{x} + v_2 \cdot \vec{y} + v_3 \cdot \vec{z}

A similar method is used for uniformly sampling within a tetrahedron, as
described in [3].

***********
Assumptions
***********

The Sampler class chooses the (x, y, z) position within a mesh volume element
with no regard for what geometry cell it lies in. Cell rejection must be
implemented within the physic-code-specific wrapper script. 

**********************
Sample Calculations
**********************

This section provides the sample calculations to justify the results in the
nosetests: test_uniform, test_bias, test_bias_spatial.

Consider a mesh with three mesh volume elements with volumes (3, 0.5, 2). The
source on the mesh has two energy groups. The source density distribution is:

.. math::
     q' = ((2000, 1000), (9000, 3000), (1, 2))

The source intensity is found by multiplying by the volumes:

.. math::
     q = ((6000, 3000), (4500, 1500), (2, 4))

Normalizing yields the analog PDF:

.. math::
     p = ((1000/2501, 500/2501), (750/2501, 250/2501) (1/7503, 2/7503) )

or approximately:

.. math::
     p = ((0.39984, 0.199920), (0.29988, 0.09996), (0.0001333, 0.0002666) )

Case 1: Uniform Sampling
------------------------

For uniform sampling the biased source density distribution is created by
normalizing the source density to 1 within each mesh volume element:

.. math::
     \hat{q}' = ((2/3, 1/3), (3/4, 1/4), (1/3, 2/3))

The biased source intensity is found by multiplying by the volumes:

.. math::
     \hat{q} = ((2, 1), (3/8, 1/8), (2/3, 4/3))

Normalizing yields the biased PDF:

.. math::
     \hat{p} = ((4/11, 2/11), (3/44, 1/44), (4/33, 8/33) )

or appoximately:

.. math::
     \hat{p} = ((0.363636, 0.181818), (0.0681818, 0.0227272), (0.121212, 0.242424) )

The weights of particle born from these phase space bins should then be the
ratio of the unbiased to biased PDF values:

.. math::
     w = ((1.09956, 1.09956), (4.3982, 4.3982), (0.0010996, 0.0010996))

Case 2: Uniform Sampling with Threshold
---------------------------------------

If a threshold of 5 is added to uniform sampling, the biased source density
distribution is created by normalizing the source density to 1 within each
mesh volume element that has a total strength above that threshold:

.. math::
     \hat{q}' = ((2/3, 1/3), (3/4, 1/4), (1 , 2))  hmmmm

The biased source intensity is found by multiplying by the volumes:

.. math::
     \hat{q} = ((2, 1), (3/8, 1/8), (2, 4))

Normalizing yields the biased PDF:

.. math::
     \hat{p} = ((4/11, 2/11), (3/44, 1/44), (4/33, 8/33) )

or appoximately:

.. math::
     \hat{p} = ((0.363636, 0.181818), (0.0681818, 0.0227272), (0.121212, 0.242424) )

The weights of particle born from these phase space bins should then be the
ratio of the unbiased to biased PDF values:

.. math::
     w = ((1.09956, 1.09956), (4.3982, 4.3982), (0.0010996, 0.0010996))

Case 3: User-Specified Biasing
------------------------------
Now consider some user-specified bias source density distribution:

.. math::
     \hat{q}' = ((1, 2), (3, 3))

The biased source intensity is found by multiplying by the volumes:

.. math::
     \hat{q} = ((3, 6), (1.5, 1.5))

Normalizing yields the biased PDF:

.. math::
     \hat{p} = ((0.25, 0.5), (0.125, 0.125)

The weights of particle born from these phase space bins should then be the
ratio of the unbiased to biased PDF values:

.. math::
     w = ((1.6, 0.4), (2.4, 0.8))




**********
References
**********

[1] M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)

[2] A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)

[3] C. Rocchini and P. Cignoni, "Generating Random Points in a Tetrahedron," 
    Journal of Graphics Tools, 5, 200–202 (2001).

***************
Further Reading
***************

[4] E. Biondo, A. Davis, A. Scopatz, P. P.H. Wilson, "Rigorous Two-Step 
    Activation for Fusion Systems with PyNE,” Proc. of the 18th Topical 
    Meeting of the Radiation Protection \& Shielding Division of ANS, Knoxville, 
    TN (2014).

[5] Relson, E. "Improved Methods For Sampling Mesh-Based Volumetric Sources In 
    Monte Carlo Transport." MS thesis University of Wisconsin, Madison WI, 2013.
