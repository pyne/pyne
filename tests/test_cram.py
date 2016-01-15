import nose
import numpy as np
from nose.plugins.skip import Skip, SkipTest
from numpy.testing import assert_almost_equal
from pyne import nucname, alara, data, material

def test_cram14():
    test_nuclides = ['H3', 'C14', 'V50', 'Rb87', 'In115', 'Sr90', 'Te123', 'Te130', 'I131', 'Cs137', 'La138', 'Nd144', 'Sm147', 'Lu176', 'Re187', 'Os186']
    error14 = []
    for species in test_nuclides:
        nuclides = [nucname.id(species)]

        children = data.decay_children(species)
        children = list(children)
        while children:
            nuclides += children
            parent = children[0]
            children = data.decay_children(parent)
            children = list(children)

        nuclides = sorted(nuclides)
        half_life = data.half_life(species) # Half life of species

        n_0 = np.zeros(len(nuclides))
        position = nuclides.index(nucname.id(species))
        n_0[position] = 1 # Set initial amount of material

        # Compute CRAM 14 Solution
        cram14 = alara.cram(nuclides, half_life, n_0, 14)

        # Compute PyNE decay solution
        initial = material.Material({species: 1.0})
        final = initial.decay(half_life)
        pdecay = list(final.values())
        error_14 = np.abs(cram14 - pdecay)*2 / (cram14 + pdecay)
        error14.append(np.mean(error_14))

    assert_almost_equal(np.mean(error14), 0, 2)

def test_cram16():
    test_nuclides = ['H3', 'C14', 'V50', 'Rb87', 'In115', 'Sr90', 'Te123', 'Te130', 'I131', 'Cs137', 'La138', 'Nd144', 'Sm147', 'Lu176', 'Re187', 'Os186']
    error16 = []
    for species in test_nuclides:
        nuclides = [nucname.id(species)]

        children = data.decay_children(species)
        children = list(children)
        while children:
            nuclides += children
            parent = children[0]
            children = data.decay_children(parent)
            children = list(children)

        nuclides = sorted(nuclides)
        half_life = data.half_life(species) # Half life of species

        n_0 = np.zeros(len(nuclides))
        position = nuclides.index(nucname.id(species))
        n_0[position] = 1 # Set initial amount of material

        # Compute CRAM 16 Solution
        cram16 = alara.cram(nuclides, half_life, n_0, 14)

        # Compute PyNE decay solution
        initial = material.Material({species: 1.0})
        final = initial.decay(half_life)
        pdecay = list(final.values())

        error_16 = np.abs(cram16 - pdecay)*2 / (cram16 + pdecay)
        error16.append(np.mean(error_16))

    assert_almost_equal(np.mean(error16), 0, 2)

