"""PyNE Simple Simulation Input Definitions tests"""

import os
import unittest

import numpy as np

from pyne.simplesim import definition, cards, inputfile

class TestSurfaces(unittest.TestCase):
    """Tests the ISurface class and its subclasses, all of which are subclassed
    from ICard and are in the :py:mod:`pyne.simplesim.cards` module."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_Cell(self):
        # test density_units exception.
        pass

    def test_AxisCylinder(self):
        """Checks AxisCylinder's methods, properties, and exceptions."""

        ## __init__()
        cyl = cards.AxisCylinder('mycyl', 'z', 0.4)
        self.assertEquals(cyl.name, 'mycyl')
        self.assertEquals(cyl.cartesian_axis, 'z')
        self.assertEquals(cyl.radius, 0.4)
        self.assertEquals(cyl.reflecting, None)
        self.assertEquals(cyl.white, None)
        # Z can be capitalized in input.
        cyl = cards.AxisCylinder('mycyl', 'Z', 0.4)
        self.assertEquals(cyl.cartesian_axis, 'z')
        cyl = cards.AxisCylinder('mycyl', 'Y', 0.4)
        self.assertEquals(cyl.cartesian_axis, 'y')
        cyl = cards.AxisCylinder('mycyl', 'X', 0.4)
        self.assertEquals(cyl.cartesian_axis, 'x')

        ## comment()
        self.assertEquals(cyl.comment(),
                "Axis cylinder mycyl: aligned and centered on x axis, "
                "with radius 0.4000 (diameter 0.8000).")

        ## shift()
        # z-aligned
        cyl = cards.AxisCylinder('mycyl', 'z', 0.4)
        cyl.shift([0, 0, 3])
        # Radius should not change.
        self.assertEquals(cyl.radius, 0.4)
        # Using a numpy array should work.
        cyl.shift(np.array([0, 0, 3]))
        # The axis cylinder cannot be shifted perpendicular to its axis.
        self.assertRaises(ValueError, cyl.shift, [3, 0, 0])
        self.assertRaises(ValueError, cyl.shift, [0, 3, 3])
        # y-aligned
        cyl = cards.AxisCylinder('mycyl', 'y', 0.4)
        cyl.shift([0, 3, 0])
        # Radius should not change.
        self.assertEquals(cyl.radius, 0.4)
        # Using a numpy array should work.
        cyl.shift(np.array([0, 3, 0]))
        # The axis cylinder cannot be shifted perpendicular to its axis.
        self.assertRaises(ValueError, cyl.shift, [3, 0, 0])
        self.assertRaises(ValueError, cyl.shift, [0, 3, 3])
        # x-aligned
        cyl = cards.AxisCylinder('mycyl', 'x', 0.4)
        cyl.shift([3, 0, 0])
        # Radius should not change.
        self.assertEquals(cyl.radius, 0.4)
        # Using a numpy array should work.
        cyl.shift(np.array([3, 0, 0]))
        # The axis cylinder cannot be shifted perpendicular to its axis.
        self.assertRaises(ValueError, cyl.shift, [0, 0, 3])
        self.assertRaises(ValueError, cyl.shift, [3, 3, 0])

        ## stretch()
        # z-aligned
        cyl = cards.AxisCylinder('mycyl', 'z', 0.4)
        # Stretch along axis has no effect.
        cyl.stretch([0, 0, 2])
        self.assertEquals(cyl.radius, 0.4)
        # Uniform perpendicular stretch.
        cyl.stretch([3, 3, 0])
        self.assertAlmostEqual(cyl.radius, 1.2)
        # Combined stretch.
        cyl.stretch([3, 3, 2])
        self.assertAlmostEqual(cyl.radius, 3.6)
        # Non-uniform perpendicular stretches.
        self.assertRaises(ValueError, cyl.stretch, [0, 3, 0])
        self.assertRaises(ValueError, cyl.stretch, [2, 3, 1])
        # y-aligned
        cyl = cards.AxisCylinder('mycyl', 'y', 0.4)
        cyl.stretch([0, 2, 0])
        self.assertEquals(cyl.radius, 0.4)
        cyl.stretch([3, 0, 3])
        self.assertAlmostEqual(cyl.radius, 1.2)
        cyl.stretch([3, 2, 3])
        self.assertAlmostEqual(cyl.radius, 3.6)
        self.assertRaises(ValueError, cyl.stretch, [0, 0, 3])
        self.assertRaises(ValueError, cyl.stretch, [2, 1, 3])
        # x-aligned
        cyl = cards.AxisCylinder('mycyl', 'x', 0.4)
        cyl.stretch([2, 0, 0])
        self.assertEquals(cyl.radius, 0.4)
        cyl.stretch([0, 3, 3])
        self.assertAlmostEqual(cyl.radius, 1.2)
        cyl.stretch([2, 3, 3])
        self.assertAlmostEqual(cyl.radius, 3.6)
        self.assertRaises(ValueError, cyl.stretch, [0, 3, 0])
        self.assertRaises(ValueError, cyl.stretch, [1, 3, 2])



