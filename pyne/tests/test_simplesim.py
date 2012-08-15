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
        """Tests AxisCylinder's methods, properties, and exceptions."""

        ## __init__()
        cyl = cards.AxisCylinder('mycyl', 'z', 0.4)
        self.assertEquals(cyl.name, 'mycyl')
        self.assertEquals(cyl.cartesian_axis, 'z')
        self.assertEquals(cyl.radius, 0.4)
        self.assertIsNone(cyl.reflecting)
        self.assertIsNone(cyl.white)
        # Set reflecting.
        cyl = cards.AxisCylinder('mycyl', 'z', 0.4, reflecting=True)
        self.assertEquals(cyl.reflecting, True)
        self.assertIsNone(cyl.white)
        # Set white.
        cyl = cards.AxisCylinder('mycyl', 'z', 0.4, white=True)
        self.assertIsNone(cyl.reflecting)
        self.assertEquals(cyl.white, True)
        # Set reflecting and white.
        self.assertRaises(ValueError, cards.AxisCylinder, 'mycyl', 'z', 0.4,
                reflecting=True, white=True)
        # Test the exception message.
        try:
            cyl = cards.AxisCylinder('mycyl', 'z', 0.4, reflecting=True,
                    white=True)
        except ValueError as e:
            self.assertEquals(e.message, "The user set the surface to be "
                    "reflecting AND white, but can only be neither or "
                    "one of the two.")

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
                "with radius 0.4000 cm (diameter 0.8000 cm).")

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
        # Test the exception message.
        try:
            cyl.shift([3, 3, 0])
        except ValueError as e:
            self.assertEquals(e.message, "A cylinder aligned with the x axis "
                    "cannot be shifted in the y or z directions.")

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
        # Test the exception message.
        try:
            cyl.stretch([1, 3, 2])
        except ValueError as e:
            self.assertEquals(e.message, "Stretches perpendicular to the "
                    "axis must be uniform in the two perpendicular "
                    "directions. User provided y stretch 3.0000 and z "
                    "stretch 2.0000 for a x-aligned cylinder.")

    def test_AxisPlane(self):
        """Tests Plane's methods, properties, and exceptions."""

        ## __init__()
        plane = cards.AxisPlane('myplane', 'x', 3, reflecting=True) 
        self.assertEquals(plane.name, 'myplane')
        self.assertEquals(plane.cartesian_axis, 'x')
        self.assertEquals(plane.position, 3)
        # test_AxisCylinder() checks these booleans more thoroughly.
        self.assertEquals(plane.reflecting, True)
        self.assertIsNone(plane.white)

        ## comment()
        self.assertEquals(plane.comment(), "Axis plane myplane: x = 3.0000 cm")

        ## shift()
        plane = cards.AxisPlane('myplane', 'x', 3)
        plane.shift([3, 0, 0])
        self.assertEquals(plane.position, 6)
        plane.shift([0, 3, 2])
        self.assertEquals(plane.position, 6)
        plane = cards.AxisPlane('myplane', 'y', 3)
        plane.shift([0, 3, 0])
        self.assertEquals(plane.position, 6)
        plane.shift([3, 0, 2])
        self.assertEquals(plane.position, 6)
        plane = cards.AxisPlane('myplane', 'z', 3)
        plane.shift([0, 0, 3])
        self.assertEquals(plane.position, 6)
        plane.shift([2, 3, 0])
        self.assertEquals(plane.position, 6)

        ## stretch()
        plane = cards.AxisPlane('myplane', 'x', 3)
        plane.stretch([3, 0, 0])
        self.assertEquals(plane.position, 9)
        plane.stretch([0, 3, 2])
        self.assertEquals(plane.position, 9)
        plane = cards.AxisPlane('myplane', 'y', 3)
        plane.stretch([0, 3, 0])
        self.assertEquals(plane.position, 9)
        plane.stretch([3, 0, 2])
        self.assertEquals(plane.position, 9)
        plane = cards.AxisPlane('myplane', 'z', 3)
        plane.stretch([0, 0, 3])
        self.assertEquals(plane.position, 9)
        plane.stretch([2, 3, 0])
        self.assertEquals(plane.position, 9)

    def test_Parallelepiped(self):
        """Tests Parallelepiped's methods, properties, and exceptions."""

        ## __init__()
        pp = cards.Parallelepiped('mypp', -2, 3, -4, 5, -6, 7, reflecting=True,
                white=False)
        self.assertEquals(pp.name, 'mypp')
        self.assertTrue((pp.xlims == np.array([-2, 3])).all())
        self.assertTrue((pp.ylims == np.array([-4, 5])).all())
        self.assertTrue((pp.zlims == np.array([-6, 7])).all())
        self.assertEquals(pp.reflecting, True)
        self.assertEquals(pp.white, False)

        ## comment()
        self.assertEquals(pp.comment(), "Parallelepiped mypp: "
                "[-2.0000, 3.0000] x [-4.0000, 5.0000] x "
                "[-6.0000, 7.0000] cm")

        ## shift()
        pp.shift([2, 1, -1])
        self.assertTrue((pp.xlims == np.array([0, 5])).all())
        self.assertTrue((pp.ylims == np.array([-3, 6])).all())
        self.assertTrue((pp.zlims == np.array([-7, 6])).all())
            
        ## stretch()
        pp.stretch([4, 2, 3])
        self.assertTrue((pp.xlims == np.array([0, 20])).all())
        self.assertTrue((pp.ylims == np.array([-6, 12])).all())
        self.assertTrue((pp.zlims == np.array([-21, 18])).all())
        # Reflect.
        pp.stretch([-2, -1, -2])
        self.assertTrue((pp.xlims == np.array([-40, 0])).all())
        self.assertTrue((pp.ylims == np.array([-12, 6])).all())
        self.assertTrue((pp.zlims == np.array([-36, 42])).all())


