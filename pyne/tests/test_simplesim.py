"""PyNE Simple Simulation Input Definitions tests"""

import os
import unittest

import numpy as np

from pyne import material
from pyne.simplesim import definition, cards, inputfile

# TODO test the exception that setting a name to '' is not aight.


class TestCells(unittest.TestCase):
    """Tests the :py:class:`cards.Cell` class and its subclasses, all of which
    are subclassed from :py:class:`cards.ICard` and are in the :py:mod:`cards`
    module.

    """
        # test density_units exception.
        # test '-inf' test.
        # test dxtran sphere reference.

class TestSurfaces(unittest.TestCase):
    """Tests the :py:class:`cards.ISurface` class and its
    subclasses, all of which are subclassed from :py:class:`cards.ICard` and
    are in the :py:mod:`cards` module.
    
    """
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_AxisCylinder(self):
        """Tests :py:class:`cards.AxisCylinder`'s methods, properties, and
        exceptions.
        
        """
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
                "Axis cylinder 'mycyl': aligned and centered on x axis, "
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
        """Tests :py:class:`cards.AxisPlane`'s methods, properties, and
        exceptions.
        
        """
        ## __init__()
        plane = cards.AxisPlane('myplane', 'x', 3, reflecting=True) 
        self.assertEquals(plane.name, 'myplane')
        self.assertEquals(plane.cartesian_axis, 'x')
        self.assertEquals(plane.position, 3)
        # test_AxisCylinder() checks these booleans more thoroughly.
        self.assertEquals(plane.reflecting, True)
        self.assertIsNone(plane.white)

        ## comment()
        self.assertEquals(plane.comment(), "Axis plane 'myplane': x = 3.0000 cm.")

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
        """Tests :py:class:`cards.Parallelepiped`'s methods, properties, and
        exceptions.
        
        """
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
        self.assertEquals(pp.comment(), "Parallelepiped 'mypp': "
                "[-2.0000, 3.0000] x [-4.0000, 5.0000] x "
                "[-6.0000, 7.0000] cm.")

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
    

class TestOptions(unittest.TestCase):
    """Tests the :py:class:`cards.IOptions` class and its subclasses, all of which are subclassed
    from :py:class:`cards.ICard` and are in the :py:mod:`cards` module.
    
    """
    def tests_Criticality(self):
        """Tests :py:class:`cards.Criticality`'s methods, properties, and
        exceptions.
        
        """
        ## __init__()
        # Default arguments.
        critsrc = cards.Criticality()
        self.assertEquals(critsrc.name, 'criticality')
        # Test trying to change the name.
        self.assertRaises(StandardError, setattr, critsrc, 'name', 'testname')
        try:
            critsrc.name = 'testname'
        except StandardError as e:
            self.assertEquals(e.message, "This is a unique card, meaning "
                    "only one card of this type can be found in a "
                    "``definition``. Accordingly, the name is read-only.")
        self.assertEquals(critsrc.n_histories, 1000)
        self.assertEquals(critsrc.keff_guess, 1)
        self.assertEquals(critsrc.n_skip_cycles, 30)
        self.assertEquals(critsrc.n_cycles, 130)
        # Testing keyword argument functionality.
        critsrc = cards.Criticality(n_cycles=300)
        self.assertEquals(critsrc.n_histories, 1000)
        self.assertEquals(critsrc.keff_guess, 1)
        self.assertEquals(critsrc.n_skip_cycles, 30)
        self.assertEquals(critsrc.n_cycles, 300)
        # Testing documentation example.
        critsrc = cards.Criticality(2000, 1.5, 30, 300)
        self.assertEquals(critsrc.n_histories, 2000)
        self.assertEquals(critsrc.keff_guess, 1.5)
        self.assertEquals(critsrc.n_skip_cycles, 30)
        self.assertEquals(critsrc.n_cycles, 300)

        # Test exceptions for each of the properties.
        self.assertRaises(ValueError, setattr, critsrc, 'n_histories', 0.5)
        self.assertRaises(ValueError, setattr, critsrc, 'n_histories', -1)
        self.assertRaises(ValueError, setattr, critsrc, 'keff_guess', -1)
        self.assertRaises(ValueError, setattr, critsrc, 'n_skip_cycles', 0.5)
        self.assertRaises(ValueError, setattr, critsrc, 'n_skip_cycles', -1)
        self.assertRaises(ValueError, setattr, critsrc, 'n_cycles', 0.5)
        self.assertRaises(ValueError, setattr, critsrc, 'n_cycles', -1)
        # Test the exception messages.
        try:
            critsrc.n_histories = 0.5
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_histories`` "
                    "must be an integer. User provided 0.5000.")
        try:
            critsrc.n_histories = -1
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_histories`` "
                    "must be positive. User provided -1.")
        try:
            critsrc.keff_guess = -1
        except ValueError as e:
            self.assertEquals(e.message, "The property ``keff_guess`` "
                    "must be non-negative. User provided -1.0000.")
        try:
            critsrc.n_skip_cycles = 0.5
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_skip_cycles`` "
                    "must be an integer. User provided 0.5000.")
        try:
            critsrc.n_skip_cycles = -1
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_skip_cycles`` "
                    "must be positive. User provided -1.")
        try:
            critsrc.n_cycles = 0.5
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_cycles`` "
                    "must be an integer. User provided 0.5000.")
        try:
            critsrc.n_cycles = -1
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_cycles`` "
                    "must be equal to or greater than ``n_skip_cycles``. "
                    "User provided -1.")

        ## comment()
        critsrc = cards.Criticality()
        self.assertEquals(critsrc.comment(), "Criticality source "
                "'criticality': n_histories: 1000, keff_guess: 1.0000"
                ", n_skip_cycles: 30, n_cycles: 130.")
            
    def test_CriticalityPoints(self):
        """Tests :py:class:`cards.CriticalityPoints`'s methods, properties, and
        exceptions.
        
        """
        ## __init__()
        # Default arguments.
        critpts = cards.CriticalityPoints()
        self.assertEquals(critpts.name, 'criticalitypoints')
        self.assertTrue(critpts.points == [[0, 0, 0]])
        critpts = cards.CriticalityPoints([[1, 2, 3],
                                        np.array([np.pi, np.e, 0])])
        self.assertEquals(critpts.comment(), "Criticality points "
                "'criticalitypoints': (1.0000, 2.0000, 3.0000), "
                "(3.1416, 2.7183, 0.0000).")
        self.assertRaises(ValueError, setattr, critpts, 'points', 
                [[0, 0, 0], [1]])
        try:
            critpts.points = [[0, 0, 0], [1]]
        except ValueError as e:
            self.assertEquals(e.message, "Length of all point lists/arrays "
                    "must be 3. User provided a point [1].")


class TestSystemDefinition(unittest.TestCase):
    """Tests the :py:class:`definition.SystemDefinition` class."""

    def setUp(self):
        self.uo2 = material.from_atom_frac({'U235': 0.05,
                            'U238': 0.95,
                            'O16' : 2.0}, name='UO2')
        self.h2o = material.from_atom_frac({'H1' : 2.0,
                            'O16': 1.0}, name='H2O')
        
        # Surfaces.
        radius = 0.40 # cm
        self.pin = cards.AxisCylinder('fuelpin', 'X', radius)
        pitch = 1.2 # cm
        self.cellbound = cards.Parallelepiped('bound',
                -pitch / 2, pitch / 2, -pitch / 2, pitch / 2, 0, 0,
                reflecting=True)
        
        # Cells.
        self.fuel = cards.CellMCNP('fuel', self.pin.neg, self.uo2, 11.0,
                'g/cm^3', neutron_imp=1)
        self.coolant = cards.CellMCNP('coolant', self.pin.pos &
                self.cellbound.neg, self.h2o, 1.0, 'g/cm^3', neutron_imp=1)
        self.graveyard = cards.CellVoidMCNP('graveyard', self.cellbound.pos,
                neutron_imp=0)
        
        # Create system definition from the cards above.
        self.rxr = definition.SystemDefinition(verbose=False)
        self.rxr.add_cell(self.fuel)
        self.rxr.add_cell(self.coolant)
        self.rxr.add_cell(self.graveyard)

    def test_ITally(self):
        """Tests :py:class:`cards.ITally`'s methods, properties, and
        exceptions, and those of its subclasses.

        """
        ## SurfaceCurrent
        tally = cards.SurfaceCurrent('fuel', 'electron', [self.pin,
                self.cellbound], total=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle, 'electron')
        self.assertIs(tally.cards[0], self.pin)
        self.assertIs(tally.cards[1], self.cellbound)
        self.assertTrue(len(tally.cards), 2)
        self.assertFalse(hasattr(tally, 'average'))
        self.assertTrue(tally.total)
        self.assertFalse(tally.alt_units)
        # comment()
        self.assertEquals(tally.comment(), "Surface current tally 'fuel' of "
                "electrons: surfaces 'fuelpin'; 'bound'; and total of all "
                "provided.")
        tally = cards.SurfaceCurrent('fuel', 'photon', [[self.pin,
                self.cellbound]], alt_units=True)
        self.assertEquals(tally.particle, 'photon')
        self.assertTrue(len(tally.cards), 1)
        self.assertTrue(len(tally.cards[0]), 2)
        self.assertIs(tally.cards[0][0], self.pin)
        self.assertIs(tally.cards[0][1], self.cellbound)
        self.assertFalse(tally.total)
        self.assertTrue(tally.alt_units)
        # comment()
        self.assertEquals(tally.comment(), "Surface current tally 'fuel' "
                "of photons: total in 'fuelpin', 'bound'.")

        ## SurfaceFlux
        tally = cards.SurfaceFlux('fuel', 'electron', [self.pin,
                self.cellbound], average=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle, 'electron')
        self.assertIs(tally.cards[0], self.pin)
        self.assertIs(tally.cards[1], self.cellbound)
        self.assertTrue(len(tally.cards), 2)
        self.assertTrue(tally.average)
        self.assertFalse(tally.alt_units)
        # comment()
        self.assertEquals(tally.comment(), "Surface flux tally 'fuel' "
                "of electrons: surfaces 'fuelpin'; 'bound'; and avg. "
                "of all provided.")
        tally = cards.SurfaceFlux('fuel', 'proton', [[self.pin, self.cellbound]],
                    alt_units=True)
        self.assertEquals(tally.particle, 'proton')
        self.assertTrue(len(tally.cards), 1)
        self.assertTrue(len(tally.cards[0]), 2)
        self.assertIs(tally.cards[0][0], self.pin)
        self.assertIs(tally.cards[0][1], self.cellbound)
        self.assertFalse(tally.average)
        self.assertTrue(tally.alt_units)
        # comment()
        self.assertEquals(tally.comment(), "Surface flux tally 'fuel' "
                "of protons: avg. in 'fuelpin', 'bound'.")

        ## CellFlux
        # One cell.
        tally = cards.CellFlux('fuel', 'neutron', self.fuel)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle, 'neutron')
        # Test ``particle`` property.
        tally = cards.CellFlux('fuel', 'photon', self.fuel)
        self.assertEquals(tally.particle, 'photon')
        tally = cards.CellFlux('fuel', 'electron', self.fuel)
        self.assertEquals(tally.particle, 'electron')
        tally = cards.CellFlux('fuel', 'proton', self.fuel)
        self.assertEquals(tally.particle, 'proton')
        # Test exception on ``particle`` property.
        self.assertRaises(ValueError, setattr, tally, 'particle', 'test')
        try:
            tally.particle = 'test'
        except ValueError as e:
            self.assertEquals(e.message, 
                    "The property ``particle`` must be 'neutron', "
                    "'photon', 'electron', or 'proton'. "
                    "User provided 'test'.")
        self.assertEquals(tally.comment(), "Cell flux tally 'fuel' "
                "of protons: cell 'fuel'.")
        self.assertIs(tally._unique_card_list()[0], self.fuel)
        # Two individual cells.
        tally = cards.CellFlux('both', 'neutron', [self.fuel, self.coolant])
        self.assertEquals(tally.comment(), "Cell flux tally 'both' "
                "of neutrons: cells 'fuel'; 'coolant'.")
        self.assertIs(tally._unique_card_list()[0], self.fuel)
        self.assertIs(tally._unique_card_list()[1], self.coolant)
        # Two individual cells, with average over all.
        tally = cards.CellFlux('withavg', 'neutron', [self.fuel, self.coolant],
                average=True)
        self.assertEquals(tally.comment(), "Cell flux tally 'withavg' "
                "of neutrons: cells 'fuel'; 'coolant'; and avg. of all "
                "provided.")
        self.assertIs(tally._unique_card_list()[0], self.fuel)
        self.assertIs(tally._unique_card_list()[1], self.coolant)
        # Two individual cells, and an averaging, with an average over all.
        tally = cards.CellFlux('withavg', 'neutron', [self.fuel,
                [self.fuel, self.coolant], self.coolant], average=True)
        self.assertEquals(tally.comment(), "Cell flux tally 'withavg' "
                "of neutrons: cells 'fuel'; avg. in 'fuel', 'coolant'; "
                "'coolant'; and avg. of all provided.")
        self.assertIs(tally._unique_card_list()[0], self.fuel)
        self.assertIs(tally._unique_card_list()[1], self.coolant)
        self.assertTrue(len(tally._unique_card_list()) == 2)

        ## CellEnergyDeposition
        tally = cards.CellEnergyDeposition('energy', 'neutron', self.fuel)
        self.assertEquals(tally.name, 'energy')
        self.assertEquals(tally.particle, 'neutron')
        self.assertIs(tally.cards, self.fuel)
        self.assertEquals(tally.comment(), "Energy deposition tally "
                "'energy' of neutrons: cell 'fuel'.")
        tally = cards.CellEnergyDeposition('energy', ['neutron', 'proton'],
                self.fuel)
        self.assertIs(type(tally.particle), list)
        self.assertEquals(tally.particle[0], 'neutron')
        self.assertEquals(tally.particle[1], 'proton')
        self.assertEquals(tally.comment(), "Energy deposition tally " 
                "'energy' of neutrons, protons: cell 'fuel'.")
        tally = cards.CellEnergyDeposition('energy', 'all', self.fuel)
        self.assertEquals(tally.comment(), "Energy deposition tally "
                "'energy' of all: cell 'fuel'.")
        # Test exceptions.
        self.assertRaises(ValueError, cards.CellEnergyDeposition, 'energy',
                ['neutron', 'all'], self.fuel)
        self.assertRaises(ValueError, cards.CellEnergyDeposition,
                'energy', 'all', self.fuel, alt_units=True)

        ## CellFissionEnergyDeposition
        tally = cards.CellFissionEnergyDeposition('fuel', [self.fuel,
                self.coolant], average=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle, 'neutron')
        self.assertIs(type(tally.cards), list)
        self.assertIs(tally.cards[0], self.fuel)
        self.assertIs(tally.cards[1], self.coolant)
        self.assertTrue(tally.average)
        self.assertFalse(tally.alt_units)
        self.assertEquals(tally.comment(), "Fission energy deposition tally "
                "'fuel' of neutrons: cells 'fuel'; 'coolant'; and avg. of "
                "all provided.")
        tally = cards.CellFissionEnergyDeposition('fuel', [[self.fuel,
                self.coolant]], alt_units=True)
        self.assertTrue(len(tally.cards), 1)
        self.assertTrue(len(tally.cards[0]), 2)
        self.assertIs(tally.cards[0][0], self.fuel)
        self.assertIs(tally.cards[0][1], self.coolant)
        self.assertFalse(tally.average)
        self.assertTrue(tally.alt_units)

        ## CellPulseHeight
        tally = cards.CellPulseHeight('fuel', ['proton', 'electron'], [self.fuel,
                self.coolant], alt_units=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle[0], 'proton')
        self.assertEquals(tally.particle[1], 'electron')
        self.assertIs(tally.cards[0], self.fuel)
        self.assertIs(tally.cards[1], self.coolant)
        self.assertFalse(tally.average)
        self.assertTrue(tally.alt_units)
        tally.average = True
        self.assertTrue(tally.average)

        ## CellChargeDeposition
        tally = cards.CellChargeDeposition('fuel', ['proton', 'electron'],
                [self.fuel, self.coolant])
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle[0], 'proton')
        self.assertEquals(tally.particle[1], 'electron')
        self.assertIs(tally.cards[0], self.fuel)
        self.assertIs(tally.cards[1], self.coolant)
        self.assertFalse(tally.average)
        self.assertFalse(tally.alt_units)
        tally.average = True
        self.assertTrue(tally.average)

        ## PointDetector
        det = cards.PointDetector('point', 'neutron', ([0, 0, 0], 0))
        self.assertEquals(det.name, 'point')
        self.assertEquals(det.particle, 'neutron')
        self.assertEquals(det.spec[0], [0, 0, 0])
        self.assertEquals(det.spec[1], 0)
        self.assertEquals(det.sep_direct, True)
        self.assertEquals(det.comment(), "Point detector tally 'point' of "
                "neutrons: point (0.0000, 0.0000, 0.0000) cm, "
                "radius 0.0000 cm; direct contrib is separate.")
        det = cards.PointDetector('point', 'neutron', (np.array([1, 2, 3]), 4))
        self.assertTrue((det.spec[0] == [1, 2, 3]).all())
        self.assertEquals(det.spec[1], 4)
        self.assertEquals(det.sep_direct, True)
        self.assertEquals(det.comment(), "Point detector tally 'point' of "
                "neutrons: point (1.0000, 2.0000, 3.0000) cm, "
                "radius 4.0000 cm; direct contrib is separate.")
        det = cards.PointDetector('point', 'neutron', ([1, 0, 0], -3))
        self.assertEquals(det.comment(), "Point detector tally 'point' of "
                "neutrons: point (1.0000, 0.0000, 0.0000) cm, "
                "radius 3.0000 mfp; direct contrib is separate.")
        det = cards.PointDetector('point', 'photon', [([0, 0, 0],  0),
                                                      ([1, 0, 0], -3)])
        self.assertEquals(det.spec[0][0], [0, 0, 0])
        self.assertEquals(det.spec[0][1], 0)
        self.assertEquals(det.spec[1][0], [1, 0, 0])
        self.assertEquals(det.spec[1][1], -3)
        self.assertEquals(det.comment(), "Point detector tally 'point' of "
                "photons: point (0.0000, 0.0000, 0.0000) cm, "
                "radius 0.0000 cm; "
                "point (1.0000, 0.0000, 0.0000) cm, radius 3.0000 mfp; "
                "direct contrib is separate.")
        det = cards.PointDetector('point', 'photon', ([0, 0, 0], 0),
                sep_direct=False)
        self.assertFalse(det.sep_direct)

        ## RingDetector
        det = cards.RingDetector('ring', 'neutron', ('x', 10.0, 2.0,  1.0))
        self.assertEquals(det.name, 'ring')
        self.assertEquals(det.particle, 'neutron')
        self.assertTrue(det.sep_direct)
        self.assertEquals(det.comment(), "Ring detector tally 'ring' of "
                "neutrons: ring x = 10.0000 cm, radius 2.0000 cm, s.o.e. "
                "radius 1.0000 cm; direct contrib is separate.")
        det = cards.RingDetector('ring', 'neutron', ('x', 10.0, 2.0, -1.0))
        self.assertEquals(det.comment(), "Ring detector tally 'ring' of "
                "neutrons: ring x = 10.0000 cm, radius 2.0000 cm, s.o.e. "
                "radius 1.0000 mfp; direct contrib is separate.")
        det = cards.RingDetector('ring', 'neutron', [('x', 10.0, 2.0, -1.0),
                                                     ('y', 20.0, 3.0, 1.0)])
        self.assertEquals(det.comment(), "Ring detector tally 'ring' of "
                "neutrons: ring x = 10.0000 cm, radius 2.0000 cm, s.o.e. "
                "radius 1.0000 mfp; ring y = 20.0000 cm, radius 3.0000 "
                "cm, s.o.e. radius 1.0000 cm; direct contrib is separate.")
        det = cards.RingDetector('ring', 'neutron', ('x', 10.0, 2.0, -1.0), 
                sep_direct=False)
        self.assertEquals(det.comment(), "Ring detector tally 'ring' of "
                "neutrons: ring x = 10.0000 cm, radius 2.0000 cm, s.o.e. "
                "radius 1.0000 mfp; direct contrib is not separate.")

        ## EnergyGrid
        egrid = cards.EnergyGrid('grid0', None, [1e-4, 1, 100e3, 10e6])
        self.assertEquals(egrid.comment(), "Energy grid 'grid0' for "
                "all tallies: 4 groups.")
        egrid = cards.EnergyGrid('grid0', None, np.array([1e-4, 1, 100e3,
                    10e6]))
        self.assertEquals(egrid.comment(), "Energy grid 'grid0' for "
                "all tallies: 4 groups.")


class TestSimulationDefinition(unittest.TestCase):
    """Tests the :py:class:`definition.SimulationDefinition` class."""
    # The system definition is complete.
    
    #sim = definition.MCNPSimulation(rxr)
    #sim.add_card(cards.Criticality())
    #sim.add_card(cards.CriticalityPoints())
