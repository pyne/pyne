"""PyNE Simple Simulation Input Definitions tests"""

import os
import unittest

import numpy as np

from pyne import material
from pyne.simplesim import definition, cards, inputfile

# TODO test the exception that setting a name to '' is not aight.
# TODO use assertRaisesRegexp

class TestCells(unittest.TestCase):
    """Tests the :py:class:`cards.Cell` class and its subclasses, all of which
    are subclassed from :py:class:`cards.ICard` and are in the :py:mod:`cards`
    module.

    """
        # test density_units exception.
        # test '-inf' test.
        # test dxtran sphere reference.
        # test mcnp particle designators.
    pass

class TestSurfaces(unittest.TestCase):
    """Tests the :py:class:`cards.ISurface` class and its
    subclasses, all of which are subclassed from :py:class:`cards.ICard` and
    are in the :py:mod:`cards` module.
    
    """
    def setUp(self):
        pass

    def tearDown(self):
        pass
    

class TestOptions(unittest.TestCase):
    """Tests the :py:class:`cards.IOptions` class and its subclasses, all of which are subclassed
    from :py:class:`cards.ICard` and are in the :py:mod:`cards` module.
    
    """


class TestSystemDefinition(unittest.TestCase):
    """Tests the :py:class:`definition.SystemDefinition` class."""

    def setUp(self):
        self.uo2 = cards.Material(name='UO2')
        self.uo2.from_atom_frac({'U235': 0.05,
                                 'U238': 0.95,
                                 'O16' : 2.0})
        self.h2o = cards.Material(name='H2O')
        self.h2o.from_atom_frac({'H1' : 2.0,
                                 'O16': 1.0})
        
        # Surfaces.
        radius = 0.40 # cm
        self.pin = cards.AxisCylinder('fuelpin', 'X', radius)
        pitch = 1.2 # cm
        self.cellbound = cards.Parallelepiped('bound',
                -pitch / 2, pitch / 2, -pitch / 2, pitch / 2, 0, 0,
                reflecting=True)
        
        # Cells.
        self.fuel = cards.CellMCNP('fuel', self.pin.neg, self.uo2, 11.0,
                'g/cm^3', importance=('neutron', 1))
        self.coolant = cards.CellMCNP('coolant', self.pin.pos &
                self.cellbound.neg, self.h2o, 1.0, 'g/cm^3', 
                importance=('neutron', 1))
        self.graveyard = cards.CellMCNP('graveyard', self.cellbound.pos,
                importance=('neutron', 0))
        
        # Create system definition from the cards above.
        self.rxr = definition.SystemDefinition(verbose=False)
        self.rxr.add_cell(self.fuel)
        self.rxr.add_cell(self.coolant)
        self.rxr.add_cell(self.graveyard)
        self.sim = definition.MCNPSimulation(self.rxr, verbose=False)

    def test_Material(self):
        """Tests :py:class:`pyne.simplesim.cards.Material`."""

        originstory = "I found this water in a well a few years ago."
        h2o = cards.Material(name='water', description=originstory)
        h2o.from_atom_frac({10010: 1.0, 'O16': 2.0})
        h2o.tables = {10010: '71c'}
        self.sim.sys.add_material(h2o)
        self.assertEquals(h2o.comment(), "Material 'water': "
                "I found this water in a well a few years ago.")
        self.assertEquals(h2o.mcnp('%g', self.sim), "M3\n"
        "       1001.71c  1 $ H1\n"
        "       8016      2 $ O16")

        h2o = cards.Material(name='water', tables={10010: '71c'})
        h2o.from_atom_frac({10010: 1.0, 'O16': 2.0})
        self.assertEquals(h2o.mcnp('%g', self.sim), "M3\n"
        "       1001.71c  1 $ H1\n"
        "       8016      2 $ O16")

        h2o = cards.Material(name='water', tables={'H1': '71c'})
        h2o.from_atom_frac({10010: 1.0, 'O16': 2.0})
        self.assertEquals(h2o.mcnp('%g', self.sim), "M3\n"
        "       1001.71c  1 $ H1\n"
        "       8016      2 $ O16")

        ## ScatteringLaw
        sl = cards.ScatteringLaw('water', {'H1': 'lwtr.16t', 80160: 'madeup'})
        self.assertEquals(sl.comment(), "Scattering law 'scatlaw-water': "
                "O16: madeup, H1: lwtr.16t.")
        self.assertEquals(sl.mcnp('%g', self.sim), "MT3 madeup lwtr.16t")

    def test_AxisCylinder(self):
        """Tests :py:class:`cards.AxisCylinder`'s methods, properties, and
        exceptions.
        
        """
        ## __init__()
        cyl = cards.AxisCylinder('mycyl', 'z', 0.4)
        self.assertEquals(cyl.name, 'mycyl')
        self.assertEquals(cyl.cartesian_axis, 'z')
        self.assertEquals(cyl.radius, 0.4)
        self.assertFalse(cyl.reflecting)
        self.assertFalse(cyl.white)
        # Set reflecting.
        cyl = cards.AxisCylinder('mycyl', 'z', 0.4, reflecting=True)
        self.assertEquals(cyl.reflecting, True)
        self.assertFalse(cyl.white)
        # Set white.
        cyl = cards.AxisCylinder('mycyl', 'z', 0.4, white=True)
        self.assertFalse(cyl.reflecting)
        self.assertTrue(cyl.white)
        ## comment()
        self.assertEquals(cyl.comment(),
                "Axis cylinder 'mycyl': white. "
                "aligned and centered on z axis, "
                "with radius 0.4 cm (diameter 0.8 cm).")

        ## mcnp()
        self.sim.sys.add_surface(cyl)
        self.assertEquals(cyl.mcnp('%g', self.sim), "+3 CZ  0.4")

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
                "with radius 0.4 cm (diameter 0.8 cm).")

        ## mcnp()
        self.sim.sys.add_surface(cyl)
        self.assertEquals(cyl.mcnp('%g', self.sim), "3  CX  0.4")

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
                    "directions. User provided y stretch 3 and z "
                    "stretch 2 for a x-aligned cylinder.")

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
        self.assertFalse(plane.white)

        ## comment()
        self.assertEquals(plane.comment(), "Axis plane 'myplane': "
                "reflecting. x = 3 cm.")

        ## mcnp()
        self.sim.sys.add_surface(plane)
        self.assertEquals(plane.mcnp('%g', self.sim), "*3 PX  3")

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
                "reflecting. "
                "[-2, 3] x [-4, 5] x "
                "[-6, 7] cm.")

        ## mcnp()
        self.sim.sys.add_surface(pp)
        self.assertEquals(pp.mcnp('% g', self.sim), 
                "*3 RPP -2  3  -4  5  -6  7")

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

    def test_Cell(self):
        """Tests :py:class:`cards.Cell`'s methods, properties, and
        exceptions.

        """

        cellA = cards.Cell('A', self.pin.neg & self.cellbound.pos)
        self.rxr.add_cell(cellA)
        self.assertEquals(cellA.comment(), "Cell 'A': region "
                "(-fuelpin & +bound), void.")
        self.assertEquals(cellA.mcnp("%.5e", self.sim), "4 0 (-1 2)")
        cellB = cards.Cell('B', self.pin.neg & self.cellbound.pos, self.uo2,
                10.0, 'g/cm^3')
        self.rxr.add_cell(cellB)
        self.assertEquals(cellB.comment(), "Cell 'B': region "
                "(-fuelpin & +bound), material 'UO2' density 10 "
                "g/cm^3.")
        self.assertEquals(cellB.mcnp("%.5e", self.sim), "5 1 -1.00000e+01 "
                "(-1 2)")
        # Incorrect density string.
        self.assertRaises(ValueError, cards.Cell, 'B', self.pin.neg &
                self.cellbound.pos, self.uo2, 10.0, 'g/cm3')
        # self.fuel is not a Material.
        self.assertRaises(ValueError, cards.Cell, 'B', self.pin.neg &
            self.cellbound.pos, self.fuel, 10.0, 'g/cm^3')
        cellC = cards.Cell('C', self.pin.neg & self.cellbound.pos, self.uo2,
                0.5, 'atoms/b/cm')
        self.assertRaises(ValueError, cards.Cell, 'B', self.pin.neg &
            self.cellbound.pos, self.uo2)
        self.assertRaises(ValueError, cards.Cell, 'B', self.pin.neg &
            self.cellbound.pos, self.uo2, density=1)
        self.assertRaises(ValueError, cards.Cell, 'B', self.pin.neg &
            self.cellbound.pos, self.uo2, density_units='g/cm^3')
        self.assertRaises(ValueError, cards.Cell, 'B', self.pin.neg &
            self.cellbound.pos, density=1)
        self.assertRaises(ValueError, cards.Cell, 'B', self.pin.neg &
            self.cellbound.pos, density_units='g/cm^3')

    def test_CellMCNP(self):
        """Tests :py:class:`cards.CellMCNP`'s methods, properties, and
        exceptions. Tests each kwarg, and test a
        few in combination.

        """
        # Without keyword arguments
        cellC = cards.Cell('C', self.pin.neg & self.cellbound.pos)
        self.rxr.add_cell(cellC)
        self.assertEquals(cellC.comment(), "Cell 'C': region "
                "(-fuelpin & +bound), void.")
        self.assertEquals(cellC.mcnp("%.5e", self.sim), "4 0 (-1 2)")
        cellD = cards.Cell('D', self.pin.neg & self.cellbound.pos, self.uo2,
                10.0, 'g/cm^3')
        self.rxr.add_cell(cellD)
        self.assertEquals(cellD.comment(), "Cell 'D': region "
                "(-fuelpin & +bound), material 'UO2' density 10 "
                "g/cm^3.")
        self.assertEquals(cellD.mcnp("%.5e", self.sim), "5 1 -1.00000e+01 "
                "(-1 2)")

        # Keywords.
        # Temperature, volume, importance.
        cellE = cards.CellMCNP('E', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                temperature=600, volume=1,
                importance=('neutron', 1))
        self.rxr.add_cell(cellE)
        self.assertEquals(cellE.comment(), "Cell 'E': region "
                "(-fuelpin & +bound), material 'UO2' density 10 "
                "g/cm^3 TMP= 600 K VOL= 1 cm^3 IMP:N= 1.")
        self.assertEquals(cellE.mcnp('%g', self.sim), "6 1 -10 "
                "(-1 2) TMP=5.17041e-08 VOL=1 IMP:N=1")
        cellF = cards.CellMCNP('F', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                importance=[('neutron', 1), ('photon', 0)])
        self.rxr.add_cell(cellF)
        self.assertEquals(cellF.comment(), "Cell 'F': region "
                "(-fuelpin & +bound), material 'UO2' density 10 "
                "g/cm^3 IMP:N= 1 IMP:P= 0.")
        self.assertEquals(cellF.mcnp('%g', self.sim), "7 1 -10 "
                "(-1 2) IMP:N=1 IMP:P=0")
        # Exponential transform
        cellG = cards.CellMCNP('G', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                exp_transform=[
                ('neutron', 'capture-to-total', 'currdir', 'toward'),
                ('photon', 0.5, 'x', 'away')])
        self.rxr.add_cell(cellG)
        self.assertEquals(cellG.comment(), "Cell 'G': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "EXT:N= stretch by capture-to-total toward currdir "
                "EXT:P= stretch by 0.5 away from x.")
        self.assertEquals(cellG.mcnp('%g', self.sim), "8 1 -10 "
                "(-1 2) EXT:N=S EXT:P=-0.5X")
        cellH = cards.CellMCNP('H', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                exp_transform=('neutron', 0.5, 'vec1', 'away'))
        vec = cards.Vector()
        vec.set('vec1', [0, 0, 0])
        self.rxr.add_cell(cellH)
        self.sim.add_misc(vec)
        self.assertEquals(cellH.comment(), "Cell 'H': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "EXT:N= stretch by 0.5 away from vec1.")
        self.assertEquals(cellH.mcnp('%g', self.sim), "9 1 -10 "
                "(-1 2) EXT:N=-0.5V0")
        # Forced collisions
        cellI = cards.CellMCNP('I', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                forced_coll=('neutron', 0.5, False))
        self.rxr.add_cell(cellI)
        self.assertEquals(cellI.comment(), "Cell 'I': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "FCL:N= prob 0.5 for entering and weight games.")
        self.assertEquals(cellI.mcnp('%g', self.sim), "10 1 -10 "
                "(-1 2) FCL:N=0.5")
        # Weight window bound
        # Must add weight window energies card.
        cellJ = cards.CellMCNP('J', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                weight_win_bound=('neutron', 3, 1, 'killall'))
        wwe = cards.WeightWindowEnergies('neutron', [1, 10, 12])
        self.rxr.add_cell(cellJ)
        self.sim.add_misc(wwe)
        self.assertEquals(cellJ.comment(), "Cell 'J': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "WWN(3,1):N= killall.")
        self.assertEquals(cellJ.mcnp('%g', self.sim), "11 1 -10 "
                "(-1 2) WWN3:N=-1")
        # DXTRAN contribution
        cellK = cards.CellMCNP('K', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                dxtran_contrib=('neutron', 'sph1', 0.5))
        dsph = cards.DXTRANSpheres('neutron', 'sph1', [1, 2, 3], 4, 5)
        self.rxr.add_cell(cellK)
        self.sim.add_misc(dsph)
        self.assertEquals(cellK.comment(), "Cell 'K': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "DXC'sph1':N= 0.5.")
        self.assertEquals(cellK.mcnp('%g', self.sim), "12 1 -10 "
                "(-1 2) DXC1:N=0.5")
        # Photon weight
        cellL = cards.CellMCNP('L', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                photon_weight=('one',))
        self.rxr.add_cell(cellL)
        self.assertEquals(cellL.comment(), "Cell 'L': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "PWT= one.")
        self.assertEquals(cellL.mcnp('%g', self.sim), "13 1 -10 (-1 2) "
                "PWT=0")
        # Exception test.
        self.assertRaises(ValueError, cards.CellMCNP, 'L', self.pin.neg &
                self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                photon_weight=('one'))
        cellL2 = cards.CellMCNP('L2', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                photon_weight=(0.5,))
        self.rxr.add_cell(cellL2)
        self.assertEquals(cellL2.comment(), "Cell 'L2': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "PWT= 0.5.")
        self.assertEquals(cellL2.mcnp('%g', self.sim), "14 1 -10 (-1 2) "
                "PWT=0.5")
        cellM = cards.CellMCNP('M', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                photon_weight=(0.5, True))
        self.rxr.add_cell(cellM)
        self.assertEquals(cellM.comment(), "Cell 'M': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "PWT= 0.5 (pre-weighted).")
        self.assertEquals(cellM.mcnp('%g', self.sim), "15 1 -10 "
                "(-1 2) PWT=-0.5")
        # Fission turnoff
        cellN = cards.CellMCNP('N', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                fission_turnoff='capture-nogamma')
        self.rxr.add_cell(cellN)
        self.assertEquals(cellN.comment(), "Cell 'N': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "NONU= capture-nogamma.")
        self.assertEquals(cellN.mcnp('%g', self.sim), "16 1 -10 "
                "(-1 2) NONU=2")
        # Detector contribution
        cellO = cards.CellMCNP('O', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                det_contrib=('det1', 0.5))
        det1 = cards.PointDetector('det1', 'neutron', [0, 0, 0], 0, 'cm')
        self.rxr.add_cell(cellO)
        self.sim.add_tally(det1)
        self.assertEquals(cellO.comment(), "Cell 'O': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "PD for tally 'det1'= 0.5.")
        self.assertEquals(cellO.mcnp('%g', self.sim), "17 1 -10 "
                "(-1 2) PD15=0.5")
        # Transformation
        cellP = cards.CellMCNP('P', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                transformation='transA')
        trans = cards.Transformation('transA', [1, 2, 3],
                [[4, 5, 6], [7, 8, 9], [10, 11, 12]])
        self.rxr.add_cell(cellP)
        self.sim.add_transformation(trans)
        self.assertEquals(cellP.comment(), "Cell 'P': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "TRCL 'transA'.")
        self.assertEquals(cellP.mcnp('%g', self.sim), "18 1 -10 "
                "(-1 2) TRCL=1")
        cellP = cards.CellMCNP('P', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                transformation=([1, 0, 0], np.eye(3)))
        self.rxr.add_cell(cellP)
        self.assertEquals(cellP.comment(), "Cell 'P': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "TRCL aux origin in main (1, 0, 0) cm, x' <1, 0, 0>, "
                "y' <0, 1, 0>, z' <0, 0, 1>.")
        self.assertEquals(cellP.mcnp('%g', self.sim), "18 1 -10 "
                "(-1 2) TRCL ( 1 0 0 1 0 0 0 1 0 0 0 1 1)")
        cellP = cards.CellMCNP('P', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                transformation=([1, 0, 0], np.eye(3), True, True))
        self.rxr.add_cell(cellP)
        self.assertEquals(cellP.comment(), "Cell 'P': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "TRCL aux origin in main (1, 0, 0) cm, "
                "x' <1, 0, 0> deg, y' <0, 1, 0> deg, "
                "z' <0, 0, 1> deg.")
        self.assertEquals(cellP.mcnp('%g', self.sim), "18 1 -10 "
                "(-1 2) *TRCL ( 1 0 0 1 0 0 0 1 0 0 0 1 1)")
        cellP = cards.CellMCNP('P', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                user_custom='EXT:N 0.7V2')
        self.rxr.add_cell(cellP)
        self.assertEquals(cellP.comment(), "Cell 'P': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "and user's custom input.")
        self.assertEquals(cellP.mcnp('%g', self.sim), "18 1 -10 "
                "(-1 2) EXT:N 0.7V2")
        #TODO set temperature to 100 or -1
        #cellE = cards.CellMCNP('E', self.pin.neg & self.cellbound.pos,
        #        self.uo2, 10.0, 'g/cm^3',
        #        temperature=100)
        #self.assertRaises(UserWarning, cellE.comment)

    def test_ExponentialTransform(self):
        """Tests :py:class:`cards.ExponentialTransform` and the related
        :py:class:`cards.Vector`.

        """
        ## Vector
        vec = cards.Vector()
        self.assertEquals(vec.name, 'vector')
        vec.set('origin', [0, 0, 0])
        self.assertEquals(vec.index('origin'), 0)
        self.assertEquals(vec.comment(), "Vector 'vector': origin: "
                "(0, 0, 0) cm.")
        # the mcnp() method doesn't actually need a sim.
        self.assertEquals(vec.mcnp('%.1e', None), "VECT V0 0.0e+00 0.0e+00 "
                "0.0e+00")
        vec.set('x-axis', np.array([1, 0, 0]))
        self.assertEquals(vec.index('origin'), 0)
        self.assertEquals(vec.index('x-axis'), 1)
        self.assertEquals(vec.comment(), "Vector 'vector': origin: "
                "(0, 0, 0) cm, x-axis: "
                "(1, 0, 0) cm.")
        # the mcnp() method doesn't actually need a sim.
        self.assertEquals(vec.mcnp('%.1e', None), "VECT V0 0.0e+00 0.0e+00 "
                "0.0e+00 V1 1.0e+00 0.0e+00 0.0e+00")
        self.sim.add_misc(vec)
        # Try modifying.
        vec.set('x-axis', [2, 1, 3])
        self.assertEquals(vec.mcnp('%.1e', None), "VECT V0 0.0e+00 0.0e+00 "
                "0.0e+00 V1 2.0e+00 1.0e+00 3.0e+00")

        ## ExponentialTransform
        extn = cards.ExponentialTransform('neutron', self.fuel, 'capture-to-total',
                'currdir', 'toward')
        self.assertEquals(extn.name, 'exptransform-neutron')
        self.assertIs(extn.cells[0], self.fuel)
        self.assertEquals(extn.stretchs[self.fuel], 'capture-to-total')
        self.assertEquals(extn.directions[self.fuel], 'currdir')
        self.assertEquals(extn.signs[self.fuel], 'toward')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "capture-to-total toward currdir.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), "EXT:N S 0 0")
        extn = cards.ExponentialTransform('neutron', self.coolant, 0.5,
                'currdir', 'toward')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'coolant' stretch by "
                "0.5 toward currdir.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), "EXT:N 0 5.0e-01 0")
        extn = cards.ExponentialTransform('neutron', self.fuel, 0.5, 'x',
                'toward')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "0.5 toward x.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), "EXT:N 5.0e-01X 0 0")
        extn = cards.ExponentialTransform('neutron', self.fuel, 0.5, 'origin',
                'away')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "0.5 away from origin.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), "EXT:N -5.0e-01V0 0 0")
        # Multiple cells.
        extn = cards.ExponentialTransform('neutron', 
                self.fuel, 'capture-to-total', 'currdir', 'toward', 
                self.coolant, 0.5, 'currdir', 'toward',
                self.graveyard, 0.5, 'x-axis', 'away')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "capture-to-total toward currdir; cell 'coolant' stretch by "
                "0.5 toward currdir; cell 'graveyard' stretch by 0.5 away "
                "from x-axis.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), 
                "EXT:N S 5.0e-01 -5.0e-01V1")
        # Using set().
        ext2 = cards.ExponentialTransform('neutron')
        ext2.set(self.fuel, 'capture-to-total', 'currdir', 'toward')
        ext2.set(self.coolant, 0.5, 'currdir', 'toward')
        ext2.set(self.graveyard, 0.5, 'x-axis', 'away')
        self.assertEquals(ext2.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "capture-to-total toward currdir; cell 'coolant' stretch by "
                "0.5 toward currdir; cell 'graveyard' stretch by 0.5 away "
                "from x-axis.")
        self.assertEquals(ext2.mcnp('%.1e', self.sim), 
                "EXT:N S 5.0e-01 -5.0e-01V1")
        # Modifying.
        ext2.set(self.graveyard, 0.7, 'x-axis', 'away')
        self.assertEquals(ext2.mcnp('%.1e', self.sim), 
                "EXT:N S 5.0e-01 -7.0e-01V1")
        # Manage order-of-cells issue/question.
        # TODO editing entries.

    def test_ForcedCollision(self):
        ## ForcedCollision
        fcl = cards.ForcedCollision('neutron', self.fuel, 0.5, True)
        self.assertEquals(fcl.name, 'forcedcoll-neutron')
        self.assertEquals(fcl.cells[0], self.fuel)
        self.assertEquals(fcl.probs[self.fuel], 0.5)
        self.assertTrue(fcl.only_enterings[self.fuel])
        self.assertEquals(fcl.comment(), "Forced collision "
                "'forcedcoll-neutron': cell 'fuel' prob 0.5 for entering only.")
        self.assertEquals(fcl.mcnp('%.1e', self.sim), "FCL:N -5.0e-01 0 0")
        fcl = cards.ForcedCollision('neutron', self.fuel, 0.5, False)
        self.assertEquals(fcl.comment(), "Forced collision "
                "'forcedcoll-neutron': cell 'fuel' prob 0.5 for entering "
                "and weight games.")
        self.assertEquals(fcl.mcnp('%.1e', self.sim), "FCL:N 5.0e-01 0 0")
        fcl = cards.ForcedCollision('neutron', self.fuel, 0.5, True,
                                               self.coolant, 0.5, False)
        self.assertEquals(fcl.comment(), "Forced collision "
                "'forcedcoll-neutron': cell 'fuel' prob 0.5 for entering "
                "only; cell 'coolant' prob 0.5 for entering and weight games.")
        self.assertEquals(fcl.mcnp('%.1e', self.sim), "FCL:N -5.0e-01 5.0e-01 0")
        fcl = cards.ForcedCollision('neutron')
        fcl.set(self.fuel, 0.5, True)
        fcl.set(self.coolant, 0.5, False)
        self.assertEquals(fcl.comment(), "Forced collision "
                "'forcedcoll-neutron': cell 'fuel' prob 0.5 for entering "
                "only; cell 'coolant' prob 0.5 for entering and weight games.")
        self.assertEquals(fcl.mcnp('%.1e', self.sim), "FCL:N -5.0e-01 5.0e-01 0")
        fcl.set(self.coolant, 0.7, False)
        self.assertEquals(fcl.comment(), "Forced collision "
                "'forcedcoll-neutron': cell 'fuel' prob 0.5 for entering "
                "only; cell 'coolant' prob 0.7 for entering and weight games.")
        self.assertEquals(fcl.mcnp('%.1e', self.sim), "FCL:N -5.0e-01 7.0e-01 0")

    def test_WeightWindows(self):
        """Tests weight window cards."""

        # Tests float conversion.
        wwe = cards.WeightWindowEnergies('photon', [1, 10])
        self.assertEquals(wwe.name, 'weightwinenergy-photon')
        self.assertEquals(wwe.comment(), "Weight window energies "
                "'weightwinenergy-photon' for photons: 2 bins.")
        self.assertEquals(wwe.mcnp('%.1g', self.sim), "WWE:P 1 1e+01")
        wwe = cards.WeightWindowEnergies('photon', [1, 10], for_gen=True)
        self.assertEquals(wwe.name, 'weightwingenenergy-photon')
        self.assertEquals(wwe.comment(), "Weight window generator energies "
                "'weightwingenenergy-photon' for photons: 2 bins.")
        self.assertEquals(wwe.mcnp('%.1g', self.sim), "WWGE:P 1 1e+01")
        wwe = cards.WeightWindowEnergies('photon', [], for_gen=True)
        self.assertEquals(wwe.name, 'weightwingenenergy-photon')
        self.assertEquals(wwe.comment(), "Weight window generator energies "
                "'weightwingenenergy-photon' for photons: default bins.")
        self.assertEquals(wwe.mcnp('%.1g', self.sim), "WWGE:P")
        wwt = cards.WeightWindowTimes('proton', [1, 1e12])
        self.assertEquals(wwt.name, 'weightwintime-proton')
        self.assertEquals(wwt.comment(), "Weight window times "
                "'weightwintime-proton' for protons: 2 intervals.")
        self.assertEquals(wwt.mcnp('%.1g', self.sim), "WWT:H 1e+08 1e+20")
        wwt = cards.WeightWindowTimes('proton', [1, 1e12], for_gen=True)
        self.assertEquals(wwt.name, 'weightwingentime-proton')
        self.assertEquals(wwt.comment(), "Weight window generator times "
                "'weightwingentime-proton' for protons: 2 intervals.")
        self.assertEquals(wwt.mcnp('%.1g', self.sim), "WWGT:H 1e+08 1e+20")
        wwt = cards.WeightWindowTimes('proton', [], for_gen=True)
        self.assertEquals(wwt.name, 'weightwingentime-proton')
        self.assertEquals(wwt.comment(), "Weight window generator times "
                "'weightwingentime-proton' for protons: default intervals.")
        self.assertEquals(wwt.mcnp('%.1g', self.sim), "WWGT:H")
        
        
        cellC = cards.Cell('C', self.pin.neg & self.cellbound.pos)
        cellD = cards.Cell('D', self.pin.neg & self.cellbound.pos)
        cellE = cards.Cell('E', self.pin.neg & self.cellbound.pos)
        self.sim.sys.add_cell(cellC)
        self.rxr.add_cell(cellD)
        self.rxr.add_cell(cellE)

        wwe = cards.WeightWindowEnergies('photon', range(1, 5), for_gen=True)
        wwt = cards.WeightWindowTimes('photon', range(1, 13), for_gen=True)
        self.sim.add_misc(wwe)
        self.sim.add_misc(wwt)

        wwn = cards.WeightWindowBound('photon', 1, 1, cellD, 0.01)
        self.assertEquals(wwn.name, 'weightwinbound-photon')
        self.assertEquals(wwn.comment(), "Weight window bounds "
                "'weightwinbound-photon' for photons: energy idx 1: "
                "time idx 1: cell 'D': 0.01.")
        self.assertEquals(wwn.mcnp('%g', self.sim), "WWN1:P 0 0 0 0 0.01 0")
        # More than one cell.
        wwn = cards.WeightWindowBound('photon' , 1, 1, cellD, 0.01 , 1, 1,
                cellE, 0.02 , 2, 1, cellE, 0.03 , 1, 3, cellD, 0.04 , 2, 3,
                cellD, 0.05 , 2, 3, cellE, 0.06)
        self.assertEquals(wwn.comment(), "Weight window bounds "
                "'weightwinbound-photon' for photons: energy idx 1: "
                "time idx 1: cell 'D': 0.01, cell 'E': 0.02, "
                "time idx 3: cell 'D': 0.04, "
                "energy idx 2: time idx 1: cell 'E': 0.03, "
                "time idx 3: cell 'D': 0.05, cell 'E': 0.06.")
        self.assertEquals(wwn.mcnp('%g', self.sim), 
                "WWN1:P 0 0 0 0 0.01 0.02\n"
                "WWN2:P 0 0 0 0 0 0.03\n"
                "WWN9:P 0 0 0 0 0.04 0\n"
                "WWN10:P 0 0 0 0 0.05 0.06")
        # Make sure set() works.
        wwn = cards.WeightWindowBound('photon')
        wwn.set(1, 1, cellD, 0.01)
        wwn.set(1, 1, cellE, 0.02)
        wwn.set(2, 1, cellE, 0.03)
        wwn.set(1, 3, cellD, 0.04)
        wwn.set(2, 3, cellD, 0.05)
        wwn.set(2, 3, cellE, 0.06)
        self.assertEquals(wwn.comment(), "Weight window bounds "
                "'weightwinbound-photon' for photons: energy idx 1: "
                "time idx 1: cell 'D': 0.01, cell 'E': 0.02, "
                "time idx 3: cell 'D': 0.04, "
                "energy idx 2: time idx 1: cell 'E': 0.03, "
                "time idx 3: cell 'D': 0.05, cell 'E': 0.06.")
        self.assertEquals(wwn.mcnp('%g', self.sim), 
                "WWN1:P 0 0 0 0 0.01 0.02\n"
                "WWN2:P 0 0 0 0 0 0.03\n"
                "WWN9:P 0 0 0 0 0.04 0\n"
                "WWN10:P 0 0 0 0 0.05 0.06")

        # Make sure using WWGT default times gives the proper linear index.
        wwe = cards.WeightWindowEnergies('photon', [], for_gen=True)
        # Overwriting a card.
        self.sim.add_misc(wwe)
        #self.assertRaises(UserWarning, self.sim.add_misc, wwe)
        wwn = cards.WeightWindowBound('photon', 8, 6, cellD, 0.01)
        self.assertEquals(wwn.comment(), "Weight window bounds "
                "'weightwinbound-photon' for photons: energy idx 8: "
                "time idx 6: cell 'D': 0.01.")
        self.assertEquals(wwn.mcnp('%g', self.sim), "WWN58:P 0 0 0 0 0.01 0")

        # Test killall
        wwn = cards.WeightWindowBound('photon', 1, 1, cellD, 'killall')
        self.assertEquals(wwn.comment(), "Weight window bounds "
                "'weightwinbound-photon' for photons: energy idx 1: "
                "time idx 1: cell 'D': killall.")
        self.assertEquals(wwn.mcnp('%g', self.sim), "WWN1:P 0 0 0 0 -1 0")

        # Check exception when using a particle type for which WWE/WWT cards
        # are not defined.
        wwn = cards.WeightWindowBound('neutron', 1, 1, cellD, 0.01)
        self.assertRaises(Exception, wwn.mcnp,'%g', self.sim)
        try:
            wwn.mcnp('%g', self.sim)
        except Exception as e:
            self.assertEquals(e.message, "No WWGT:N or WWT:N card found in "
                    "the simulation.")
        ## Importance
        impn = cards.Importance('neutron', self.fuel, 1, self.coolant, 2)
        self.assertEquals(impn.name, 'importance-neutron')
        self.assertEquals(impn.comment(), "Importance 'importance-neutron': "
                "cell 'fuel' 1; cell 'coolant' 2.")
        self.assertEquals(impn.mcnp('%g', self.sim), "IMP:N 1 2 0 0 0 0")
        impn.set(self.graveyard, 0)
        self.assertEquals(impn.comment(), "Importance 'importance-neutron': "
                "cell 'fuel' 1; cell 'coolant' 2; cell 'graveyard' 0.")
        self.assertEquals(impn.mcnp('%g', self.sim), "IMP:N 1 2 0 0 0 0")
        # Modifying.
        impn.set(self.graveyard, 1)
        self.assertEquals(impn.comment(), "Importance 'importance-neutron': "
                "cell 'fuel' 1; cell 'coolant' 2; cell 'graveyard' 1.")
        self.assertEquals(impn.mcnp('%g', self.sim), "IMP:N 1 2 1 0 0 0")
        # Using args.
        args = [self.fuel, 1, self.coolant, 2, self.graveyard, 3]
        impn = cards.Importance('neutron', *args)
        self.assertEquals(impn.comment(), "Importance 'importance-neutron': "
                "cell 'fuel' 1; cell 'coolant' 2; cell 'graveyard' 3.")
        self.assertEquals(impn.mcnp('%g', self.sim), "IMP:N 1 2 3 0 0 0")

        ## Volume
        vol = cards.Volume(self.fuel, 1)
        self.assertEquals(vol.name, 'volume')
        self.assertEquals(vol.comment(), "Volume 'volume': cell 'fuel' 1 cm^3.")
        self.assertEquals(vol.mcnp('%g', self.sim), "VOL 1 5J")
        vol = cards.Volume(cellD, 1)
        self.assertEquals(vol.mcnp('%g', self.sim), "VOL 4J 1 1J")
        # Multiple cells
        vol = cards.Volume(self.fuel, 1, self.coolant, 2, manual=True)
        self.assertEquals(vol.comment(), "Volume 'volume': (all manual) "
                "cell 'fuel' 1 cm^3, cell 'coolant' 2 cm^3.")
        self.assertEquals(vol.mcnp('%g', self.sim), "VOL NO 1 2 4J")
        # Using set()
        vol.set(self.coolant, 3)
        self.assertEquals(vol.mcnp('%g', self.sim), "VOL NO 1 3 4J")

        ## Area
        are = cards.Area(self.fuel, 10)
        self.assertEquals(are.name, 'area')
        self.assertEquals(are.comment(), "Area 'area': cell 'fuel' 10 cm^2.")
        self.assertEquals(are.mcnp('%g', self.sim), "AREA 10 5J")
        are = cards.Area(cellD, 10)
        self.assertEquals(are.mcnp('%g', self.sim), "AREA 4J 10 1J")
        # Multiple cells
        are = cards.Area(self.fuel, 10, self.coolant, 20)
        self.assertEquals(are.comment(), "Area 'area': "
                "cell 'fuel' 10 cm^2, cell 'coolant' 20 cm^2.")
        self.assertEquals(are.mcnp('%g', self.sim), "AREA 10 20 4J")
        # Using set()
        are.set(self.coolant, 30)
        self.assertEquals(are.mcnp('%g', self.sim), "AREA 10 30 4J")

        ## TemperatureTimes
        thtme = cards.TemperatureTimes([1e10, 2e10])
        self.assertEquals(thtme.name, 'temptimes')
        self.assertEquals(thtme.comment(), "Temperature times 'temptimes' "
                "(in MeV): 1e+10 2e+10.")
        self.assertEquals(thtme.mcnp('%g', self.sim), "THTME 1e+18 2e+18")

        ## Temperature
        temp = cards.Temperature(self.fuel, 600)
        self.assertEquals(temp.name, 'temperature')
        self.assertEquals(temp.comment(), "Temperatures 'temperature': "
                "cell 'fuel' 600 K.")
        self.assertEquals(temp.mcnp('%g', self.sim),
                "TMP 5.17041e-08 5J")
        temp = cards.Temperature(self.fuel, 600, cellE, 900, index=2)
        self.assertEquals(temp.comment(), "Temperatures for time index "
                "2 'temperature-idx2': cell 'fuel' 600 K, "
                "cell 'E' 900 K.")
        self.assertEquals(temp.mcnp('%g', self.sim),
                "TMP2 5.17041e-08 4J 7.75561e-08")
        # set()
        temp = cards.Temperature(index=2)
        temp.set(self.fuel, 600)
        temp.set(cellE, 900)
        self.assertEquals(temp.comment(), "Temperatures for time index "
                "2 'temperature-idx2': cell 'fuel' 600 K, "
                "cell 'E' 900 K.")
        self.assertEquals(temp.mcnp('%g', self.sim),
                "TMP2 5.17041e-08 4J 7.75561e-08")
        # Modifying.
        temp.set(cellE, 950)
        self.assertEquals(temp.comment(), "Temperatures for time index "
                "2 'temperature-idx2': cell 'fuel' 600 K, "
                "cell 'E' 950 K.")
        self.assertEquals(temp.mcnp('%g', self.sim),
                "TMP2 5.17041e-08 4J 8.18648e-08")

        ## DXTRANSpheres
        dsph = cards.DXTRANSpheres('neutron', 'sph1', [1, 2, 3], 4, 5)
        self.assertEquals(dsph.name, 'dxtranspheres-neutron')
        self.assertEquals(dsph.comment(), "DXTRAN spheres "
                "'dxtranspheres-neutron': up. cut. default, low cut. default"
                ", min photon wgt. default. sphere 'sph1' at (1, 2, 3) cm, "
                "in. rad. 4 cm, out. rad. 5 cm.")
        self.assertEquals(dsph.mcnp('%g', self.sim),
                "DXT:N  1 2 3 4 5  0 0 0")
        dsph = cards.DXTRANSpheres('neutron', 'sph1', [1, 2, 3], 4, 5,
                                              'sph2', [4, 5, 6], 7, 8,
                                   upper_cutoff=0.1, lower_cutoff=0.05,
                                   min_photon_weight=0.5)
        self.assertEquals(dsph.comment(), "DXTRAN spheres "
                "'dxtranspheres-neutron': up. cut. 0.1, low cut. 0.05"
                ", min photon wgt. 0.5. sphere 'sph1' at (1, 2, 3) cm, "
                "in. rad. 4 cm, out. rad. 5 cm; sphere 'sph2' at (4, 5, 6) cm,"
                " in. rad. 7 cm, out. rad. 8 cm.")
        self.assertEquals(dsph.mcnp('%g', self.sim),
                "DXT:N  1 2 3 4 5  4 5 6 7 8  0.1 0.05 0.5")
        dsph = cards.DXTRANSpheres('neutron', upper_cutoff=0.1,
                                   lower_cutoff=0.05,
                                   min_photon_weight=0.5)
        dsph.set('sph1', [1, 2, 3], 4, 5)
        dsph.set('sph2', [4, 5, 6], 7, 8)
        self.assertEquals(dsph.comment(), "DXTRAN spheres "
                "'dxtranspheres-neutron': up. cut. 0.1, low cut. 0.05"
                ", min photon wgt. 0.5. sphere 'sph1' at (1, 2, 3) cm, "
                "in. rad. 4 cm, out. rad. 5 cm; sphere 'sph2' at (4, 5, 6) cm,"
                " in. rad. 7 cm, out. rad. 8 cm.")
        self.assertEquals(dsph.mcnp('%g', self.sim),
                "DXT:N  1 2 3 4 5  4 5 6 7 8  0.1 0.05 0.5")
        dsph.set('sph2', [4.5, 5.5, 6.5], 8.5, 9.5)
        self.assertEquals(dsph.mcnp('%g', self.sim),
                "DXT:N  1 2 3 4 5  4.5 5.5 6.5 8.5 9.5  0.1 0.05 0.5")
        self.sim.add_misc(dsph)

        # DXTRANContribution
        dxc = cards.DXTRANContribution('neutron', None, self.fuel, 0.5)
        self.assertEquals(dxc.name, 'dxtrancont-neutron')
        self.assertEquals(dxc.comment(), "DXTRAN contribution "
                "for all spheres 'dxtrancont-neutron': cell 'fuel' 0.5.")
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC:N 0.5 5J")
        dxc = cards.DXTRANContribution('neutron', 'sph2', self.fuel, 0.5)
        self.assertEquals(dxc.comment(), "DXTRAN contribution "
                "for sphere 'sph2' 'dxtrancont-sph2-neutron': "
                "cell 'fuel' 0.5.")
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC2:N 0.5 5J")
        dxc = cards.DXTRANContribution('neutron', 'sph2', self.fuel, 0.5,
                                                       self.coolant, 0.7)
        self.assertEquals(dxc.comment(), "DXTRAN contribution "
                "for sphere 'sph2' 'dxtrancont-sph2-neutron': "
                "cell 'fuel' 0.5, cell 'coolant' 0.7.")
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC2:N 0.5 0.7 4J")
        dxc = cards.DXTRANContribution('neutron', 'sph2')
        dxc.set(self.fuel, 0.5)
        dxc.set(self.coolant, 0.7)
        self.assertEquals(dxc.comment(), "DXTRAN contribution "
                "for sphere 'sph2' 'dxtrancont-sph2-neutron': "
                "cell 'fuel' 0.5, cell 'coolant' 0.7.")
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC2:N 0.5 0.7 4J")
        dxc.set(self.coolant, 0.8)
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC2:N 0.5 0.8 4J")

        ## FissionTurnoff
        # Default.
        fto = cards.FissionTurnoff()
        self.assertEquals(fto.name, 'fissionturnoff')
        self.assertEquals(fto.comment(), "Fission turnoff 'fissionturnoff': "
                "capture-gamma for all cells.")
        self.assertEquals(fto.mcnp('%g', self.sim), "NONU")
        # With an argument.
        fto = cards.FissionTurnoff(self.fuel, 'real-gamma',
                                   self.coolant, 'capture-nogamma')
        self.assertEquals(fto.comment(), "Fission turnoff 'fissionturnoff': "
                "cell 'fuel' real-gamma, cell 'coolant' capture-nogamma.")
        self.assertEquals(fto.mcnp('%g', self.sim), "NONU 1 2 4J")
        # set()
        fto = cards.FissionTurnoff()
        fto.set(self.fuel, 'real-gamma')
        fto.set(self.coolant, 'capture-nogamma')
        self.assertEquals(fto.comment(), "Fission turnoff 'fissionturnoff': "
                "cell 'fuel' real-gamma, cell 'coolant' capture-nogamma.")
        self.assertEquals(fto.mcnp('%g', self.sim), "NONU 1 2 4J")
        # Modifying.
        fto.set(self.coolant, 'capture-gamma')
        self.assertEquals(fto.mcnp('%g', self.sim), "NONU 1 0 4J")

        ## DetectorContribution
        det = cards.PointDetector('point', 'neutron', [0, 0, 0], 0, 'cm')
        det2 = cards.PointDetector('point2', 'neutron', [0, 0, 0], 0, 'cm')
        self.sim.add_tally(det)
        self.sim.add_tally(det2)
        dc = cards.DetectorContribution('point2', self.fuel, 0.5)
        self.assertEquals(dc.name, 'detcontrib-point2')
        self.assertEquals(dc.comment(), "Detector contribution "
                "'detcontrib-point2': cell 'fuel' 0.5.")
        self.assertEquals(dc.mcnp('%g', self.sim), "PD25 0.5 5J")
        dc = cards.DetectorContribution('point2', self.fuel, 0.5, self.coolant,
                0.6)
        self.assertEquals(dc.comment(), "Detector contribution "
                "'detcontrib-point2': cell 'fuel' 0.5, cell 'coolant' 0.6.")
        self.assertEquals(dc.mcnp('%g', self.sim), "PD25 0.5 0.6 4J")
        dc = cards.DetectorContribution('point2')
        dc.set(self.fuel, 0.5)
        dc.set(self.coolant, 0.6)
        self.assertEquals(dc.comment(), "Detector contribution "
                "'detcontrib-point2': cell 'fuel' 0.5, cell 'coolant' 0.6.")
        self.assertEquals(dc.mcnp('%g', self.sim), "PD25 0.5 0.6 4J")
        dc.set(self.coolant, 0.7)
        self.assertEquals(dc.mcnp('%g', self.sim), "PD25 0.5 0.7 4J")

        ## PhotonWeight
        pw = cards.PhotonWeight()
        self.assertEquals(pw.name, 'photonweight')
        pw.set(self.fuel, 'off')
        self.assertEquals(pw.comment(), "Photon weight thresholds "
                "'photonweight': cell 'fuel' off.")
        self.assertEquals(pw.mcnp('%g', self.sim), "PWT -1.0E6 5J")
        pw.set(self.coolant, 'one')
        self.assertEquals(pw.comment(), "Photon weight thresholds "
                "'photonweight': cell 'fuel' off, "
                "cell 'coolant' one.")
        self.assertEquals(pw.mcnp('%g', self.sim), "PWT -1.0E6 0 4J")
        pw.set(self.coolant, 0.5)
        self.assertEquals(pw.comment(), "Photon weight thresholds "
                "'photonweight': cell 'fuel' off, "
                "cell 'coolant' 0.5.")
        self.assertEquals(pw.mcnp('%g', self.sim), "PWT -1.0E6 0.5 4J")
        pw.set(self.coolant, 0.5, pre_weight=True)
        self.assertEquals(pw.comment(), "Photon weight thresholds "
                "'photonweight': cell 'fuel' off, "
                "cell 'coolant' 0.5 (pre-weighted).")
        self.assertEquals(pw.mcnp('%g', self.sim), "PWT -1.0E6 -0.5 4J")

    def test_Transformation(self):
        """Tests :py:class:`cards.Transformation`."""

        trans = cards.Transformation('trans1', [1, 2, 3],
                [[4, 5, 6], [7, 8, 9], [10, 11, 12]])
        self.sim.add_transformation(trans)
        self.assertEquals(trans.comment(), "Transformation 'trans1': "
                "pos. of aux origin in main (1, 2, 3) cm, "
                "x' <4, 7, 10>, y' <5, 8, 11>, z' <6, 9, 12>.")
        self.assertEquals(trans.mcnp('%.1f', self.sim), "TR1 "
                "1.0 2.0 3.0 "
                "4.0 5.0 6.0 "
                "7.0 8.0 9.0 "
                "10.0 11.0 12.0 1")
        trans = cards.Transformation('trans1', [1, 2, 3],
                np.array([[4, 5, 6], [7, 8, 9], [10, 11, 12]]))
        self.assertEquals(trans.comment(), "Transformation 'trans1': "
                "pos. of aux origin in main (1, 2, 3) cm, "
                "x' <4, 7, 10>, y' <5, 8, 11>, z' <6, 9, 12>.")
        self.assertEquals(trans.mcnp('%.1f', self.sim), "TR1 "
                "1.0 2.0 3.0 "
                "4.0 5.0 6.0 "
                "7.0 8.0 9.0 "
                "10.0 11.0 12.0 1")
        trans = cards.Transformation('trans1', [1, 2, 3],
                np.matrix([[4, 5, 6], [7, 8, 9], [10, 11, 12]]))
        self.assertEquals(trans.comment(), "Transformation 'trans1': "
                "pos. of aux origin in main (1, 2, 3) cm, "
                "x' <4, 7, 10>, y' <5, 8, 11>, z' <6, 9, 12>.")
        self.assertEquals(trans.mcnp('%.1f', self.sim), "TR1 "
                "1.0 2.0 3.0 "
                "4.0 5.0 6.0 "
                "7.0 8.0 9.0 "
                "10.0 11.0 12.0 1")
        trans = cards.Transformation('trans1', [1, 2, 3],
                [[4, 5, 6], [7, 8, 9], [10, 11, 12]],
                aux_in_main=False)
        self.assertFalse(trans.aux_in_main)
        self.assertEquals(trans.comment(), "Transformation 'trans1': "
                "pos. of main origin in aux (1, 2, 3) cm, "
                "x' <4, 7, 10>, y' <5, 8, 11>, z' <6, 9, 12>.")
        self.assertEquals(trans.mcnp('%.1f', self.sim), "TR1 "
                "1.0 2.0 3.0 "
                "4.0 5.0 6.0 "
                "7.0 8.0 9.0 "
                "10.0 11.0 12.0 -1")
        trans = cards.Transformation('trans1', [1, 2, 3],
                [[4, 5, 6], [7, 8, 9], [10, 11, 12]],
                degrees=True)
        self.assertEquals(trans.comment(), "Transformation 'trans1': "
                "pos. of aux origin in main (1, 2, 3) cm, "
                "x' <4, 7, 10> deg, y' <5, 8, 11> deg, z' <6, 9, 12> deg.")
        self.assertEquals(trans.mcnp('%.1f', self.sim), "*TR1 "
                "1.0 2.0 3.0 "
                "4.0 5.0 6.0 "
                "7.0 8.0 9.0 "
                "10.0 11.0 12.0 1")

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
                    "must be an integer. User provided 0.5.")
        try:
            critsrc.n_histories = -1
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_histories`` "
                    "must be positive. User provided -1.")
        try:
            critsrc.keff_guess = -1
        except ValueError as e:
            self.assertEquals(e.message, "The property ``keff_guess`` "
                    "must be non-negative. User provided -1.")
        try:
            critsrc.n_skip_cycles = 0.5
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_skip_cycles`` "
                    "must be an integer. User provided 0.5.")
        try:
            critsrc.n_skip_cycles = -1
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_skip_cycles`` "
                    "must be positive. User provided -1.")
        try:
            critsrc.n_cycles = 0.5
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_cycles`` "
                    "must be an integer. User provided 0.5.")
        try:
            critsrc.n_cycles = -1
        except ValueError as e:
            self.assertEquals(e.message, "The property ``n_cycles`` "
                    "must be equal to or greater than ``n_skip_cycles``. "
                    "User provided -1.")

        ## comment()
        critsrc = cards.Criticality()
        self.assertEquals(critsrc.comment(), "Criticality source "
                "'criticality': n_histories: 1000, keff_guess: 1"
                ", n_skip_cycles: 30, n_cycles: 130.")

        ## mcnp()
        self.assertEquals(critsrc.mcnp('%g', self.sim), 
                "KCODE 1000 1 30 130")

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
                                        np.array([3.1416, 2.7183, 0])])

        ## comment()
        self.assertEquals(critpts.comment(), "Criticality points "
                "'criticalitypoints': (1, 2, 3) cm, "
                "(3.1416, 2.7183, 0) cm.")
        self.assertRaises(ValueError, setattr, critpts, 'points', 
                [[0, 0, 0], [1]])
        try:
            critpts.points = [[0, 0, 0], [1]]
        except ValueError as e:
            self.assertEquals(e.message, "Length of all point lists/arrays "
                    "must be 3. User provided a point [1].")

        ## mcnp()
        self.assertEquals(critpts.mcnp('%g', self.sim),
                "KSRC 1 2 3\n"
                "     3.1416 2.7183 0")

    def test_ITally(self):
        """Tests :py:class:`cards.ITally`'s methods, properties, and
        exceptions, and those of its subclasses.

        """
        ## SurfaceCurrent
        tally = cards.SurfaceCurrent('fuel', 'electron', [self.pin,
                self.cellbound], total=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle.name, 'electron')
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
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim),
                "F11:E  1 2 T")

        tally = cards.SurfaceCurrent('fuel2', 'photon', [[self.pin,
                self.cellbound]], alt_units=True)
        self.assertEquals(tally.particle.name, 'photon')
        self.assertTrue(len(tally.cards), 1)
        self.assertTrue(len(tally.cards[0]), 2)
        self.assertIs(tally.cards[0][0], self.pin)
        self.assertIs(tally.cards[0][1], self.cellbound)
        self.assertFalse(tally.total)
        self.assertTrue(tally.alt_units)
        # comment()
        self.assertEquals(tally.comment(), "Surface current tally 'fuel2' "
                "of photons (in alt. units): total in 'fuelpin', 'bound'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim),
                "*F21:P  ( 1 2)")

        ## SurfaceFlux
        tally = cards.SurfaceFlux('fuel', 'electron', [self.pin,
                self.cellbound], average=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle.name, 'electron')
        self.assertIs(tally.cards[0], self.pin)
        self.assertIs(tally.cards[1], self.cellbound)
        self.assertTrue(len(tally.cards), 2)
        self.assertTrue(tally.average)
        self.assertFalse(tally.alt_units)
        # comment()
        self.assertEquals(tally.comment(), "Surface flux tally 'fuel' "
                "of electrons: surfaces 'fuelpin'; 'bound'; and avg. "
                "of all provided.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F12:E  1 2 T")

        tally = cards.SurfaceFlux('fuel2', 'proton', [[self.pin, self.cellbound]],
                    alt_units=True)
        self.assertEquals(tally.particle.name, 'proton')
        self.assertTrue(len(tally.cards), 1)
        self.assertTrue(len(tally.cards[0]), 2)
        self.assertIs(tally.cards[0][0], self.pin)
        self.assertIs(tally.cards[0][1], self.cellbound)
        self.assertFalse(tally.average)
        self.assertTrue(tally.alt_units)
        # comment()
        self.assertEquals(tally.comment(), "Surface flux tally 'fuel2' "
                "of protons (in alt. units): avg. in 'fuelpin', 'bound'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "*F22:H  ( 1 2)")

        ## CellFlux
        # One cell.
        tally = cards.CellFlux('fuel', 'neutron', self.fuel)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle.name, 'neutron')
        self.assertEquals(tally.comment(), "Cell flux tally 'fuel' "
                "of neutrons: cell 'fuel'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F14:N  1")
        self.assertIs(tally._unique_card_list()[0], self.fuel)
        # Two individual cells.
        tally = cards.CellFlux('both', 'neutron', [self.fuel, self.coolant])
        self.assertEquals(tally.comment(), "Cell flux tally 'both' "
                "of neutrons: cells 'fuel'; 'coolant'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F24:N  1 2")
        self.assertIs(tally._unique_card_list()[0], self.fuel)
        self.assertIs(tally._unique_card_list()[1], self.coolant)
        # Two individual cells, with average over all.
        tally = cards.CellFlux('withavg', 'neutron', [self.fuel, self.coolant],
                average=True)
        self.assertEquals(tally.comment(), "Cell flux tally 'withavg' "
                "of neutrons: cells 'fuel'; 'coolant'; and avg. of all "
                "provided.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F34:N  1 2 T")
        self.assertIs(tally._unique_card_list()[0], self.fuel)
        self.assertIs(tally._unique_card_list()[1], self.coolant)
        # Two individual cells, and an averaging, with an average over all.
        tally = cards.CellFlux('withavg', 'neutron', [self.fuel,
                [self.fuel, self.coolant], self.coolant], average=True)
        self.assertEquals(tally.comment(), "Cell flux tally 'withavg' "
                "of neutrons: cells 'fuel'; avg. in 'fuel', 'coolant'; "
                "'coolant'; and avg. of all provided.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F34:N  1 ( 1 2) 2 T")
        self.assertIs(tally._unique_card_list()[0], self.fuel)
        self.assertIs(tally._unique_card_list()[1], self.coolant)
        self.assertTrue(len(tally._unique_card_list()) == 2)

        ## CellEnergyDeposition
        tally = cards.CellEnergyDeposition('energy', 'neutron', self.fuel)
        self.assertEquals(tally.name, 'energy')
        self.assertEquals(tally.particle.name, 'neutron')
        self.assertIs(tally.cards, self.fuel)
        # comment()
        self.assertEquals(tally.comment(), "Energy deposition tally "
                "'energy' of neutrons: cell 'fuel'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F16:N  1")
        tally = cards.CellEnergyDeposition('energy2', ['neutron', 'proton'],
                self.fuel)
        self.assertIs(type(tally.particle), list)
        self.assertEquals(tally.particle[0].name, 'neutron')
        self.assertEquals(tally.particle[1].name, 'proton')
        self.assertEquals(tally.comment(), "Energy deposition tally " 
                "'energy2' of neutrons, protons: cell 'fuel'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F26:N,H  1")
        tally = cards.CellEnergyDeposition('energy3', 'all', self.fuel)
        self.assertEquals(tally.comment(), "Energy deposition tally "
                "'energy3' of all: cell 'fuel'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "+F36  1")
        # Test exceptions.
        self.assertRaises(ValueError, cards.CellEnergyDeposition, 'energy',
                ['neutron', 'all'], self.fuel)
        self.assertRaises(ValueError, cards.CellEnergyDeposition,
                'energy', 'all', self.fuel, alt_units=True)

        ## CellFissionEnergyDeposition
        tally = cards.CellFissionEnergyDeposition('fuel', [self.fuel,
                self.coolant], average=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle.name, 'neutron')
        self.assertIs(type(tally.cards), list)
        self.assertIs(tally.cards[0], self.fuel)
        self.assertIs(tally.cards[1], self.coolant)
        self.assertTrue(tally.average)
        self.assertFalse(tally.alt_units)
        self.assertEquals(tally.comment(), "Fission energy deposition tally "
                "'fuel' of neutrons: cells 'fuel'; 'coolant'; and avg. of "
                "all provided.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F17:N  1 2 T")
        tally = cards.CellFissionEnergyDeposition('fuel', [[self.fuel,
                self.coolant]], alt_units=True)
        self.assertTrue(len(tally.cards), 1)
        self.assertTrue(len(tally.cards[0]), 2)
        self.assertIs(tally.cards[0][0], self.fuel)
        self.assertIs(tally.cards[0][1], self.coolant)
        self.assertFalse(tally.average)
        self.assertTrue(tally.alt_units)
        self.assertEquals(tally.comment(), "Fission energy deposition tally "
                "'fuel' of neutrons (in alt. units): avg. in 'fuel', "
                "'coolant'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "*F17:N  ( 1 2)")

        ## CellPulseHeight
        tally = cards.CellPulseHeight('fuel', ['proton', 'electron'], [self.fuel,
                self.coolant], alt_units=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle[0].name, 'proton')
        self.assertEquals(tally.particle[1].name, 'electron')
        self.assertIs(tally.cards[0], self.fuel)
        self.assertIs(tally.cards[1], self.coolant)
        self.assertFalse(tally.average)
        self.assertTrue(tally.alt_units)
        tally.average = True
        self.assertTrue(tally.average)
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "*F18:H,E  1 2 T")

        ## CellChargeDeposition
        tally = cards.CellChargeDeposition('fuel', ['proton', 'electron'],
                [self.fuel, self.coolant])
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle[0].name, 'proton')
        self.assertEquals(tally.particle[1].name, 'electron')
        self.assertIs(tally.cards[0], self.fuel)
        self.assertIs(tally.cards[1], self.coolant)
        self.assertFalse(tally.average)
        self.assertFalse(tally.alt_units)
        tally.average = True
        self.assertTrue(tally.average)
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "+F18:H,E  1 2 T")

        ## PointDetector
        det = cards.PointDetector('point', 'neutron', [0, 0, 0], 0, 'cm')
        self.assertEquals(det.name, 'point')
        self.assertEquals(det.particle.name, 'neutron')
        self.assertEquals(det.detectors[0].pos, [0, 0, 0])
        self.assertEquals(det.detectors[0].soe_rad, 0)
        self.assertEquals(det.sep_direct, True)
        self.assertEquals(det.comment(), "Point detector tally 'point' of "
                "neutrons: point (0, 0, 0) cm, "
                "radius 0 cm; direct contrib is separate.")
        # mcnp()
        self.sim.add_tally(det)
        self.assertEquals(det.mcnp('%.5g', self.sim), "F15:N  0 0 0 0")
        det = cards.PointDetector('point2', 'neutron', np.array([1, 2, 3]), 4,
                'cm')
        self.assertTrue((det.detectors[0].pos == np.array([1, 2, 3])).all())
        self.assertEquals(det.detectors[0].soe_rad, 4)
        self.assertEquals(det.sep_direct, True)
        self.assertEquals(det.comment(), "Point detector tally 'point2' of "
                "neutrons: point (1, 2, 3) cm, "
                "radius 4 cm; direct contrib is separate.")
        # mcnp()
        self.sim.add_tally(det)
        self.assertEquals(det.mcnp('%.5g', self.sim), "F25:N  1 2 3 4")
        det = cards.PointDetector('point', 'neutron', [1, 0, 0], 3, 'mfp')
        self.assertEquals(det.comment(), "Point detector tally 'point' of "
                "neutrons: point (1, 0, 0) cm, "
                "radius 3 mfp; direct contrib is separate.")
        # mcnp()
        self.sim.add_tally(det)
        self.assertEquals(det.mcnp('%.5g', self.sim), "F15:N  1 0 0 -3")
        det = cards.PointDetector('point', 'photon', [0, 0, 0], 0, 'cm',
                                                     [1, 0, 0], 3, 'mfp')
        self.assertEquals(det.detectors[0].pos, [0, 0, 0])
        self.assertEquals(det.detectors[0].soe_rad, 0)
        self.assertEquals(det.detectors[1].pos, [1, 0, 0])
        self.assertEquals(det.detectors[1].soe_rad, 3)
        self.assertEquals(det.comment(), "Point detector tally 'point' of "
                "photons: point (0, 0, 0) cm, "
                "radius 0 cm; "
                "point (1, 0, 0) cm, radius 3 mfp; "
                "direct contrib is separate.")
        # mcnp()
        self.sim.add_tally(det)
        self.assertEquals(det.mcnp('%.5g', self.sim),
                "F15:P  0 0 0 0\n     1 0 0 -3")
        det = cards.PointDetector('point', 'photon', [0, 0, 0], 0, 'cm',
                sep_direct=False)
        self.assertFalse(det.sep_direct)
        # mcnp()
        self.sim.add_tally(det)
        self.assertEquals(det.mcnp('%.5g', self.sim),
                "F15:P  0 0 0 0 ND")

        ## RingDetector
        det = cards.RingDetector('ring', 'neutron', 'x', 10.0, 2.0, 1.0, 'cm')
        self.assertEquals(det.name, 'ring')
        self.assertEquals(det.particle.name, 'neutron')
        self.assertTrue(det.sep_direct)
        self.assertEquals(det.comment(), "Ring detector tally 'ring' of "
                "neutrons: along x axis. ring x = 10 cm, radius 2 cm, s.o.e. "
                "radius 1 cm; direct contrib is separate.")
        # mcnp()
        self.sim.add_tally(det)
        self.assertEquals(det.mcnp('%.5g', self.sim),
                "F35X:N  10 2 1")
        det = cards.RingDetector('ring', 'neutron', 'x', 10.0, 2.0, 1.0, 'mfp')
        self.assertEquals(det.comment(), "Ring detector tally 'ring' of "
                "neutrons: along x axis. ring x = 10 cm, radius 2 cm, s.o.e. "
                "radius 1 mfp; direct contrib is separate.")
        # mcnp()
        self.sim.add_tally(det)
        self.assertEquals(det.mcnp('%.5g', self.sim),
                "F35X:N  10 2 -1")
        det = cards.RingDetector('ring', 'neutron', 'x', 10.0, 2.0, 1.0, 'mfp',
                                                         20.0, 3.0, 1.0, 'cm')
        self.assertEquals(det.comment(), "Ring detector tally 'ring' of "
                "neutrons: along x axis. ring x = 10 cm, radius 2 cm, s.o.e. "
                "radius 1 mfp; ring x = 20 cm, radius 3 "
                "cm, s.o.e. radius 1 cm; direct contrib is separate.")
        # mcnp()
        self.sim.add_tally(det)
        self.assertEquals(det.mcnp('%.5g', self.sim),
                "F35X:N  10 2 -1\n     20 3 1")
        det = cards.RingDetector('ring', 'neutron', 'x', 10.0, 2.0, 1.0, 'mfp',
                sep_direct=False)
        self.assertEquals(det.comment(), "Ring detector tally 'ring' of "
                "neutrons: along x axis. ring x = 10 cm, radius 2 cm, s.o.e. "
                "radius 1 mfp; direct contrib is not separate.")
        # mcnp()
        self.sim.add_tally(det)
        self.assertEquals(det.mcnp('%.5g', self.sim),
                "F35X:N  10 2 -1 ND")

    def test_EnergyGrid(self):
        """Tests :py:class:`cards.EnergyGrid`."""

        det = cards.RingDetector('ring', 'neutron', 'x', 10.0, 2.0, 1.0, 'mfp',
                sep_direct=False)
        ## EnergyGrid
        egrid = cards.EnergyGrid('grid0', None, [1e-4, 1, 100e3, 10e6])
        self.assertEquals(egrid.comment(), "Energy grid 'grid0' for "
                "all tallies: 4 groups.")
        self.assertEquals(egrid.mcnp('%.5g', self.sim),
                "E0 0.0001 1 1e+05 1e+07")
        egrid = cards.EnergyGrid('grid0', None, np.array([1e-4, 1, 100e3,
                    10e6]))
        self.assertEquals(egrid.comment(), "Energy grid 'grid0' for "
                "all tallies: 4 groups.")
        self.assertEquals(egrid.mcnp('%.5g', self.sim),
                "E0 0.0001 1 1e+05 1e+07")
        # For a specific tally.
        egrid = cards.EnergyGrid('grid1', det, np.array([1, 2]))
        self.sim.add_tally(det)
        self.assertEquals(egrid.comment(), "Energy grid 'grid1' for "
                "tally ring: 2 groups.")
        self.assertEquals(egrid.mcnp('%.5g', self.sim),
                "E1 1 2")


class TestMisc(unittest.TestCase):
    """Tests subclasses of :py:class:`cards.IMisc`."""
    pass


class TestSimulationDefinition(unittest.TestCase):
    """Tests the :py:class:`definition.SimulationDefinition` class."""
    # The system definition is complete.
    
    #sim = definition.MCNPSimulation(rxr)
    #sim.add_card(cards.Criticality())
    #sim.add_card(cards.CriticalityPoints())
    pass

class TestMCNPInput(unittest.TestCase):

    def setUp(self):
        self.sys = definition.SystemDefinition(verbose=False)
        self.sim = definition.MCNPSimulation(self.sys, verbose=False)

    def tearDown(self):
        #os.unlink('inptest')
        pass

    def test_InfLattice(self):
        """Tests the input file for an infinite lattice reactor. Checks
        generated output against the text file `inflattice_compare`.

        """
        fname = 'simplesim_inflattice'

        # Define system.
        # Materials.
        uo2 = cards.Material(name='UO2')
        uo2.from_atom_frac({'U235': 0.05,
                            'U238': 0.95,
                            'O16' : 2.00})
        h2o = cards.Material(name='H2O')
        h2o.from_atom_frac({'H1' : 2.0,
                            'O16': 1.0})
        # Surfaces.
        radius = 0.40
        pin = cards.AxisCylinder('pin', 'X', radius)
        pitch = 1.2
        cellbound = cards.Parallelepiped('bound',
                -pitch / 2, pitch / 2, -pitch / 2, pitch / 2, 0, 0,
                reflecting=True)
        # Cells.
        fuel = cards.CellMCNP('fuel', pin.neg, uo2,
                11.0, 'g/cm^3',
                importance=('neutron', 1),
                volume=1)
        coolant = cards.CellMCNP('coolant', pin.pos & cellbound.neg, h2o,
                1.0, 'g/cm^3',
                importance=('neutron', 1),
                volume=1)
        graveyard = cards.CellMCNP('graveyard', cellbound.pos,
                importance=('neutron', 0))

        # Add cells to the system.
        self.sys.add_cell(fuel)
        self.sys.add_cell(coolant)
        self.sys.add_cell(graveyard)

        # Add source, tallies to simulation.
        self.sim.add_misc(cards.ScatteringLaw('H2O', {'H1': 'lwtr'}))
        self.sim.add_source(cards.Criticality())
        self.sim.add_source(cards.CriticalityPoints())
        self.sim.add_tally(cards.CellFlux('flux', 'neutron', [fuel, coolant]))
        self.sim.add_misc(cards.EnergyGrid('egrid0', None,
                10**np.arange(-9.9, 1.1, 0.1)))

        # Create input file.
        inp = inputfile.MCNPInput(self.sim, heading="Infinite lattice.")
        inp.write(fname + '_default')
        # Check against text file.
        self.assertEquals(
                open(fname + '_default').readlines(),
                open(fname + '_default_compare').readlines())
        os.unlink(fname + '_default')

        # Test the & line continuation.
        inp = inputfile.MCNPInput(self.sim, heading="Infinite lattice.",
                cont_by_amp=True)
        inp.write(fname + '_amp')
        self.assertEquals(
                open(fname + '_amp').readlines(),
                open(fname + '_amp_compare').readlines())
        os.unlink(fname + '_amp')

    #@mock.patch_object(warnings, 'warn')
    def test_bypass_wrap(self): #, mock_warn):
        """Test that a warning is raised when a card requesting wrapping to be
        bypassed ever violates the 80-column rule.

        """
        # TODO expecting a warning.
        fname = 'simplesim_bypass_wrap'
        self.sim.add_source(cards.CriticalityPoints([[np.pi, np.pi, 0]]))
        inp = inputfile.MCNPInput(self.sim, float_format="%.50e")
        inp.write(fname)
        #self.assertTrue(mock_warn.called)
        os.unlink(fname)


class TestModifying(unittest.TestCase):
    """Ensures that the code is amenable to modifications made after initial
    creation.
    
    """
    pass

## The following tests are redundant, but are to make sure that the examples in
# the documentation function as expected.
class TestGodiva(unittest.TestCase):
    pass
