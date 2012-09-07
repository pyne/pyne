"""PyNE Simple Simulation Input Definitions tests

.. moduleauthor:: Chris Dembia <cld72@cornell.edu>

"""

import os
import pickle
import unittest

import numpy as np

from pyne import material
from pyne.simplesim import definition, cards, inputfile
import pyne.simplesim.nestedgeom as ng

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
        uo2 = material.from_atom_frac({'U235': 0.05, 'U238': 0.95, 'O16' : 2.00})
        self.uo2 = cards.Material(uo2, name='UO2')
        h2o = material.from_atom_frac({'H1' : 2.0, 'O16': 1.0}, attrs={'name': 'H2O'})
        self.h2o = cards.Material(h2o)
        
        # Surfaces.
        radius = 0.40 # cm
        self.pin = cards.AxisCylinder('fuelpin', 'Z', radius)
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

    def test_saving(self):
        """Tests saving the simulation."""
        fname = 'simplesim_saving'
        fid = open(fname, 'w')
        #pickle.dump(self.rxr, fid)
        fid.close()
        os.unlink(fname)

    def test_CardNumReturn(self):
        """Tests that a card number is returned when adding a card to a
        definition.
        
        """
        play1 = cards.Cell('play1', self.pin.neg)
        self.assertEquals(self.rxr.add_cell(play1), 4)
        
        radius = 0.40 # cm
        surf1 = cards.AxisCylinder('1', 'Z', radius)
        self.assertEquals(self.rxr.add_surface(surf1), 3)
        
        mat1 = cards.Material(material.Material(), name='H2O2')
        self.assertEquals(self.rxr.add_material(mat1), 3)

        sd = cards.Distribution('distA', [-2, 2], [0, 1])
        self.assertEquals(self.sim.add_dist(sd), 1)
        
        tr = cards.Transformation('source', [1, 0, 0], np.eye(3))
        self.assertEquals(self.sim.add_transformation(tr), 1)


    def test_SystemManipulation(self):
        """Tests the add and remove methods of
        :py:class:`pyne.simplesim.definition.SystemDefinition.
        
        """
        play1 = cards.Cell('play1', self.pin.neg)
        play2 = cards.Cell('play2', self.pin.neg)
        play3 = cards.Cell('play3', self.pin.neg)
        # Just making sure an exception is NOT raised.
        self.rxr.add_cell(play1, play2, play3)
        ncells = len(self.rxr.cells)

        radius = 0.40 # cm
        surf1 = cards.AxisCylinder('1', 'Z', radius)
        surf2 = cards.AxisCylinder('2', 'Z', radius)
        self.rxr.add_surface(surf1, surf2)
        nsurfs = len(self.rxr.surfaces)

        mat1 = cards.Material(material.Material(), name='H2O2')
        mat2 = cards.Material(material.Material(), name='H2O3')
        self.rxr.add_material(mat1, mat2)
        nmats = len(self.rxr.materials)

        self.rxr._register_universe('harhar')

        # Removal. TODO this check isn't thorough at all.
        self.rxr.remove_cell('play1')
        self.assertEquals(len(self.rxr.cells), ncells - 1)

        self.rxr.remove_surface('1')
        self.assertEquals(len(self.rxr.surfaces), nsurfs - 1)

        self.rxr.remove_material('H2O2')
        self.assertEquals(len(self.rxr.materials), nmats - 1)

        self.assertEquals(self.rxr.remove_universe('harhar'), 'harhar')
        self.assertEquals(len(self.rxr.universes), 0)

    def test_SimManipulation(self):
        """Tests the add and remove methods of
        :py:class:`pyne.simplesim.definition.SimulationDefinition.
        
        """
        cs = cards.Criticality()
        self.sim.add_source(cs)
        self.assertEquals(self.sim.remove_source('criticality'), cs)
        self.assertEquals(len(self.sim.source), 0)

        egrid = cards.EnergyGrid('g1', None, [0, 1])
        self.sim.add_misc(egrid)
        self.assertEquals(self.sim.remove_misc('g1'), egrid)
        self.assertEquals(len(self.sim.misc), 0)

        sd = cards.Distribution('distA', [-2, 2], [0, 1])
        self.sim.add_dist(sd)
        self.assertEquals(self.sim.remove_dist('distA'), sd)
        self.assertEquals(len(self.sim.dists), 0)
        
        tr = cards.Transformation('source', [1, 0, 0], np.eye(3))
        self.sim.add_transformation(tr)
        self.assertEquals(self.sim.remove_transformation('source'), tr)
        self.assertEquals(len(self.sim.transformations), 0)

        tally = cards.CellFlux('fuel', 'neutron', 'fuel')
        self.sim.add_tally(tally)

        self.assertEquals(self.sim.remove_tally('fuel'), tally)
        self.assertEquals(len(self.sim.tally), 0)

    def test_nestedgeom(self):
        """Tests the :py:mod:`pyne.simplesim.nestedgeom` module."""

        unit = ng.Surf('fuelpin')
        self.assertEquals(unit.comment(), " surf 'fuelpin'")
        self.assertEquals(unit.mcnp('', self.sim), " 1")
        
        unit = ng.Cell('coolant')
        self.assertEquals(unit.comment(), " cell 'coolant'")
        self.assertEquals(unit.mcnp('', self.sim), " 2")

        unit = ng.Univ('ha')
        self.sim.sys._register_universe('ha')
        self.assertEquals(unit.comment(), " univ 'ha'")
        self.assertEquals(unit.mcnp('', self.sim), " U=1")
        
        unit = ng.Union(ng.Surf('fuelpin'), ng.Surf('bound'))
        self.assertEquals(unit.comment(),
                " union of ( surf 'fuelpin', surf 'bound')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 2)")

        unit = ng.Union(ng.Cell('fuel'), ng.Cell('coolant'),
                ng.Cell('graveyard'))
        self.assertEquals(unit.comment(), 
                " union of ( cell 'fuel', cell 'coolant', cell 'graveyard')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 2 3)")

        unit = ng.Union(ng.Univ('ha'))
        self.assertEquals(unit.comment(), " union of ( univ 'ha')")
        self.assertEquals(unit.mcnp('', self.sim), " ( U=1)")

        # nesting
        unit = ng.Cell('fuel') < ng.FCell('coolant')
        self.assertEquals(unit.comment(), " ( cell 'fuel' in cell 'coolant')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < 2)")

        unit = ng.Cell('fuel').of(ng.FCell('coolant'))
        self.assertEquals(unit.comment(), " ( cell 'fuel' in cell 'coolant')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < 2)")

        unit = ng.Cell('fuel') < ng.FCell('coolant') < ng.FCell('graveyard')
        self.assertEquals(unit.comment(),
                " ( cell 'fuel' in cell 'coolant' in cell 'graveyard')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < 2 < 3)")

        unit = ng.Cell('fuel').of(ng.FCell('coolant').of(ng.FCell('graveyard')))
        self.assertEquals(unit.comment(),
                " ( cell 'fuel' in cell 'coolant' in cell 'graveyard')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < 2 < 3)")

        # unions and nesting, various syntax.
        unit = ng.Surf('fuelpin').of(ng.Union(ng.FCell('fuel'),
                                              ng.FCell('coolant')))
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant'))")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < ( 1 2))")

        unit = ng.Surf('fuelpin').of(ng.FCell('fuel').union(
            ng.FCell('coolant')))
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant'))")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < ( 1 2))")

        unit = ng.Surf('fuelpin').of(ng.FCell('fuel') | ng.FCell('coolant'))
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant'))")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < ( 1 2))")

        unit = ng.Surf('fuelpin').of(
                ng.FCell('fuel').union(ng.FCell('coolant')).of(
                    ng.FCell('fuel')))
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant') in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < ( 1 2) < 1)")

        unit = ng.Surf('fuelpin').of((ng.FCell('fuel') | 
            ng.FCell('coolant')).of(
            ng.FCell('fuel')))
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant') in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < ( 1 2) < 1)")

        unit = ng.Surf('fuelpin') < (ng.FCell('fuel') | ng.FCell('coolant')) < \
                ng.FCell('fuel')
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant') in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < ( 1 2) < 1)")

        # nesting with universe, and union
        unit = ng.Surf('fuelpin') < ng.Univ('ha') < ng.FCell('graveyard')
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in univ 'ha' in cell 'graveyard')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < U=1 < 3)")

        unit = ng.Surf('fuelpin') < ng.Union(ng.Univ('ha')) < \
                ng.FCell('graveyard')
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( univ 'ha') "
                "in cell 'graveyard')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < ( U=1) < 3)")

        # vectorized
        unit = ng.Vec(ng.Surf('fuelpin'), ng.Surf('bound')) < \
                ng.FCell('fuel')
        self.assertEquals(unit.comment(), 
                " ( over ( surf 'fuelpin', surf 'bound') in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 2 < 1)")

        usedlater = ng.Surf('fuelpin') < ng.Vec(ng.FCell('fuel'), 
                ng.FCell('coolant'))
        self.assertEquals(usedlater.comment(), 
                " ( surf 'fuelpin' in over ( cell 'fuel', cell 'coolant'))")
        self.assertEquals(usedlater.mcnp('', self.sim), " ( 1 < 1 2)")

        # using overloaded operator
        unit = (ng.Surf('fuelpin') & ng.Surf('bound')) < \
                ng.FCell('fuel')
        self.assertEquals(unit.comment(), 
                " ( over ( surf 'fuelpin', surf 'bound') in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 2 < 1)")

        unit = ng.Surf('fuelpin') < (ng.FCell('fuel') & ng.FCell('coolant'))
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in over ( cell 'fuel', cell 'coolant'))")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < 1 2)")

        # lattice
        unit = ng.Surf('fuelpin') < (ng.FCell('fuel') | ng.FCell('coolant',
            ng.Lin(5))) < ng.FCell('fuel')
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant'-lat linear idx 5) in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < ( 1 2[5]) < 1)")

        unit = ng.Surf('fuelpin') < (ng.FCell('fuel') |
                ng.FCell('coolant').lat(ng.Lin(5))) < ng.FCell('fuel')
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant'-lat linear idx 5) in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), " ( 1 < ( 1 2[5]) < 1)")

        unit = ng.Surf('fuelpin') < (ng.FCell('fuel') | ng.FCell('coolant',
            ng.Rng([0,1], [3,5], np.array([-1,4])))) < ng.FCell('fuel')
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant'-lat x range 0:1, y range 3:5, z range -1:4) "
                "in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), 
                " ( 1 < ( 1 2[0:1 3:5 -1:4]) < 1)")

        unit = ng.Surf('fuelpin') < (ng.FCell('fuel') | ng.FCell('coolant',
            ng.Cor([1, 3, 2]))) < ng.FCell('fuel')
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant'-lat coords (1, 3, 2)) "
                "in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), 
                " ( 1 < ( 1 2[ 1 3 2]) < 1)")

        unit = ng.Surf('fuelpin') < (ng.FCell('fuel') | ng.FCell('coolant',
            ng.Cor([np.array([1, 3, 2])]))) < ng.FCell('fuel')
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant'-lat coords (1, 3, 2)) "
                "in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), 
                " ( 1 < ( 1 2[ 1 3 2]) < 1)")

        unit = ng.Surf('fuelpin') < (ng.FCell('fuel') | ng.FCell('coolant',
            ng.Cor([[1, 3, 2], [-1, -2, -3]]))) < ng.FCell('fuel')
        self.assertEquals(unit.comment(), 
                " ( surf 'fuelpin' in union of ( cell 'fuel', cell "
                "'coolant'-lat coords (1, 3, 2), (-1, -2, -3)) "
                "in cell 'fuel')")
        self.assertEquals(unit.mcnp('', self.sim), 
                " ( 1 < ( 1 2[ 1 3 2, -1 -2 -3]) < 1)")

        tally = cards.SurfaceFlux('fuel', 'neutron', unit)
        self.sim.add_tally(tally)
        self.assertEquals(tally.comment(), "Surface flux tally 'fuel' "
                "of neutrons:  ( surf 'fuelpin' in union of ( cell "
                "'fuel', cell 'coolant'-lat coords (1, 3, 2), (-1, -2, -3)"
                ") in cell 'fuel').")
        self.assertEquals(tally.mcnp('%.5g', self.sim), 
                "F12:N  ( 1 < ( 1 2[ 1 3 2, -1 -2 -3]) < 1)")

        tally = cards.SurfaceFlux('fuel', 'neutron', ['bound', unit])
        self.assertEquals(tally.comment(), "Surface flux tally 'fuel' "
                "of neutrons: surfaces 'bound';  ( surf 'fuelpin' in "
                "union of ( cell 'fuel', cell 'coolant'-lat coords "
                "(1, 3, 2), (-1, -2, -3)) in cell 'fuel').")
        self.assertEquals(tally.mcnp('%.5g', self.sim), 
                "F12:N  2 ( 1 < ( 1 2[ 1 3 2, -1 -2 -3]) < 1)")

        tally = cards.SurfaceFlux('fuel', 'neutron', [unit, usedlater])
        self.assertEquals(tally.comment(), "Surface flux tally 'fuel' "
                "of neutrons: surfaces  ( surf 'fuelpin' in "
                "union of ( cell 'fuel', cell 'coolant'-lat coords "
                "(1, 3, 2), (-1, -2, -3)) in cell 'fuel');  "
                "( surf 'fuelpin' in over ( cell 'fuel', cell 'coolant')).")
        self.assertEquals(tally.mcnp('%.5g', self.sim), 
                "F12:N  ( 1 < ( 1 2[ 1 3 2, -1 -2 -3]) < 1) ( 1 < 1 2)")


    def test_Region(self):
        """Tests the shifting and stretching of surfaces in a Region."""

        reg1 = self.pin.pos & self.cellbound.neg
        reg1.stretch([2, 2, 4])
        self.assertEquals(self.pin.radius, 0.8)
        self.assertEquals(self.cellbound.xlims.tolist(), [-1.2, 1.2])
        self.assertEquals(self.cellbound.ylims.tolist(), [-1.2, 1.2])
        self.assertEquals(self.cellbound.zlims.tolist(), [0, 0])

#    def test_Facet(self):
#        """Tests the usage of macrobody :py:class:`cards.Facet`'s."""
#        myreg = self.cellbound.facet('east').neg
#        cell = cards.Cell('facettest', myreg)
#        self.rxr.add_cell(cell)
#        self.assertEquals(cell.comment(),
#                "Cell 'facettest': region -bound, void.")
#        self.assertEquals(cell.mcnp('%.5g', self.sim), "4 0 -2.1")

    def test_Universes(self):
        """Tests :py:class:`cards.Universes`."""
        uni = cards.Universes('fuel', 'unitcell', True, 
                            'coolant', 'unitcell', True,
                            'graveyard', 'grave', False)
        self.sim.add_misc(uni)
        self.assertEquals(uni.comment(), "Universes 'universes': "
                "cell 'fuel' unitcell (truncated), "
                "cell 'coolant' unitcell (truncated), "
                "cell 'graveyard' grave (not truncated).")
        self.assertEquals(uni.mcnp('%.5g', self.sim), 
                "U 1 1 -2")

        uni = cards.Universes()
        uni.set('fuel', 'unitcell', True)
        uni.set('coolant', 'unitcell', True)
        self.sim.add_misc(uni)
        self.assertEquals(uni.comment(), "Universes 'universes': "
                "cell 'fuel' unitcell (truncated), "
                "cell 'coolant' unitcell (truncated).")
        self.assertEquals(uni.mcnp('%.5g', self.sim), 
                "U 1 1 1J")

    def test_Fill(self):
        """Tests :py:class:`cards.Fill`."""
        uni = cards.Universes()
        uni.set('fuel', 'unitcell', True)
        uni.set('coolant', 'otheruniv', True)
        self.sim.add_misc(uni) 

        uni = cards.Fill('fuel', 'unitcell',
                         'coolant', 'otheruniv')
        self.assertEquals(uni.comment(), "Fill 'fill': "
                "cell 'fuel' unitcell, "
                "cell 'coolant' otheruniv.")
        self.assertEquals(uni.mcnp('%.5g', self.sim), 
                "FILL 1 2 1J")

        uni = cards.Fill()
        uni.set('fuel', 'unitcell')
        uni.set('coolant', 'otheruniv')
        self.assertEquals(uni.comment(), "Fill 'fill': "
                "cell 'fuel' unitcell, "
                "cell 'coolant' otheruniv.")
        self.assertEquals(uni.mcnp('%.5g', self.sim), 
                "FILL 1 2 1J")

    def test_Lattice(self):
        """Tests :py:class:`Lattice`."""

        lat = cards.Lattice('fuel', 'hexahedra', 'coolant', 'hexagonal')
        self.assertEquals(lat.comment(), "Lattice 'lattice': "
                "cell 'fuel' hexahedra, cell 'coolant' hexagonal.")
        self.assertEquals(lat.mcnp('%.5g', self.sim), "LAT 1 2 1J")

        lat = cards.Lattice()
        lat.set('fuel', 'hexahedra')
        lat.set('coolant', 'hexagonal')
        self.assertEquals(lat.comment(), "Lattice 'lattice': "
                "cell 'fuel' hexahedra, cell 'coolant' hexagonal.")
        self.assertEquals(lat.mcnp('%.5g', self.sim), "LAT 1 2 1J")

    def test_Distribution(self):
        """Tests the :py:class:`cards.Distribution` class."""

        sd = cards.Distribution('distA', [-2, 2], [0, 1])
        self.sim.add_dist(sd)
        self.assertEquals(sd.comment(), "Source distribution 'distA': "
                "key setting default, val setting default, "
                "KEYS: -2, 2, VALS: 0, 1.")
        self.assertEquals(sd.mcnp('%.5g', self.sim), 
                "SI1 -2 2\nSP1 0 1")
        sd = cards.Distribution('distA', [-2, 2], [0, 1], 'histogram', 'prob')
        self.sim.add_dist(sd)
        self.assertEquals(sd.comment(), "Source distribution 'distA': "
                "key setting histogram, val setting prob, "
                "KEYS: -2, 2, VALS: 0, 1.")
        self.assertEquals(sd.mcnp('%.5g', self.sim), 
                "SI1 H -2 2\nSP1 D 0 1")

        sd = cards.Distribution('distA', 
                ['neutron', 'photon', 'spont-fiss-by-hist', 92238],
                [1, 1, 'fuel', 2],
                key_setting='discrete',
                val_setting='partint')
        self.sim.add_dist(sd)
        self.assertEquals(sd.comment(), "Source distribution 'distA': "
                "key setting discrete, val setting partint, "
                "KEYS: neutron, photon, spont-fiss-by-hist, 92238, "
                "VALS: 1, 1, fuel, 2.")
        self.assertEquals(sd.mcnp('%.5g', self.sim),
                "SI1 L N P -SF 92238\nSP1 W 1 1 -1 2")

        sd = cards.Distribution('distA', ['distB', 'distC'], [0.3, 0.7],
                'dist')

        sdB = cards.Distribution('distB', [-2, 2], [0, 1])
        sdC = cards.Distribution('distC', [-2, 2], [0, 1])
        self.sim.add_dist(sd)
        self.sim.add_dist(sdB)
        self.sim.add_dist(sdC)
        self.assertEquals(sd.comment(), "Source distribution 'distA': "
                "key setting dist, val setting default, "
                "KEYS: distB, distC, "
                "VALS: 0.3, 0.7.")
        self.assertEquals(sd.mcnp('%.5g', self.sim),
                "SI1 S 2 3\nSP1 0.3 0.7")

        sd = cards.Distribution('distA', [], [], 'analytic', 'maxwell')
        self.sim.add_dist(sd)
        self.assertEquals(sd.comment(), "Source distribution 'distA': "
                "key setting analytic, val setting maxwell, "
                "KEYS: default, VALS: default.")
        self.assertEquals(sd.mcnp('%.5g', self.sim), "SP1 -2")

        sd = cards.Distribution('distA', [], [], key_setting='analytic',
                                                 val_setting='maxwell')
        self.sim.add_dist(sd)
        self.assertEquals(sd.comment(), "Source distribution 'distA': "
                "key setting analytic, val setting maxwell, "
                "KEYS: default, VALS: default.")
        self.assertEquals(sd.mcnp('%.5g', self.sim),
                "SP1 -2")

        sd = cards.Distribution('distA', [0, 10], [1, 3], 'analytic', 'watt')
        self.sim.add_dist(sd)
        self.assertEquals(sd.comment(), "Source distribution 'distA': "
                "key setting analytic, val setting watt, "
                "KEYS: 0, 10, VALS: 1, 3.")
        self.assertEquals(sd.mcnp('%.5g', self.sim),
                "SI1 0 10\nSP1 -3 1 3")

    def test_GeneralSource(self):
        """Tests the :py:class:`cards.GeneralSource` class."""

        gs = cards.GeneralSource()
        self.assertEquals(gs.name, 'gensource')
        self.assertEquals(gs.comment(), "General source 'gensource'.")
        self.assertEquals(gs.mcnp('%.5g', self.sim), "SDEF")
        gs = cards.GeneralSource(particle='photon', cell='fuel')
        self.assertEquals(gs.comment(), "General source 'gensource': "
                "particle=photon, cell=fuel.")
        self.assertEquals(gs.mcnp('%.5g', self.sim), "SDEF PAR=P CEL=1")

        gs = cards.GeneralSource(particle='partdist', cell='celldist')
        sdA = cards.Distribution('partdist', [-2, 2], [0, 1])
        sdB = cards.Distribution('celldist', [-2, 2], [0, 1])
        self.sim.add_dist(sdA)
        self.sim.add_dist(sdB)
        self.assertEquals(gs.comment(), "General source 'gensource': "
                "particle=partdist, cell=celldist.")
        self.assertEquals(gs.mcnp('%.5g', self.sim), "SDEF PAR=D1 CEL=D2")

        gs = cards.GeneralSource(energy='energydist', time=2e10)
        sdA = cards.Distribution('energydist', [-2, 2], [0, 1])
        self.sim.add_dist(sdA)
        self.assertEquals(gs.comment(), "General source 'gensource': "
                "energy=energydist, time=20000000000.0.")
        self.assertEquals(gs.mcnp('%.5g', self.sim), "SDEF ERG=D3 TME=2e+18")

        gs = cards.GeneralSource(ref_dir=[1, 0, 0], cosine=0.5)
        self.assertEquals(gs.comment(), "General source 'gensource': "
                "ref. dir=[1, 0, 0], cosine=0.5.")
        self.assertEquals(gs.mcnp('%.5g', self.sim), "SDEF VEC=1 0 0 DIR=0.5")

        gs = cards.GeneralSource(x='xdist', y='ydist', z='zdist')
        sdA = cards.Distribution('xdist', [-2, 2], [0, 1])
        sdB = cards.Distribution('ydist', [-2, 2], [0, 1])
        sdC = cards.Distribution('zdist', [-2, 2], [0, 1])
        self.sim.add_dist(sdA)
        self.sim.add_dist(sdB)
        self.sim.add_dist(sdC)
        self.assertEquals(gs.comment(), "General source 'gensource': "
                "x=xdist, y=ydist, z=zdist.")
        self.assertEquals(gs.mcnp('%.5g', self.sim), "SDEF X=D4 Y=D5 Z=D6")
        
        tr = cards.Transformation('source', [1, 0, 0], np.eye(3))
        self.sim.add_transformation(tr)
        gs = cards.GeneralSource(transformation='source')
        self.assertEquals(gs.comment(), "General source 'gensource': "
                "trans.=source.")
        self.assertEquals(gs.mcnp('%.5g', self.sim), "SDEF TR=1")

        gs = cards.GeneralSource(beam_emit=(1, 2, 3), beam_aper=(4, 5))
        self.assertEquals(gs.comment(), "General source 'gensource': "
                "beam emit.=(1, 2, 3), beam aper.=(4, 5).")
        self.assertEquals(gs.mcnp('%.5g', self.sim), 
                "SDEF BEM=1 2 3 BAP=4 5")

    def test_Saving(self):
        """Tests saving a system definition to a JSON file."""
        self.sim.save('simplesim_sys.json')
        os.remove('simplesim_sys.json')


    def test_Material(self):
        """Tests :py:class:`pyne.simplesim.cards.Material`."""

        originstory = "I found this water in a well a few years ago."
        h2o = cards.Material(material.from_atom_frac({10010: 1.0, 'O16': 2.0}), 
                             name='water', description=originstory)
        h2o.tables = {'10010': '71c'}
        self.sim.sys.add_material(h2o)
        self.assertEquals(h2o.comment(), "Material 'water': "
                "I found this water in a well a few years ago.")
        self.assertEquals(h2o.mcnp('%g', self.sim), "M3\n"
        "       1001.71c  1 $ H1\n"
        "       8016      2 $ O16")

        h2o = cards.Material(material.from_atom_frac({10010: 1.0, 'O16': 2.0}),
                             name='water', tables={'10010': '71c'})
        self.assertEquals(h2o.mcnp('%g', self.sim), "M3\n"
        "       1001.71c  1 $ H1\n"
        "       8016      2 $ O16")

        h2o = cards.Material(material.from_atom_frac({10010: 1.0, 'O16': 2.0}),
                             name='water', tables={'H1': '71c'})
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
        self.assertEquals(cellA.mcnp("%.5e", self.sim), "4 0 (-1:2)")
        cellB = cards.Cell('B', self.pin.neg & self.cellbound.pos, self.uo2,
                10.0, 'g/cm^3')
        self.rxr.add_cell(cellB)
        self.assertEquals(cellB.comment(), "Cell 'B': region "
                "(-fuelpin & +bound), material 'UO2' density 10 "
                "g/cm^3.")
        self.assertEquals(cellB.mcnp("%.5e", self.sim), "5 1 -1.00000e+01 "
                "(-1:2)")
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

        # Universes, fill, lattice.
        cellC = cards.Cell('C', self.pin.neg & self.cellbound.pos, self.uo2,
                density=1,
                density_units='g/cm^3',
                universe=('unitcell', True))
        self.rxr.add_cell(cellC)
        self.assertEquals(cellC.comment(), "Cell 'C': region "
                "(-fuelpin & +bound), "
                "material 'UO2' density 1 g/cm^3 "
                "univ unitcell (truncated).")
        self.assertEquals(cellC.mcnp('%.5g', self.sim), "6 1 -1 (-1:2) U=1")

        cellD = cards.Cell('D', self.pin.neg, self.uo2, density=2,
                density_units='g/cm^3',
                universe=('unitcell', False))
        self.rxr.add_cell(cellD)
        self.assertEquals(cellD.comment(), "Cell 'D': region "
                "-fuelpin, "
                "material 'UO2' density 2 g/cm^3 "
                "univ unitcell (not truncated).")
        self.assertEquals(cellD.mcnp('%.5g', self.sim), "7 1 -2 -1 U=-1")

        cellE = cards.Cell('E', self.pin.neg, fill='unitcell')
        self.rxr.add_cell(cellE)
        self.assertEquals(cellE.comment(), "Cell 'E': region "
                "-fuelpin, void "
                "fill unitcell.")
        self.assertEquals(cellE.mcnp('%.5g', self.sim), "8 0 -1 FILL=1")

        univ_names = [['A', 'B', 'A', 'B', 'A'],
                      ['B', 'A', 'B', 'A', 'B'],
                      ['B', 'A', 'B', 'A', 'B']]
        cellF = cards.Cell('F', self.pin.neg, lattice='hexahedra',
                fill=([0, 4], [0, 2], [], univ_names))
        self.rxr.add_cell(cellF)
        self.rxr._register_universe('A')
        self.rxr._register_universe('B')
        self.assertEquals(cellF.comment(), "Cell 'F': region "
                "-fuelpin, void lattice hexahedra fill long format.")
        self.assertEquals(cellF.mcnp('%.5g', self.sim), "9 0 -1 LAT=1\n     "
                "FILL=0:4 0:2 0:0\n     "
                "2 3 2 3 2\n     "
                "3 2 3 2 3\n     "
                "3 2 3 2 3")

        univ_names2 = [[['A',  'B',  'A'], 
                        ['B',  'A',  'B']],
                       [['A',  None, 'A'],
                        ['B',  'A',  'B']]]
        cellF = cards.Cell('F', self.pin.neg, lattice='hexagonal',
                fill=([0, 2], [0, 1], [-1, 0], univ_names2))
        self.rxr.add_cell(cellF)
        self.assertEquals(cellF.comment(), "Cell 'F': region "
                "-fuelpin, void lattice hexagonal fill long format.")
        self.assertEquals(cellF.mcnp('%.5g', self.sim), "9 0 -1 LAT=2\n     "
                "FILL=0:2 0:1 -1:0\n     "
                "$ k = -1\n     "
                "2 3 2\n     "
                "3 2 3\n     "
                "$ k = 0\n     "
                "2 0 2\n     "
                "3 2 3")

        univ_names2 = ['A',  'B',  'A']
        cellF = cards.Cell('F', self.pin.neg, lattice='hexagonal',
                fill=([0, 2], [], [], univ_names2))
        self.rxr.add_cell(cellF)
        self.assertEquals(cellF.comment(), "Cell 'F': region "
                "-fuelpin, void lattice hexagonal fill long format.")
        self.assertEquals(cellF.mcnp('%.5g', self.sim), "9 0 -1 LAT=2\n     "
                "FILL=0:2 0:0 0:0 2 3 2")

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
        self.assertEquals(cellC.mcnp("%.5e", self.sim), "4 0 (-1:2)")
        cellD = cards.Cell('D', self.pin.neg | self.cellbound.pos, self.uo2,
                10.0, 'g/cm^3')
        self.rxr.add_cell(cellD)
        self.assertEquals(cellD.comment(), "Cell 'D': region "
                "(-fuelpin | +bound), material 'UO2' density 10 "
                "g/cm^3.")
        self.assertEquals(cellD.mcnp("%.5e", self.sim), "5 1 -1.00000e+01 "
                "(-1 2)")

        # Keywords.
        # Temperature, volume, importance.
        cellE = cards.CellMCNP('E', self.pin.neg | self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                temperature=600, volume=1,
                importance=('neutron', 1))
        self.rxr.add_cell(cellE)
        self.assertEquals(cellE.comment(), "Cell 'E': region "
                "(-fuelpin | +bound), material 'UO2' density 10 "
                "g/cm^3 TMP= 600 K VOL= 1 cm^3 IMP:N= 1.")
        self.assertEquals(cellE.mcnp('%g', self.sim), "6 1 -10 "
                "(-1 2) TMP=5.17041e-08 VOL=1 IMP:N=1")
        cellF = cards.CellMCNP('F', self.pin.neg | self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                importance=[('neutron', 1), ('photon', 0)])
        self.rxr.add_cell(cellF)
        self.assertEquals(cellF.comment(), "Cell 'F': region "
                "(-fuelpin | +bound), material 'UO2' density 10 "
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
                "(-1:2) EXT:N=S EXT:P=-0.5X")
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
                "(-1:2) EXT:N=-0.5V0")
        # Forced collisions
        cellI = cards.CellMCNP('I', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                forced_coll=('neutron', 0.5, False))
        self.rxr.add_cell(cellI)
        self.assertEquals(cellI.comment(), "Cell 'I': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "FCL:N= prob 0.5 for entering and weight games.")
        self.assertEquals(cellI.mcnp('%g', self.sim), "10 1 -10 "
                "(-1:2) FCL:N=0.5")
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
                "(-1:2) WWN3:N=-1")
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
                "(-1:2) DXC1:N=0.5")
        # Photon weight
        cellL = cards.CellMCNP('L', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                photon_weight=('one',))
        self.rxr.add_cell(cellL)
        self.assertEquals(cellL.comment(), "Cell 'L': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "PWT= one.")
        self.assertEquals(cellL.mcnp('%g', self.sim), "13 1 -10 (-1:2) "
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
        self.assertEquals(cellL2.mcnp('%g', self.sim), "14 1 -10 (-1:2) "
                "PWT=0.5")
        cellM = cards.CellMCNP('M', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                photon_weight=(0.5, True))
        self.rxr.add_cell(cellM)
        self.assertEquals(cellM.comment(), "Cell 'M': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "PWT= 0.5 (pre-weighted).")
        self.assertEquals(cellM.mcnp('%g', self.sim), "15 1 -10 "
                "(-1:2) PWT=-0.5")
        # Fission turnoff
        cellN = cards.CellMCNP('N', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                fission_turnoff='capture-nogamma')
        self.rxr.add_cell(cellN)
        self.assertEquals(cellN.comment(), "Cell 'N': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "NONU= capture-nogamma.")
        self.assertEquals(cellN.mcnp('%g', self.sim), "16 1 -10 "
                "(-1:2) NONU=2")
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
                "(-1:2) PD15=0.5")
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
                "(-1:2) TRCL=1")
        cellP = cards.CellMCNP('P', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                transformation=([1, 0, 0], np.eye(3)))
        self.rxr.add_cell(cellP)
        self.assertEquals(cellP.comment(), "Cell 'P': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "TRCL aux origin in main (1, 0, 0) cm, x' <1, 0, 0>, "
                "y' <0, 1, 0>, z' <0, 0, 1>.")
        self.assertEquals(cellP.mcnp('%g', self.sim), "18 1 -10 "
                "(-1:2) TRCL ( 1 0 0 1 0 0 0 1 0 0 0 1 1)")
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
                "(-1:2) *TRCL ( 1 0 0 1 0 0 0 1 0 0 0 1 1)")
        cellP = cards.CellMCNP('P', self.pin.neg & self.cellbound.pos,
                self.uo2, 10.0, 'g/cm^3',
                user_custom='EXT:N 0.7V2')
        self.rxr.add_cell(cellP)
        self.assertEquals(cellP.comment(), "Cell 'P': region "
                "(-fuelpin & +bound), material 'UO2' density 10 g/cm^3 "
                "and user's custom input.")
        self.assertEquals(cellP.mcnp('%g', self.sim), "18 1 -10 "
                "(-1:2) EXT:N 0.7V2")
        #TODO set temperature to 100 or -1
        #cellE = cards.CellMCNP('E', self.pin.neg & self.cellbound.pos,
        #        self.uo2, 10.0, 'g/cm^3',
        #        temperature=100)
        #self.assertRaises(UserWarning, cellE.comment)

        # Universes, fill, lattice.
        cellC = cards.CellMCNP('C', self.pin.neg & self.cellbound.pos, self.uo2,
                density=1,
                density_units='g/cm^3',
                volume=1,
                universe=('unitcell', True))
        self.rxr.add_cell(cellC)
        self.assertEquals(cellC.comment(), "Cell 'C': region "
                "(-fuelpin & +bound), "
                "material 'UO2' density 1 g/cm^3 "
                "univ unitcell (truncated) VOL= 1 cm^3.")
        self.assertEquals(cellC.mcnp('%.5g', self.sim), "4  1 -1 (-1:2) "
                "U=1 VOL=1")

        cellD = cards.CellMCNP('D', self.pin.neg, self.uo2, density=2,
                density_units='g/cm^3',
                volume=1,
                universe=('unitcell', False))
        self.rxr.add_cell(cellD)
        self.assertEquals(cellD.comment(), "Cell 'D': region "
                "-fuelpin, "
                "material 'UO2' density 2 g/cm^3 "
                "univ unitcell (not truncated) VOL= 1 cm^3.")
        self.assertEquals(cellD.mcnp('%.5g', self.sim), "5  1 -2 -1 U=-1 VOL=1")

        cellE = cards.CellMCNP('E', self.pin.neg, fill='unitcell', volume=1)
        self.rxr.add_cell(cellE)
        self.assertEquals(cellE.comment(), "Cell 'E': region "
                "-fuelpin, void "
                "fill unitcell VOL= 1 cm^3.")
        self.assertEquals(cellE.mcnp('%.5g', self.sim), "6  0 -1 FILL=1 VOL=1")

        univ_names = [['A', 'B', 'A', 'B', 'A'],
                      ['B', 'A', 'B', 'A', 'B'],
                      ['B', 'A', 'B', 'A', 'B']]
        cellF = cards.CellMCNP('F', self.pin.neg, lattice='hexahedra',
                volume=1,
                fill=([0, 4], [0, 2], [], univ_names))
        self.rxr.add_cell(cellF)
        self.rxr._register_universe('A')
        self.rxr._register_universe('B')
        self.assertEquals(cellF.comment(), "Cell 'F': region "
                "-fuelpin, void lattice hexahedra fill long format "
                "VOL= 1 cm^3.")
        self.assertEquals(cellF.mcnp('%.5g', self.sim), "7  0 -1 LAT=1 VOL=1"
                "\n     "
                "FILL=0:4 0:2 0:0\n     "
                "2 3 2 3 2\n     "
                "3 2 3 2 3\n     "
                "3 2 3 2 3")

        univ_names2 = [[['A',  'B',  'A'], 
                        ['B',  'A',  'B']],
                       [['A',  None, 'A'],
                        ['B',  'A',  'B']]]
        cellF = cards.CellMCNP('F', self.pin.neg, lattice='hexagonal',
                volume=1,
                fill=([0, 2], [0, 1], [-1, 0], univ_names2))
        self.rxr.add_cell(cellF)
        self.assertEquals(cellF.comment(), "Cell 'F': region "
                "-fuelpin, void lattice hexagonal fill long format "
                "VOL= 1 cm^3.")
        self.assertEquals(cellF.mcnp('%.5g', self.sim), "7  0 -1 LAT=2 VOL=1"
                "\n     "
                "FILL=0:2 0:1 -1:0\n     "
                "$ k = -1\n     "
                "2 3 2\n     "
                "3 2 3\n     "
                "$ k = 0\n     "
                "2 0 2\n     "
                "3 2 3")

        univ_names2 = ['A',  'B',  'A']
        cellF = cards.CellMCNP('F', self.pin.neg, lattice='hexagonal',
                volume=1,
                fill=([0, 2], [], [], univ_names2))
        self.rxr.add_cell(cellF)
        self.assertEquals(cellF.comment(), "Cell 'F': region "
                "-fuelpin, void lattice hexagonal fill long format "
                "VOL= 1 cm^3.")
        self.assertEquals(cellF.mcnp('%.5g', self.sim), "7  0 -1 LAT=2 VOL=1"
                "\n     "
                "FILL=0:2 0:0 0:0 2 3 2")

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
        extn = cards.ExponentialTransform('neutron', 'fuel', 'capture-to-total',
                'currdir', 'toward')
        self.assertEquals(extn.name, 'exptransform-neutron')
        self.assertIs(extn.cells[0], 'fuel')
        self.assertEquals(extn.stretchs['fuel'], 'capture-to-total')
        self.assertEquals(extn.directions['fuel'], 'currdir')
        self.assertEquals(extn.signs['fuel'], 'toward')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "capture-to-total toward currdir.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), "EXT:N S 0 0")
        extn = cards.ExponentialTransform('neutron', 'coolant', 0.5,
                'currdir', 'toward')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'coolant' stretch by "
                "0.5 toward currdir.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), "EXT:N 0 5.0e-01 0")
        extn = cards.ExponentialTransform('neutron', 'fuel', 0.5, 'x',
                'toward')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "0.5 toward x.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), "EXT:N 5.0e-01X 0 0")
        extn = cards.ExponentialTransform('neutron', 'fuel', 0.5, 'origin',
                'away')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "0.5 away from origin.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), "EXT:N -5.0e-01V0 0 0")
        # Multiple cells.
        extn = cards.ExponentialTransform('neutron', 
                'fuel', 'capture-to-total', 'currdir', 'toward', 
                'coolant', 0.5, 'currdir', 'toward',
                'graveyard', 0.5, 'x-axis', 'away')
        self.assertEquals(extn.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "capture-to-total toward currdir; cell 'coolant' stretch by "
                "0.5 toward currdir; cell 'graveyard' stretch by 0.5 away "
                "from x-axis.")
        self.assertEquals(extn.mcnp('%.1e', self.sim), 
                "EXT:N S 5.0e-01 -5.0e-01V1")
        # Using set().
        ext2 = cards.ExponentialTransform('neutron')
        ext2.set('fuel', 'capture-to-total', 'currdir', 'toward')
        ext2.set('coolant', 0.5, 'currdir', 'toward')
        ext2.set('graveyard', 0.5, 'x-axis', 'away')
        self.assertEquals(ext2.comment(), "Exponential transform "
                "'exptransform-neutron': cell 'fuel' stretch by "
                "capture-to-total toward currdir; cell 'coolant' stretch by "
                "0.5 toward currdir; cell 'graveyard' stretch by 0.5 away "
                "from x-axis.")
        self.assertEquals(ext2.mcnp('%.1e', self.sim), 
                "EXT:N S 5.0e-01 -5.0e-01V1")
        # Modifying.
        ext2.set('graveyard', 0.7, 'x-axis', 'away')
        self.assertEquals(ext2.mcnp('%.1e', self.sim), 
                "EXT:N S 5.0e-01 -7.0e-01V1")
        # Manage order-of-cells issue/question.
        # TODO editing entries.

    def test_ForcedCollision(self):
        ## ForcedCollision
        fcl = cards.ForcedCollision('neutron', 'fuel', 0.5, True)
        self.assertEquals(fcl.name, 'forcedcoll-neutron')
        self.assertEquals(fcl.cells[0], 'fuel')
        self.assertEquals(fcl.probs['fuel'], 0.5)
        self.assertTrue(fcl.only_enterings['fuel'])
        self.assertEquals(fcl.comment(), "Forced collision "
                "'forcedcoll-neutron': cell 'fuel' prob 0.5 for entering only.")
        self.assertEquals(fcl.mcnp('%.1e', self.sim), "FCL:N -5.0e-01 0 0")
        fcl = cards.ForcedCollision('neutron', 'fuel', 0.5, False)
        self.assertEquals(fcl.comment(), "Forced collision "
                "'forcedcoll-neutron': cell 'fuel' prob 0.5 for entering "
                "and weight games.")
        self.assertEquals(fcl.mcnp('%.1e', self.sim), "FCL:N 5.0e-01 0 0")
        fcl = cards.ForcedCollision('neutron', 'fuel', 0.5, True,
                                               'coolant', 0.5, False)
        self.assertEquals(fcl.comment(), "Forced collision "
                "'forcedcoll-neutron': cell 'fuel' prob 0.5 for entering "
                "only; cell 'coolant' prob 0.5 for entering and weight games.")
        self.assertEquals(fcl.mcnp('%.1e', self.sim), "FCL:N -5.0e-01 5.0e-01 0")
        fcl = cards.ForcedCollision('neutron')
        fcl.set('fuel', 0.5, True)
        fcl.set('coolant', 0.5, False)
        self.assertEquals(fcl.comment(), "Forced collision "
                "'forcedcoll-neutron': cell 'fuel' prob 0.5 for entering "
                "only; cell 'coolant' prob 0.5 for entering and weight games.")
        self.assertEquals(fcl.mcnp('%.1e', self.sim), "FCL:N -5.0e-01 5.0e-01 0")
        fcl.set('coolant', 0.7, False)
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

        wwn = cards.WeightWindowBound('photon', 1, 1, 'D', 0.01)
        self.assertEquals(wwn.name, 'weightwinbound-photon')
        self.assertEquals(wwn.comment(), "Weight window bounds "
                "'weightwinbound-photon' for photons: energy idx 1: "
                "time idx 1: cell 'D': 0.01.")
        self.assertEquals(wwn.mcnp('%g', self.sim), "WWN1:P 0 0 0 0 0.01 0")
        # More than one cell.
        wwn = cards.WeightWindowBound('photon' , 1, 1, 'D', 0.01 , 1, 1,
                'E', 0.02 , 2, 1, 'E', 0.03 , 1, 3, 'D', 0.04 , 2, 3,
                'D', 0.05 , 2, 3, 'E', 0.06)
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
        wwn.set(1, 1, 'D', 0.01)
        wwn.set(1, 1, 'E', 0.02)
        wwn.set(2, 1, 'E', 0.03)
        wwn.set(1, 3, 'D', 0.04)
        wwn.set(2, 3, 'D', 0.05)
        wwn.set(2, 3, 'E', 0.06)
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
        wwn = cards.WeightWindowBound('photon', 8, 6, 'D', 0.01)
        self.assertEquals(wwn.comment(), "Weight window bounds "
                "'weightwinbound-photon' for photons: energy idx 8: "
                "time idx 6: cell 'D': 0.01.")
        self.assertEquals(wwn.mcnp('%g', self.sim), "WWN58:P 0 0 0 0 0.01 0")

        # Test killall
        wwn = cards.WeightWindowBound('photon', 1, 1, 'D', 'killall')
        self.assertEquals(wwn.comment(), "Weight window bounds "
                "'weightwinbound-photon' for photons: energy idx 1: "
                "time idx 1: cell 'D': killall.")
        self.assertEquals(wwn.mcnp('%g', self.sim), "WWN1:P 0 0 0 0 -1 0")

        # Check exception when using a particle type for which WWE/WWT cards
        # are not defined.
        wwn = cards.WeightWindowBound('neutron', 1, 1, 'D', 0.01)
        self.assertRaises(Exception, wwn.mcnp,'%g', self.sim)
        try:
            wwn.mcnp('%g', self.sim)
        except Exception as e:
            self.assertEquals(e.message, "No WWGT:N or WWT:N card found in "
                    "the simulation.")
        ## Importance
        impn = cards.Importance('neutron', 'fuel', 1, 'coolant', 2)
        self.assertEquals(impn.name, 'importance-neutron')
        self.assertEquals(impn.comment(), "Importance 'importance-neutron': "
                "cell 'fuel' 1; cell 'coolant' 2.")
        self.assertEquals(impn.mcnp('%g', self.sim), "IMP:N 1 2 0 0 0 0")
        impn.set('graveyard', 0)
        self.assertEquals(impn.comment(), "Importance 'importance-neutron': "
                "cell 'fuel' 1; cell 'coolant' 2; cell 'graveyard' 0.")
        self.assertEquals(impn.mcnp('%g', self.sim), "IMP:N 1 2 0 0 0 0")
        # Modifying.
        impn.set('graveyard', 1)
        self.assertEquals(impn.comment(), "Importance 'importance-neutron': "
                "cell 'fuel' 1; cell 'coolant' 2; cell 'graveyard' 1.")
        self.assertEquals(impn.mcnp('%g', self.sim), "IMP:N 1 2 1 0 0 0")
        # Using args.
        args = ['fuel', 1, 'coolant', 2, 'graveyard', 3]
        impn = cards.Importance('neutron', *args)
        self.assertEquals(impn.comment(), "Importance 'importance-neutron': "
                "cell 'fuel' 1; cell 'coolant' 2; cell 'graveyard' 3.")
        self.assertEquals(impn.mcnp('%g', self.sim), "IMP:N 1 2 3 0 0 0")

        ## Volume
        vol = cards.Volume('fuel', 1)
        self.assertEquals(vol.name, 'volume')
        self.assertEquals(vol.comment(), "Volume 'volume': cell 'fuel' 1 cm^3.")
        self.assertEquals(vol.mcnp('%g', self.sim), "VOL 1 5J")
        vol = cards.Volume('D', 1)
        self.assertEquals(vol.mcnp('%g', self.sim), "VOL 4J 1 1J")
        # Multiple cells
        vol = cards.Volume('fuel', 1, 'coolant', 2, manual=True)
        self.assertEquals(vol.comment(), "Volume 'volume': (all manual) "
                "cell 'fuel' 1 cm^3, cell 'coolant' 2 cm^3.")
        self.assertEquals(vol.mcnp('%g', self.sim), "VOL NO 1 2 4J")
        # Using set()
        vol.set('coolant', 3)
        self.assertEquals(vol.mcnp('%g', self.sim), "VOL NO 1 3 4J")

        ## Area
        are = cards.Area('fuel', 10)
        self.assertEquals(are.name, 'area')
        self.assertEquals(are.comment(), "Area 'area': cell 'fuel' 10 cm^2.")
        self.assertEquals(are.mcnp('%g', self.sim), "AREA 10 5J")
        are = cards.Area('D', 10)
        self.assertEquals(are.mcnp('%g', self.sim), "AREA 4J 10 1J")
        # Multiple cells
        are = cards.Area('fuel', 10, 'coolant', 20)
        self.assertEquals(are.comment(), "Area 'area': "
                "cell 'fuel' 10 cm^2, cell 'coolant' 20 cm^2.")
        self.assertEquals(are.mcnp('%g', self.sim), "AREA 10 20 4J")
        # Using set()
        are.set('coolant', 30)
        self.assertEquals(are.mcnp('%g', self.sim), "AREA 10 30 4J")

        ## TemperatureTimes
        thtme = cards.TemperatureTimes([1e10, 2e10])
        self.assertEquals(thtme.name, 'temptimes')
        self.assertEquals(thtme.comment(), "Temperature times 'temptimes' "
                "(in MeV): 1e+10 2e+10.")
        self.assertEquals(thtme.mcnp('%g', self.sim), "THTME 1e+18 2e+18")

        ## Temperature
        temp = cards.Temperature('fuel', 600)
        self.assertEquals(temp.name, 'temperature')
        self.assertEquals(temp.comment(), "Temperatures 'temperature': "
                "cell 'fuel' 600 K.")
        self.assertEquals(temp.mcnp('%g', self.sim),
                "TMP 5.17041e-08 5J")
        temp = cards.Temperature('fuel', 600, 'E', 900, index=2)
        self.assertEquals(temp.comment(), "Temperatures for time index "
                "2 'temperature-idx2': cell 'fuel' 600 K, "
                "cell 'E' 900 K.")
        self.assertEquals(temp.mcnp('%g', self.sim),
                "TMP2 5.17041e-08 4J 7.75561e-08")
        # set()
        temp = cards.Temperature(index=2)
        temp.set('fuel', 600)
        temp.set('E', 900)
        self.assertEquals(temp.comment(), "Temperatures for time index "
                "2 'temperature-idx2': cell 'fuel' 600 K, "
                "cell 'E' 900 K.")
        self.assertEquals(temp.mcnp('%g', self.sim),
                "TMP2 5.17041e-08 4J 7.75561e-08")
        # Modifying.
        temp.set('E', 950)
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
        dxc = cards.DXTRANContribution('neutron', None, 'fuel', 0.5)
        self.assertEquals(dxc.name, 'dxtrancont-neutron')
        self.assertEquals(dxc.comment(), "DXTRAN contribution "
                "for all spheres 'dxtrancont-neutron': cell 'fuel' 0.5.")
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC:N 0.5 5J")
        dxc = cards.DXTRANContribution('neutron', 'sph2', 'fuel', 0.5)
        self.assertEquals(dxc.comment(), "DXTRAN contribution "
                "for sphere 'sph2' 'dxtrancont-sph2-neutron': "
                "cell 'fuel' 0.5.")
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC2:N 0.5 5J")
        dxc = cards.DXTRANContribution('neutron', 'sph2', 'fuel', 0.5,
                                                       'coolant', 0.7)
        self.assertEquals(dxc.comment(), "DXTRAN contribution "
                "for sphere 'sph2' 'dxtrancont-sph2-neutron': "
                "cell 'fuel' 0.5, cell 'coolant' 0.7.")
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC2:N 0.5 0.7 4J")
        dxc = cards.DXTRANContribution('neutron', 'sph2')
        dxc.set('fuel', 0.5)
        dxc.set('coolant', 0.7)
        self.assertEquals(dxc.comment(), "DXTRAN contribution "
                "for sphere 'sph2' 'dxtrancont-sph2-neutron': "
                "cell 'fuel' 0.5, cell 'coolant' 0.7.")
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC2:N 0.5 0.7 4J")
        dxc.set('coolant', 0.8)
        self.assertEquals(dxc.mcnp('%g', self.sim), "DXC2:N 0.5 0.8 4J")

        ## FissionTurnoff
        # Default.
        fto = cards.FissionTurnoff()
        self.assertEquals(fto.name, 'fissionturnoff')
        self.assertEquals(fto.comment(), "Fission turnoff 'fissionturnoff': "
                "capture-gamma for all cells.")
        self.assertEquals(fto.mcnp('%g', self.sim), "NONU")
        # With an argument.
        fto = cards.FissionTurnoff('fuel', 'real-gamma',
                                   'coolant', 'capture-nogamma')
        self.assertEquals(fto.comment(), "Fission turnoff 'fissionturnoff': "
                "cell 'fuel' real-gamma, cell 'coolant' capture-nogamma.")
        self.assertEquals(fto.mcnp('%g', self.sim), "NONU 1 2 4J")
        # set()
        fto = cards.FissionTurnoff()
        fto.set('fuel', 'real-gamma')
        fto.set('coolant', 'capture-nogamma')
        self.assertEquals(fto.comment(), "Fission turnoff 'fissionturnoff': "
                "cell 'fuel' real-gamma, cell 'coolant' capture-nogamma.")
        self.assertEquals(fto.mcnp('%g', self.sim), "NONU 1 2 4J")
        # Modifying.
        fto.set('coolant', 'capture-gamma')
        self.assertEquals(fto.mcnp('%g', self.sim), "NONU 1 0 4J")

        ## DetectorContribution
        det = cards.PointDetector('point', 'neutron', [0, 0, 0], 0, 'cm')
        det2 = cards.PointDetector('point2', 'neutron', [0, 0, 0], 0, 'cm')
        self.sim.add_tally(det)
        self.sim.add_tally(det2)
        dc = cards.DetectorContribution('point2', 'fuel', 0.5)
        self.assertEquals(dc.name, 'detcontrib-point2')
        self.assertEquals(dc.comment(), "Detector contribution "
                "'detcontrib-point2': cell 'fuel' 0.5.")
        self.assertEquals(dc.mcnp('%g', self.sim), "PD25 0.5 5J")
        dc = cards.DetectorContribution('point2', 'fuel', 0.5, 'coolant',
                0.6)
        self.assertEquals(dc.comment(), "Detector contribution "
                "'detcontrib-point2': cell 'fuel' 0.5, cell 'coolant' 0.6.")
        self.assertEquals(dc.mcnp('%g', self.sim), "PD25 0.5 0.6 4J")
        dc = cards.DetectorContribution('point2')
        dc.set('fuel', 0.5)
        dc.set('coolant', 0.6)
        self.assertEquals(dc.comment(), "Detector contribution "
                "'detcontrib-point2': cell 'fuel' 0.5, cell 'coolant' 0.6.")
        self.assertEquals(dc.mcnp('%g', self.sim), "PD25 0.5 0.6 4J")
        dc.set('coolant', 0.7)
        self.assertEquals(dc.mcnp('%g', self.sim), "PD25 0.5 0.7 4J")

        ## PhotonWeight
        pw = cards.PhotonWeight()
        self.assertEquals(pw.name, 'photonweight')
        pw.set('fuel', 'off')
        self.assertEquals(pw.comment(), "Photon weight thresholds "
                "'photonweight': cell 'fuel' off.")
        self.assertEquals(pw.mcnp('%g', self.sim), "PWT -1.0E6 5J")
        pw.set('coolant', 'one')
        self.assertEquals(pw.comment(), "Photon weight thresholds "
                "'photonweight': cell 'fuel' off, "
                "cell 'coolant' one.")
        self.assertEquals(pw.mcnp('%g', self.sim), "PWT -1.0E6 0 4J")
        pw.set('coolant', 0.5)
        self.assertEquals(pw.comment(), "Photon weight thresholds "
                "'photonweight': cell 'fuel' off, "
                "cell 'coolant' 0.5.")
        self.assertEquals(pw.mcnp('%g', self.sim), "PWT -1.0E6 0.5 4J")
        pw.set('coolant', 0.5, pre_weight=True)
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
        tally = cards.SurfaceCurrent('fuel', 'electron', ['fuelpin',
                'bound'], total=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle.name, 'electron')
        self.assertIs(tally.cards[0], 'fuelpin')
        self.assertIs(tally.cards[1], 'bound')
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

        tally = cards.SurfaceCurrent('fuel2', 'photon', [['fuelpin',
                'bound']], alt_units=True)
        self.assertEquals(tally.particle.name, 'photon')
        self.assertTrue(len(tally.cards), 1)
        self.assertTrue(len(tally.cards[0]), 2)
        self.assertIs(tally.cards[0][0], 'fuelpin')
        self.assertIs(tally.cards[0][1], 'bound')
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
        tally = cards.SurfaceFlux('fuel', 'electron', ['fuelpin',
                'bound'], average=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle.name, 'electron')
        self.assertIs(tally.cards[0], 'fuelpin')
        self.assertIs(tally.cards[1], 'bound')
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

        tally = cards.SurfaceFlux('fuel2', 'proton', [['fuelpin', 'bound']],
                    alt_units=True)
        self.assertEquals(tally.particle.name, 'proton')
        self.assertTrue(len(tally.cards), 1)
        self.assertTrue(len(tally.cards[0]), 2)
        self.assertIs(tally.cards[0][0], 'fuelpin')
        self.assertIs(tally.cards[0][1], 'bound')
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
        tally = cards.CellFlux('fuel', 'neutron', 'fuel')
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle.name, 'neutron')
        self.assertEquals(tally.comment(), "Cell flux tally 'fuel' "
                "of neutrons: cell 'fuel'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F14:N  1")
        #self.assertIs(tally._unique_card_list()[0], 'fuel')
        # Two individual cells.
        tally = cards.CellFlux('both', 'neutron', ['fuel', 'coolant'])
        self.assertEquals(tally.comment(), "Cell flux tally 'both' "
                "of neutrons: cells 'fuel'; 'coolant'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F24:N  1 2")
        #self.assertIs(tally._unique_card_list()[0], 'fuel')
        #self.assertIs(tally._unique_card_list()[1], 'coolant')
        # Two individual cells, with average over all.
        tally = cards.CellFlux('withavg', 'neutron', ['fuel', 'coolant'],
                average=True)
        self.assertEquals(tally.comment(), "Cell flux tally 'withavg' "
                "of neutrons: cells 'fuel'; 'coolant'; and avg. of all "
                "provided.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F34:N  1 2 T")
        #self.assertIs(tally._unique_card_list()[0], 'fuel')
        #self.assertIs(tally._unique_card_list()[1], 'coolant')
        # Two individual cells, and an averaging, with an average over all.
        tally = cards.CellFlux('withavg', 'neutron', ['fuel',
                ['fuel', 'coolant'], 'coolant'], average=True)
        self.assertEquals(tally.comment(), "Cell flux tally 'withavg' "
                "of neutrons: cells 'fuel'; avg. in 'fuel', 'coolant'; "
                "'coolant'; and avg. of all provided.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F34:N  1 ( 1 2) 2 T")
        #self.assertIs(tally._unique_card_list()[0], 'fuel')
        #self.assertIs(tally._unique_card_list()[1], 'coolant')
        #self.assertTrue(len(tally._unique_card_list()) == 2)

        ## CellEnergyDeposition
        tally = cards.CellEnergyDeposition('energy', 'neutron', 'fuel')
        self.assertEquals(tally.name, 'energy')
        self.assertEquals(tally.particle.name, 'neutron')
        self.assertIs(tally.cards, 'fuel')
        # comment()
        self.assertEquals(tally.comment(), "Energy deposition tally "
                "'energy' of neutrons: cell 'fuel'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F16:N  1")
        tally = cards.CellEnergyDeposition('energy2', ['neutron', 'proton'],
                'fuel')
        self.assertIs(type(tally.particle), list)
        self.assertEquals(tally.particle[0].name, 'neutron')
        self.assertEquals(tally.particle[1].name, 'proton')
        self.assertEquals(tally.comment(), "Energy deposition tally " 
                "'energy2' of neutrons, protons: cell 'fuel'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F26:N,H  1")
        tally = cards.CellEnergyDeposition('energy3', 'all', 'fuel')
        self.assertEquals(tally.comment(), "Energy deposition tally "
                "'energy3' of all: cell 'fuel'.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "+F36  1")
        # Test exceptions.
        self.assertRaises(ValueError, cards.CellEnergyDeposition, 'energy',
                ['neutron', 'all'], 'fuel')
        self.assertRaises(ValueError, cards.CellEnergyDeposition,
                'energy', 'all', 'fuel', alt_units=True)

        ## CellFissionEnergyDeposition
        tally = cards.CellFissionEnergyDeposition('fuel', ['fuel',
                'coolant'], average=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle.name, 'neutron')
        self.assertIs(type(tally.cards), list)
        self.assertIs(tally.cards[0], 'fuel')
        self.assertIs(tally.cards[1], 'coolant')
        self.assertTrue(tally.average)
        self.assertFalse(tally.alt_units)
        self.assertEquals(tally.comment(), "Fission energy deposition tally "
                "'fuel' of neutrons: cells 'fuel'; 'coolant'; and avg. of "
                "all provided.")
        # mcnp()
        self.sim.add_tally(tally)
        self.assertEquals(tally.mcnp('%5g', self.sim), 
                "F17:N  1 2 T")
        tally = cards.CellFissionEnergyDeposition('fuel', [['fuel',
                'coolant']], alt_units=True)
        self.assertTrue(len(tally.cards), 1)
        self.assertTrue(len(tally.cards[0]), 2)
        self.assertIs(tally.cards[0][0], 'fuel')
        self.assertIs(tally.cards[0][1], 'coolant')
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
        tally = cards.CellPulseHeight('fuel', ['proton', 'electron'], ['fuel',
                'coolant'], alt_units=True)
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle[0].name, 'proton')
        self.assertEquals(tally.particle[1].name, 'electron')
        self.assertIs(tally.cards[0], 'fuel')
        self.assertIs(tally.cards[1], 'coolant')
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
                ['fuel', 'coolant'])
        self.assertEquals(tally.name, 'fuel')
        self.assertEquals(tally.particle[0].name, 'proton')
        self.assertEquals(tally.particle[1].name, 'electron')
        self.assertIs(tally.cards[0], 'fuel')
        self.assertIs(tally.cards[1], 'coolant')
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
        egrid = cards.EnergyGrid('grid1', 'ring', np.array([1, 2]))
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

    def test_CustomCards(self):
        """Tests :py:class:`pyne.simplesim.cards.Custom and subclasses."""

        fname = 'simplesim_custom'

        # Define system.
        # Materials.
        uo2 = material.from_atom_frac({'U235': 0.05, 'U238': 0.95, 'O16' : 2.00})
        uo2 = cards.Material(uo2, name='UO2')
        h2o = material.from_atom_frac({'H1' : 2.0, 'O16': 1.0}, attrs={'name': 'H2O'})
        h2o = cards.Material(h2o)
        mc = cards.MaterialCustom(mat=material.Material(), name='matc', 
                                  comment="mathey", mcnp="mathey")
        self.sim.sys.add_material(mc)

        # Surfaces.
        sc = cards.SurfaceCustom('surfc', comment="surfhey", mcnp="surfhey")
        self.sim.sys.add_surface(sc)

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
        coolant = cards.CellMCNP('coolant', pin.pos | cellbound.neg, h2o,
                1.0, 'g/cm^3',
                importance=('neutron', 1),
                volume=1)
        graveyard = cards.CellMCNP('graveyard', cellbound.pos,
                importance=('neutron', 0))
        cc = cards.CellCustom('cellc', comment="cellhey", mcnp="cellhey")

        # Add cells to the system.
        self.sys.add_cell(fuel)
        self.sys.add_cell(cc)
        self.sys.add_cell(coolant)
        self.sys.add_cell(graveyard)

        # Add source, tallies to simulation.
        self.sim.add_misc(cards.ScatteringLaw('H2O', {'H1': 'lwtr'}))
        self.sim.add_source(cards.Criticality())
        self.sim.add_source(cards.CriticalityPoints())
        self.sim.add_tally(cards.CellFlux('flux', 'neutron', 
                ['fuel', 'coolant']))
        self.sim.add_misc(cards.EnergyGrid('egrid0', None,
                10**np.arange(-9.9, 1.1, 0.1)))

        # Customs.
        sc = cards.SourceCustom('sourcec', comment="sourcehey",
            mcnp="sourcehey")
        self.sim.add_source(sc)

        dc = cards.DistributionCustom('distc', comment="disthey",
            mcnp=[('SI', "disthey"), ('SP', 'distheyp')])
        self.sim.add_dist(dc)

        # More complicated.
        t1 = cards.TallyCustom('tallyc', comment="tallyhey", mcnp="tallyhey")
        t2 = cards.TallyCustom('tally2', comment="tallyhey", mcnp="tallyhey",
                tallyclass=cards.CellFlux)
        t3 = cards.TallyCustom('tally3', comment="tallyhey", mcnp="tallyhey",
                tallyclass=cards.CellFlux,
                mcnp_pre='*')
        self.sim.add_tally(t1)
        self.sim.add_tally(t2)
        self.sim.add_tally(t3)

        mc = cards.MiscCustom('miscc', comment="mischey", mcnp="mischey")
        self.sim.add_misc(mc) 
        
        tc = cards.TransformationCustom('transc', comment="transhey",
                mcnp="transhey")
        self.sim.add_transformation(tc)

        # Create input file.
        inp = inputfile.MCNPInput(self.sim, title="Infinite lattice.")
        inp.write(fname)
        # Check against text file.
        self.assertEquals(
                open(fname).readlines(),
                open(fname + '_compare').readlines())
        os.unlink(fname)


    def test_InfLattice(self):
        """Tests the input file for an infinite lattice reactor. Checks
        generated output against the text file `inflattice_compare`.

        """
        fname = 'simplesim_inflattice'

        # Define system.
        # Materials.
        uo2 = material.from_atom_frac({'U235': 0.05, 'U238': 0.95, 'O16' : 2.00})
        uo2 = cards.Material(uo2, name='UO2')
        h2o = material.from_atom_frac({'H1' : 2.0, 'O16': 1.0}, attrs={'name': 'H2O'})
        h2o = cards.Material(h2o)

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
        coolant = cards.CellMCNP('coolant', pin.pos | cellbound.neg, h2o,
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
        self.sim.add_tally(cards.CellFlux('flux', 'neutron', 
                ['fuel', 'coolant']))
        self.sim.add_misc(cards.EnergyGrid('egrid0', None,
                10**np.arange(-9.9, 1.1, 0.1)))

        # Create input file.
        inp = inputfile.MCNPInput(self.sim, title="Infinite lattice.")
        inp.write(fname + '_default')
        # Check against text file.
        self.assertEquals(
                open(fname + '_default').readlines(),
                open(fname + '_default_compare').readlines())
        os.unlink(fname + '_default')

        # Test the & line continuation.
        inp = inputfile.MCNPInput(self.sim, title="Infinite lattice.",
                cont_by_amp=True)
        inp.write(fname + '_amp')
        self.assertEquals(
                open(fname + '_amp').readlines(),
                open(fname + '_amp_compare').readlines())
        os.unlink(fname + '_amp')

    def test_user_card(self):
        """Tests the functionality that a user can add their own card."""

        fname = 'simplesim_user'
        inp = inputfile.MCNPInput(self.sim, float_format="%.5e", title="1")
        inp.add_user_card('cell', '11 0 -5 IMP:N=0', comment='Graveyard.')
        inp.add_user_literal('cell', 'M1 1001 1\n     8016 2')
        inp.add_user_card('surface', '11 0 -5 IMP:N=0', comment='Graveyard.')
        inp.add_user_literal('surface', 'M1 1001 1\n     8016 2')
        inp.add_user_card('data', '11 0 -5 IMP:N=0', comment='Graveyard.')
        inp.add_user_literal('data', 'M1 1001 1\n     8016 2')
        inp.write(fname + '_first')
        self.assertEquals(
                open(fname + '_first').readlines(),
                open(fname + '_first_compare').readlines())
        os.unlink(fname + '_first')
        

class TestModifying(unittest.TestCase):
    """Ensures that the code is amenable to modifications made after initial
    creation.
    
    """
    pass


## The following tests are redundant, but are to make sure that the examples in
# the documentation function as expected.
class TestGodiva(unittest.TestCase):
    pass
