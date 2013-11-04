.. _simplesim_ex_inflattice:

========================
Infinite lattice example
========================
    
We create an MCNP input file named `ex_simplesim_inflattice` for a 2-cell (fuel
and moderator) infinite lattice reactor.

This example can also be found in `pyne/examples/simplesim.py`.

************
Starting out
************

It'll be helpful to wrap our code in our own class::


    class InfLattice(object):
        def __init__(self):

We define the system::

        ## Define the system: materials, surfaces, regions, cells.
        self.sys = definition.SystemDefinition(verbose=False)

First come materials::

        ## Materials.
        # Must provide a name as a keyword argument for material cards. See the
        # documentation for :py:mod:`pyne.material` for more information.
        uo2 = material.from_atom_frac({'U235': 0.05, 'U238': 0.95, 'O16' : 2.00})
        self.uo2 = cards.Material(uo2, name='UO2')

        h2o = material.from_atom_frac({'H1' : 2.0, 'O16': 1.0}, attrs={'name': 'H2O'})
        self.h2o = cards.Material(h2o)

Then surfaces and regions::

        ## Surfaces.
        # There are two surfaces: one for the pin and one for the unit cell
        # boundary.
        radius = 0.40

        # This creates an axis-aligned and axis-centered cylinder along the z
        # axis, with radius 0.40 cm.
        self.pin = cards.AxisCylinder('pin', 'Z', radius)

        # The Parallelepiped is a macrobody. The surface is reflecting,
        # creating an infinte geometry. The surface is infinite in the z
        # direction.
        pitch = 1.2
        self.cellbound = cards.Parallelepiped('bound',
                -pitch / 2, pitch / 2, -pitch / 2, pitch / 2, 0, 0,
                reflecting=True)

Then we create cells::

        ## Cells.
        # We combine the materials and surfaces above into cells. We use MCNP
        # cells in order to specify particle importances and volumes directly
        # on the cell card. We could alternatively use the
        # :py:class:`Importance` and :py:class:`Volume` cards.

        # fuel cell.
        # The fuel is the region of space inside the pin, pin.neg. 
        self.fuelregion = self.pin.neg

        # The neutron importance is 1, and the user-provided volume is 1 cm^3.
        self.fuel = cards.CellMCNP('fuel', self.fuelregion, self.uo2,
                11.0, 'g/cm^3',
                importance=('neutron', 1),
                volume=1)

        # coolant cell.
        # The region is between the pin and the unit cell boundary.
        self.coolantregion = self.pin.pos | self.cellbound.neg
        self.coolant = cards.CellMCNP('coolant', self.coolantregion, self.h2o,
                1.0, 'g/cm^3',
                importance=('neutron', 1),
                volume=1)

        # graveyard cell: where particles go to die.
        # The region is everything beyond the unit cell boundary.
        self.graveyardregion = self.cellbound.pos

        # This is a void cell, meaning it does not have a material.
        self.graveyard = cards.CellMCNP('graveyard', self.graveyardregion,
                importance=('neutron', 0))

        # We add the cells to the system. The order we add them is the order
        # they are printed in the input file.
        self.sys.add_cell(self.fuel)

        # We can add multiple cells at once.
        self.sys.add_cell(self.coolant, self.graveyard)
       
Now we define the simulation of the system::

        ## Define the simulation: sources, tallies, misc. Don't clutter the
        # command window.
        self.sim = definition.MCNPSimulation(self.sys, verbose=False)

We have source and tally cards::

        # Specify a thermal scattering law for the H2O material. This is a
        # unique card per material.
        self.sim.add_misc(cards.ScatteringLaw('H2O', {'H1': 'lwtr'}))

        # Add a criticality source, use default values. This is a unique card,
        # so we do not provide a card name.
        self.sim.add_source(cards.Criticality())

        # Add points at which to start neutrons; use default point (0, 0, 0).
        self.sim.add_source(cards.CriticalityPoints())

        # Tally neutron flux in both the fuel and coolant cells.
        self.sim.add_tally(cards.CellFlux('flux', 'neutron', 
                                          ['fuel', 'coolant']))

        # The energy grid on which to tally neutrons, applied to all tallies.
        self.sim.add_misc(cards.EnergyGrid('egrid0', None,
                                           10**np.arange(-9.9, 1.1, 0.1)))

That's it for the constructor. In our class we define the following method that
actually creates the input::

    def write(self):
        """Writes the input to 'ex_simplesim_inflattice'."""

        # Create input file, specifying the title of the input.
        self.inp = inputfile.MCNPInput(self.sim, title="Infinite lattice.")
        self.inp.write('ex_simplesim_inflattice')

That's it for the class. We then use the class::

    # Create all relevant objects for the infinite lattice example.
    inflat = InfLattice()
    # Write to a file.
    inflat.write()

See below for what this generates.


************************
Playing around afterward
************************
We might want to do a second run of the code for different parameters. 
It's super easy to update our input file. First let's change the radius 
of the pin and write the input again::

        inflat.pin.radius = 0.45
        inflat.write()

If you open the input you'll find the radius has been updated. Sometimes 
we don't have the card object, and so to change a card we need to access 
it from its place in a dictionary in the definition. In this case we need 
to know its category, and its name. Here is an example of how we access a 
source card whose name we know, and how we can modify it::

        inflat.sim.source['criticality'].keff_guess = 1.5
        inflat.write()

The input is updated as we hoped.


*****************************
The output (ie an input deck)
*****************************

This is what is generated::

    Infinite lattice.
      C Generated with the Python package PyNE (pyne.github.com).
      C ==========
      C Cell Cards
      C ==========
      C Cell 'fuel': region -pin, material 'UO2' density 11 g/cm^3 VOL= 1 cm^3
      c     IMP:N= 1.
    1 1 -11 -1 VOL= 1 IMP:N=1
      C
      C Cell 'coolant': region (+pin | -bound), material 'H2O' density 1 g/cm^3 VOL=
      c     1 cm^3 IMP:N= 1.
    2 2 -1 (1 -2) VOL= 1 IMP:N=1
      C
      C Cell 'graveyard': region +bound, void IMP:N= 0.
    3 0 2 IMP:N=0
      C
    
      C =============
      C Surface Cards
      C =============
      C Axis cylinder 'pin': aligned and centered on z axis, with radius 0.4 cm
      c     (diameter 0.8 cm).
    1  CZ   0.4
      C
      C Parallelepiped 'bound': reflecting. [-0.6, 0.6] x [-0.6, 0.6] x [0, 0] cm.
    *2 RPP -0.6  0.6  -0.6  0.6   0  0
      C
    
      C ==========
      C Data Cards
      C ==========
      C
      C **************
      C Material Cards
      C **************
      C Material 'UO2'.
    M1
           8016    2 $ O16
          92235    0.05 $ U235
          92238    0.95 $ U238
      C
      C Material 'H2O'.
    M2
           1001    2 $ H1
           8016    1 $ O16
      C
      C
      C ************
      C Source Cards
      C ************
      C Criticality source 'criticality': n_histories: 1000, keff_guess: 1,
      c     n_skip_cycles: 30, n_cycles: 130.
    KCODE 1000  1 30 130
      C
      C Criticality points 'criticalitypoints': (0, 0, 0) cm.
    KSRC  0  0  0
      C
      C
      C ***********
      C Tally Cards
      C ***********
      C Cell flux tally 'flux' of neutrons: cells 'fuel'; 'coolant'.
    F14:N  1 2
      C
      C
      C *******************
      C Miscellaneous Cards
      C *******************
      C Scattering law 'scatlaw-H2O': H1: lwtr.
    MT2 lwtr
      C
      C Energy grid 'egrid0' for all tallies: 110 groups.
    E0  1.2589e-10  1.5849e-10  1.9953e-10  2.5119e-10  3.1623e-10  3.9811e-10
         5.0119e-10  6.3096e-10  7.9433e-10  1e-09  1.2589e-09  1.5849e-09
         1.9953e-09  2.5119e-09  3.1623e-09  3.9811e-09  5.0119e-09  6.3096e-09
         7.9433e-09  1e-08  1.2589e-08  1.5849e-08  1.9953e-08  2.5119e-08
         3.1623e-08  3.9811e-08  5.0119e-08  6.3096e-08  7.9433e-08  1e-07
         1.2589e-07  1.5849e-07  1.9953e-07  2.5119e-07  3.1623e-07  3.9811e-07
         5.0119e-07  6.3096e-07  7.9433e-07  1e-06  1.2589e-06  1.5849e-06
         1.9953e-06  2.5119e-06  3.1623e-06  3.9811e-06  5.0119e-06  6.3096e-06
         7.9433e-06  1e-05  1.2589e-05  1.5849e-05  1.9953e-05  2.5119e-05
         3.1623e-05  3.9811e-05  5.0119e-05  6.3096e-05  7.9433e-05  0.0001
         0.00012589  0.00015849  0.00019953  0.00025119  0.00031623  0.00039811
         0.00050119  0.00063096  0.00079433  0.001  0.0012589  0.0015849  0.0019953
         0.0025119  0.0031623  0.0039811  0.0050119  0.0063096  0.0079433  0.01
         0.012589  0.015849  0.019953  0.025119  0.031623  0.039811  0.050119
         0.063096  0.079433  0.1  0.12589  0.15849  0.19953  0.25119  0.31623
         0.39811  0.50119  0.63096  0.79433  1  1.2589  1.5849  1.9953  2.5119
         3.1623  3.9811  5.0119  6.3096  7.9433  10
      C


If we don't want all the comments, we can use a keyword argument on the input
file initialization::

    self.inp = inputfile.MCNPInput(self.sim, title="Infinite lattice.",
            comments=False)


