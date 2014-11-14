.. currentmodule:: pyne.particle

.. _usersguide_particle:

==========================
Particle Naming Conventions
==========================
Transport physics codes like MCNP, Fluka, Geant4 etc each identify 
particles, in their own unique way, MCNP for example uses single letter 
mneumonics like 'n' for neutron, Fluka uses the full name of the particle 
or its own particle id number like 'NEUTRON'. The particle class aims to 
have a single unified naming scheme, a Camel Case like name for humans, 
like "MuonNeutrino", the fundamental particles map to Berkley Particle
Data Center numbering schemes, and some physics codes allow the scoring 
of so called 'Heavy Ions', which is handled using the pyne nucname module. 
So the particle identification methods use:

1. **name**: The human readable Camel Case form of the name.
2. **id**: The PDC id number.

Thus by using either the particle name or the PDC number, you can be 
guarenteed to map to the same particle. There are also helper functions 
like, is_valid to determine if your partice id is a valid one.

Valid particle names can be output using the code specific output 
functions, if unknown so an invalid particle name is generated:

1. **mcnp**: The MCNP4/5 equivalent of the particle name, if unknown returns '?'
2. **mcnp6**: The MCNP6 equivalent of the particle name, if unknown returns '?'
3. **fluka**: The Fluka equiavlent of the particle name, if unknown returns 
 '????????'
4. **geant4**: The Geant4 equivalent of the particle name, if unknown returns  
 '????????'


---------------
Examples of Use
---------------

.. code-block:: ipython

    In [1]: from pyne import particle

    In [2]: particle.name('Proton')
    Out[2]: Proton

    In [3]: particle.name('Hydrogen')
    Out[3]: Proton

    In [4]: particle.name('Protium')
    Out[4]: Proton

    In [5]: particle.name('Neutron')
    Out[5]: Neutron

    In [6]: particle.is_valid('AM242M')
    Out[6]: True

    In [7]: particle.is_valid['Clearly not a particle name']
    Out[7]: False

    In [8]: particle.is_heavy_ion['AM242M]
    Out[8]: True

    In [9]: particle.is_heavy_ion('Hydrogen')
    Out[9]: False

    In [10]: x = particle.name('Hydrogen')
    In [11]: print x.mcnp()
    Out[11]: ?


Further information on the particle module may be seen in the library reference 
:ref:`pyne_particle`.
