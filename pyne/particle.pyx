"""Python wrapper for particle library.

The particle module provides a string & integer based canonical form for particle
idenfification, and a more human readable form"

  **Particle Names:** The names are chosen to be concatenated camel case form of
  regular particle names. For example, "Neutron" or "AntiProton".

  Some common examples & their alternative names
  
  ========  =======
  Particle  Altname  
  ========  =======
  Neutron      - 
  Proton    Hydrogen, Protium, H1 etc
  Photon    Gamma, X-Ray
  Electron  Beta-, Beta

  All alternaive particle names map back to the cannonical form, e.g. Hydrogen 
  maps  back to Proton

  A valid particle name is ANY particle from the names map, i.e. those that MAP 
 to the Berkley Particle Data Centre numbering scheme [PDC]_, and any valid 
 nucid/name, like H1 or 10010000, this is to support the high energy physics 
 codes that some users are starting to employ.

  **PDC Numbers:** The PDC numbers map from those specified by the Berkley 
  Particle Data Centre and only have valid numbers for those in the name map, 
  some valid particle names, like H1, have a PDC number of 0.

  **Doc strings:** The doc strings are human readable versions of the name, 
  with grammar and spaces, like what humans enjoy.
  
The particle module has a variety of functions to convert between forms, for 
example the name() function will take either the full name or the pdc number,
or a valid nucid and return in its true name form. The same is true of the 
id() function, returning the valid PDC number of the particle.

------------------------
Adding more particles!
------------------------

Adding more particles will be nessesary in time, but it is very straight
forward to add more.

1. add to the define NUM_PARTICLES variable, the number of particles you 
 are going to add
2. add the name to the ``_names`` array, add the pdc number to the 
 ``_pdcids`` array.
3. in the _fill_maps() function add the human readable form to the 
 ``_docs`` array.
4. in the _fill_maps() function add to the altnames map if there are other forms.
5. in the _fill_maps() function add the mcnp and fluka forms to part2* arrays.

------------------------

.. [PDC] http://pdg.lbl.gov/2014/mcdata/mass_width_2014.mcd

------------------------

"""
from __future__ import unicode_literals

# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from libc.string cimport const_char
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libcpp.string cimport string as std_string

from pyne.utils import QA_warn

# local imports 
cimport extra_types
cimport pyne.cpp_utils
cimport pyne.pyne_config
import pyne.pyne_config
cimport cpp_nucname
cimport cpp_particle
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv


QA_warn(__name__)

# names
cdef conv._SetStr names_proxy = conv.SetStr(False)
names_proxy.set_ptr = &cpp_particle.names
names = names_proxy

"""
# docs
cdef conv._MapStrStr docs_proxy = conv.SetStrStr(False)
docs_proxy.set_ptr = &cpp_particle.docs
docs = docs_proxy
"""

# altnames
cdef conv._MapStrInt altnames_proxy = conv.MapStrInt(False)
altnames_proxy.map_ptr = &cpp_particle.altnames
altnames = altnames_proxy

# id_name
cdef conv._MapIntStr id_name_proxy = conv.MapIntStr(False)
id_name_proxy.map_ptr = &cpp_particle.id_name
id_name = id_name_proxy

# name_id
cdef conv._MapStrInt name_id_proxy = conv.MapStrInt(False)
name_id_proxy.map_ptr = &cpp_particle.name_id
name_id = name_id_proxy

def name(x):
    """name(x)
    
    Gives the unique particle name.

    Parameters
    ----------
    x : str, int, or char*

    Returns
    -------
    n : str Unique particle name
    """
    if isinstance(x, basestring):
        x_bytes = x.encode()
        cn = cpp_particle.name(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.name(<int> long(x))

    n = bytes(<char *> cn.c_str()).decode()
    
    return n

def mcnp(x):
    """mcnp(x)
    
    Gives the unique mcnp name for particle.

    Parameters
    ----------
    x : str, int, or char*

    Returns
    -------
    n : str Unique particle name
    """
    if isinstance(x, basestring):
        x_bytes = x.encode()
        cn = cpp_particle.mcnp(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.mcnp(<int> long(x))

    n = bytes(<char *> cn.c_str()).decode()
    
    return n


def mcnp6(x):
    """mcnp6(x)
    
    Gives the unique mcnp6 name for particle.

    Parameters
    ----------
    x : str, int, or char*

    Returns
    -------
    n : str Unique particle name
    """
    if isinstance(x, basestring):
        x_bytes = x.encode()
        cn = cpp_particle.mcnp6(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.mcnp6(<int> long(x))

    n = bytes(<char *> cn.c_str()).decode()
    
    return n


def fluka(x):
    """fluka(x)
    
    Gives the unique fluka name for particle.

    Parameters
    ----------
    x : str, int, or char*

    Returns
    -------
    n : str Unique particle name
    """
    if isinstance(x, basestring):
        x_bytes = x.encode()
        cn = cpp_particle.fluka(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.fluka(<int> long(x))

    n = bytes(<char *> cn.c_str()).decode()
    
    return n

def geant4(x):
    """geant4(x)
    
    Gives the unique geant name for particle.

    Parameters
    ----------
    x : str, int, or char*

    Returns
    -------
    n : str Unique particle name
    """
    if isinstance(x, basestring):
        x_bytes = x.encode()
        cn = cpp_particle.geant4(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.geant4(<int> long(x))

    n = bytes(<char *> cn.c_str()).decode()
    
    return n


def describe(x):
    """describe(x)
    
    Gives the unique particle name description.

    Parameters
    ----------
    x : str, int, or char*

    Returns
    -------
    n : str Unique particle description
    """

    if isinstance(x, basestring):
        x_bytes = x.encode()
        cn = cpp_particle.describe(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.describe(<int> long(x))

    n = bytes(<char *> cn.c_str()).decode()
    
    return n

def id(x):
    """id(x)
    
    Gives the unique particle id number

    Parameters
    ----------
    x : str, int, or char*

    Returns
    -------
    n : int Unique PDC number
    """

    if isinstance(x, basestring):
        x_bytes = x.encode()
        cn = cpp_particle.id(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.id(<int> long(x))

    n = cn
    
    return n

def is_valid(x):
    """is_valid(x)
    
    Is the argument a valid particle id

    Parameters
    ----------
    x : str, int, or char*

    Returns
    -------
    n : bool true/false 
    """

    if isinstance(x, basestring):
        x_bytes = x.encode()
        cn = cpp_particle.is_valid(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.is_valid(<int> long(x))

    n = cn
    
    return n

def is_heavy_ion(x):
    """is_heavy_ion(x)
    
    Is the argument a valid particle id

    Parameters
    ----------
    x : str, int, or char*

    Returns
    -------
    n : bool true/false 
    """

    if isinstance(x, basestring):
        x_bytes = x.encode()
        cn = cpp_particle.is_heavy_ion(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.is_heavy_ion(<int> long(x))

    n = cn
    
    return n


