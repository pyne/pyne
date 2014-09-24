"""Python wrapper for particle library.

"""
from __future__ import unicode_literals

# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from libc.string cimport const_char
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libcpp.string cimport string as std_string

from warnings import warn
from pyne.utils import VnVWarning

# local imports 
cimport extra_types
cimport pyne.cpp_utils
cimport pyne.pyne_config
import pyne.pyne_config
cimport cpp_nucname
cimport cpp_particle
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv


warn(__name__ + " is not yet V&V compliant.", VnVWarning)

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

def pdc_number(x):
    """pdc_number(x)
    
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
        cn = cpp_particle.pdc_number(std_string(<char *> x_bytes))
    elif isinstance(x, int):
        cn = cpp_particle.pdc_number(<int> long(x))

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


