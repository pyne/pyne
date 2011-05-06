"""Python wrapper for nucname library."""
# Python imports 
from collections import Iterable

# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std
cimport cpp_nucname
cimport stlconverters as conv

#
# Conversion dictionaries
#

LLzz = conv.map_to_dict_str_int(cpp_nucname.LLzz)
zzLL = conv.map_to_dict_int_str(cpp_nucname.zzLL)


#
# Elemental string sets
#

LAN = conv.cpp_to_py_set_str(cpp_nucname.LAN)
ACT = conv.cpp_to_py_set_str(cpp_nucname.ACT)
TRU = conv.cpp_to_py_set_str(cpp_nucname.TRU)
MA = conv.cpp_to_py_set_str(cpp_nucname.MA)
FP = conv.cpp_to_py_set_str(cpp_nucname.FP)


#
# Elemental integer sets
#

lan = conv.cpp_to_py_set_int(cpp_nucname.lan)
act = conv.cpp_to_py_set_int(cpp_nucname.act)
tru = conv.cpp_to_py_set_int(cpp_nucname.tru)
ma = conv.cpp_to_py_set_int(cpp_nucname.ma)
fp = conv.cpp_to_py_set_int(cpp_nucname.fp)


class NucTypeError(Exception):
    def __init__(self, nuc=None):
        self.nuc = nuc

    def __str__(self):
        msg = "Nuclide type not an int, str, or iterable"
        if self.nuc is not None:
            msg += ": " + repr(self.nuc) 
        return msg

#
# Current Form Function
#

def current_form(nuc):
    """Find the current form of a nuclide.

    Args:
        * nuc (int, str, or iterable): Input nuclide(s).

    Returns:
        * form_flag (str or iterable): The form identifier string from ["zzaaam", "LLAAAM", "MCNP"].
    """
    cdef std.string cpp_curr_form 

    if isinstance(nuc, basestring):
        cpp_curr_form = cpp_nucname.current_form(nuc)
        form_flag = cpp_current_form.c_str()
    elif isinstance(nuc, int) or isinstance(nuc, long):
        cpp_curr_form = cpp_nucname.current_form(<int> nuc)
        form_flag = cpp_current_form.c_str()
    elif isinstance(nuc, Iterable):
        try:
            form_flag = nuc.__class__(current_form(n) for n in nuc)
        except TypeError:
            form_flag = [current_form(n) for n in nuc]
    else:
        raise NucTypeError(nuc)

#    return cpp_current_form.c_str()
    return form_flag


#
# zzaaam Functions
#

def zzaaam(nuc):
    """Converts a nuclide to its zzaaam (int) form. 

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * newnuc (int): Output nuclide in zzaaam form.
    """

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.zzaaam(std.string(nuc))
    elif isinstance(nuc, int) or isinstance(nuc, long):
        newnuc = cpp_nucname.zzaaam(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc


def LLAAAM(nuc):
    """Converts an arbitrary nuclide and its LLAAAM form. 

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * newnuc (int): Output nuclide in LLAAAM form.
    """
    cdef std.string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.LLAAAM(std.string(nuc))
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.LLAAAM(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc.c_str()


def mcnp(nuc):
    """Converts an arbitrary nuclide and its MCNP form. 

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * newnuc (int): Output nuclide in MCNP form.
    """

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.mcnp(std.string(nuc))
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.mcnp(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc


#
# Helper Functions
#

def nuc_weight(nuc):
    """Calculates the weight of a nuclide in [amu].

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * weight (float): Atomic weight of this nuclide [amu].
    """
    if isinstance(nuc, basestring):
        weight = cpp_nucname.nuc_weight(std.string(nuc))
    elif isinstance(nuc, int):
        weight = cpp_nucname.nuc_weight(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return weight



#
# isovec_keys_2_* Functions
#

def isovec_keys_2_zzaaam(dict isovec):
    """Converts all keys of an isotopic vector dictionary to zzaaam form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of zzaaam (int) form.
    """
#    newvec = {}
#
#    for iso in isovec.keys():
#        newvec[mixed_2_zzaaam(iso)] = isovec[iso]
#
#    return newvec


def isovec_keys_2_LLAAAM(dict isovec):
    """Converts all keys of an isotopic vector dictionary to LLAAAM form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of LLAAAM (str) form.
    """
#    newvec = {}
#
#    for iso in isovec.keys():
#        newvec[mixed_2_LLAAAM(iso)] = isovec[iso]
#
#    return newvec


def isovec_keys_2_MCNP(dict isovec):
    """Converts all keys of an isotopic vector dictionary to MCNP form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of MCNP (int) form.
    """
#    newvec = {}

#    for iso in isovec.keys():
#        newvec[mixed_2_MCNP(iso)] = isovec[iso]
#
#    return newvec

