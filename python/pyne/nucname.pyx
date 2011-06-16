"""Python wrapper for nucname library."""
# Python imports 
#from collections import Iterable

# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
#from cython cimport pointer

# local imports 
cimport std
cimport cpp_nucname
cimport pyne.stlconverters as conv
import pyne.stlconverters as conv

#
# Conversion dictionaries
#

LLzz = conv.map_to_dict_str_int(cpp_nucname.LLzz)
zzLL = conv.map_to_dict_int_str(cpp_nucname.zzLL)


#
# Elemental string sets
#

# Issues with deallocating using the below :(
#LAN = conv.SetProxyStr()
#(<conv._SetProxyStr> LAN).init(&cpp_nucname.LAN)

cdef conv._SetProxyStr LAN_proxy = conv.SetProxyStr()
LAN_proxy.init(&cpp_nucname.LAN)
LAN = LAN_proxy

cdef conv._SetProxyStr ACT_proxy = conv.SetProxyStr()
ACT_proxy.init(&cpp_nucname.ACT)
ACT = ACT_proxy

cdef conv._SetProxyStr TRU_proxy = conv.SetProxyStr()
TRU_proxy.init(&cpp_nucname.TRU)
TRU = TRU_proxy

cdef conv._SetProxyStr MA_proxy = conv.SetProxyStr()
MA_proxy.init(&cpp_nucname.MA)
MA = MA_proxy

cdef conv._SetProxyStr FP_proxy = conv.SetProxyStr()
FP_proxy.init(&cpp_nucname.FP)
FP = FP_proxy


#
# Elemental integer sets
#

# Issues with deallocating using the below :(
#lan = conv.SetProxyInt()
#(<conv._SetProxyInt> lan).init(&cpp_nucname.lan)

cdef conv._SetProxyInt lan_proxy = conv.SetProxyInt()
lan_proxy.init(&cpp_nucname.lan)
lan = lan_proxy

cdef conv._SetProxyInt act_proxy = conv.SetProxyInt()
act_proxy.init(&cpp_nucname.act)
act = act_proxy

cdef conv._SetProxyInt tru_proxy = conv.SetProxyInt()
tru_proxy.init(&cpp_nucname.tru)
tru = tru_proxy

cdef conv._SetProxyInt ma_proxy = conv.SetProxyInt()
ma_proxy.init(&cpp_nucname.ma)
ma = ma_proxy

cdef conv._SetProxyInt fp_proxy = conv.SetProxyInt()
fp_proxy.init(&cpp_nucname.fp)
fp = fp_proxy




class NucTypeError(Exception):
    def __init__(self, nuc=None):
        self.nuc = nuc

    def __str__(self):
        msg = "Nuclide type not an int or str"
        if self.nuc is not None:
            msg += ": " + repr(self.nuc) 
        return msg

#
# Current Form Function
#

def current_form(nuc):
    """Find the current form of a nuclide.

    Parameters
    ----------
    nuc : int or str 
        Input nuclide(s).

    Returns
    -------
    form_flag : str
        The form identifier string from ["zzaaam", "LLAAAM", "MCNP"].
    """
    cdef std.string cpp_curr_form 

    if isinstance(nuc, basestring):
        cpp_curr_form = cpp_nucname.current_form(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        cpp_curr_form = cpp_nucname.current_form(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return cpp_curr_form.c_str()


#
# zzaaam Functions
#

def zzaaam(nuc):
    """Converts a nuclide to its zzaaam form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : int 
        Output nuclide in zzaaam form.
    """

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.zzaaam(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        newnuc = cpp_nucname.zzaaam(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc


def LLAAAM(nuc):
    """Converts a nuclide to its LLAAAM form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : str 
        Output nuclide in LLAAAM form.
    """
    cdef std.string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.LLAAAM(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.LLAAAM(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc.c_str()


def mcnp(nuc):
    """Converts a nuclide to its MCNP form (int). 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : int 
        Output nuclide in MCNP form.
    """

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.mcnp(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.mcnp(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc


def serpent(nuc):
    """Converts a nuclide to its Serepnt form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : str 
        Output nuclide in serpent form.
    """
    cdef std.string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.serpent(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.serpent(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc.c_str()


#
# Helper Functions
#

def nuc_weight(nuc):
    """Calculates the weight of a nuclide in [amu].

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    weight : float
        Atomic weight of this nuclide [amu].
    """
    if isinstance(nuc, basestring):
        weight = cpp_nucname.nuc_weight(<char *> nuc)
    elif isinstance(nuc, int):
        weight = cpp_nucname.nuc_weight(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return weight


