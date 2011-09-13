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
cimport pyne.cpp_pyne
cimport pyne.pyne_config
import pyne.pyne_config

cimport cpp_nucname
cimport pyne.stlconverters as conv
import pyne.stlconverters as conv

#
# Conversion dictionaries
#

cdef conv._MapStrInt name_zz_proxy = conv.MapStrInt(False)
name_zz_proxy.map_ptr = &cpp_nucname.name_zz
name_zz = name_zz_proxy

cdef conv._MapIntStr zzname_proxy = conv.MapIntStr(False)
zzname_proxy.map_ptr = &cpp_nucname.zz_name
zz_name = zzname_proxy


#
# Elemental string sets
#

cdef conv._SetStr LAN_proxy = conv.SetStr(False)
LAN_proxy.set_ptr = &cpp_nucname.LAN
LAN = LAN_proxy

cdef conv._SetStr ACT_proxy = conv.SetStr(False)
ACT_proxy.set_ptr = &cpp_nucname.ACT
ACT = ACT_proxy

cdef conv._SetStr TRU_proxy = conv.SetStr(False)
TRU_proxy.set_ptr = &cpp_nucname.TRU
TRU = TRU_proxy

cdef conv._SetStr MA_proxy = conv.SetStr(False)
MA_proxy.set_ptr = &cpp_nucname.MA
MA = MA_proxy

cdef conv._SetStr FP_proxy = conv.SetStr(False)
FP_proxy.set_ptr = &cpp_nucname.FP
FP = FP_proxy


#
# Elemental integer sets
#

cdef conv._SetInt lan_proxy = conv.SetInt(False)
lan_proxy.set_ptr = &cpp_nucname.lan
lan = lan_proxy

cdef conv._SetInt act_proxy = conv.SetInt(False)
act_proxy.set_ptr = &cpp_nucname.act
act = act_proxy

cdef conv._SetInt tru_proxy = conv.SetInt(False)
tru_proxy.set_ptr = &cpp_nucname.tru
tru = tru_proxy

cdef conv._SetInt ma_proxy = conv.SetInt(False)
ma_proxy.set_ptr = &cpp_nucname.ma
ma = ma_proxy

cdef conv._SetInt fp_proxy = conv.SetInt(False)
fp_proxy.set_ptr = &cpp_nucname.fp
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
        The form identifier string from ["zzaaam", "name", "MCNP"].
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
    """Converts a nuclide to its zzaaam form (952420). 

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


def name(nuc):
    """Converts a nuclide to its name form ('AM242M'). 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : str 
        Output nuclide in name form.
    """
    cdef std.string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.name(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.name(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc.c_str()


def mcnp(nuc):
    """Converts a nuclide to its MCNP form (95642). 

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
    """Converts a nuclide to its Serepnt form ('Am-242m'). 

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



def nist(nuc):
    """Converts a nuclide to NIST form ('242Am').

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : str 
        Output nuclide in nist form.
    """
    cdef std.string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.nist(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.nist(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc.c_str()



#
# C++ Helper Functions
#

cdef cpp_set[int] zzaaam_set(object nuc_sequence):
    cdef int nuc_zz
    cdef cpp_set[int] nuc_set = cpp_set[int]()

    for nuc in nuc_sequence:
        nuc_zz = zzaaam(nuc)
        nuc_set.insert(nuc_zz)

    return nuc_set


