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
include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string
cimport pyne.cpp_pyne
cimport pyne.pyne_config
import pyne.pyne_config

from pyne cimport cpp_nucname
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv

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
# Is Nuclide Function
#

def isnuclide(nuc):
    """Test if nuc is a valid nuclide.

    Parameters
    ----------
    nuc : int or str 
        Input nuclide(s).

    Returns
    -------
    flag : bool

    """
    if isinstance(nuc, basestring):
        flag = cpp_nucname.isnuclide(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        flag = cpp_nucname.isnuclide(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return flag



def id(nuc):
    """Converts a nuclide to its identifier form (952420000).

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : int
        Output nuclide id.

    """
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.id(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        newnuc = cpp_nucname.id(<int> nuc)
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
    cdef std_string newnuc
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.name(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.name(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return <char *> newnuc.c_str()


def znum(nuc):
    """Retrieves a nuclide's charge number (95).

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    z : int
        The number of protons in the nucleus.

    """
    if isinstance(nuc, basestring):
        z = cpp_nucname.znum(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        z = cpp_nucname.znum(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return z


def anum(nuc):
    """Retrieves a nuclide's nucleon number (95).

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    a : int
        The number of protons and neutrons in the nucleus.

    """
    if isinstance(nuc, basestring):
        a = cpp_nucname.anum(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        a = cpp_nucname.anum(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return a


def snum(nuc):
    """Retrieves a nuclide's excitation number (95).

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    s : int
        The excitation level the nucleus.

    """
    if isinstance(nuc, basestring):
        s = cpp_nucname.snum(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        s = cpp_nucname.snum(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return s


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


def zzaaam_to_id(nuc):
    """Converts a nuclide directly from ZZAAAM form (952420) to
    the canonical identifier form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide in ZZAAAM form.

    Returns
    -------
    newnuc : int 
        Output nuclide in identifier form.

    """
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.zzaaam_to_id(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        newnuc = cpp_nucname.zzaaam_to_id(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return newnuc


def mcnp(nuc):
    """Converts a nuclide to its MCNP form (92636). 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : int 
        Output nuclide in MCNP form.

    Notes
    -----
    Most metastables in this form add 300 + 100*m where 
    m is the isomeric state (U-236m = 92636).  However,
    MCNP special cases Am-242 and Am-242m by switching 
    the meaning. Thus Am-242m = 95242 and Am-242 = 95642.

    """

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.mcnp(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.mcnp(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc


def mcnp_to_id(nuc):
    """Converts a nuclide directly from MCNP form (92636) to
    the canonical identifier form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide in MCNP form.

    Returns
    -------
    newnuc : int 
        Output nuclide in identifier form.

    """
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.mcnp_to_id(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        newnuc = cpp_nucname.mcnp_to_id(<int> nuc)
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
    cdef std_string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.serpent(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.serpent(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return <char *> newnuc.c_str()


def serpent_to_id(nuc):
    """Converts a nuclide directly from Serpent form ('Am-242m') to
    the canonical identifier form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide in Serpent form.

    Returns
    -------
    newnuc : int 
        Output nuclide in identifier form.

    """
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.serpent_to_id(<char *> nuc)
    #elif isinstance(nuc, int) or isinstance(nuc, long):
    #    newnuc = cpp_nucname.serpent_to_id(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return newnuc


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
    cdef std_string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.nist(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.nist(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return <char *> newnuc.c_str()


def nist_to_id(nuc):
    """Converts a nuclide directly from NIST form ('242Am') to
    the canonical identifier form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide in NIST form.

    Returns
    -------
    newnuc : int 
        Output nuclide in identifier form.

    """
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.nist_to_id(<char *> nuc)
    #elif isinstance(nuc, int) or isinstance(nuc, long):
    #    newnuc = cpp_nucname.nist_to_id(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return newnuc


def cinder(nuc):
    """Converts a nuclide to its CINDER (aaazzzm) form (2420951). 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : int 
        Output nuclide in CINDER (aaazzzm) form.

    """

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.cinder(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.cinder(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return newnuc


def cinder_to_id(nuc):
    """Converts a nuclide directly from Cinder form (2420951) to
    the canonical identifier form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide in Cinder form.

    Returns
    -------
    newnuc : int 
        Output nuclide in identifier form.

    """
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.cinder_to_id(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        newnuc = cpp_nucname.cinder_to_id(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return newnuc


def alara(nuc):
    """Converts a nuclide to its ALARA form ('am:242'). 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : str 
        Output nuclide in name form.

    """
    cdef std_string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.alara(<char *> nuc)
    elif isinstance(nuc, int):
        newnuc = cpp_nucname.alara(<int> nuc)
    else:
        raise NucTypeError(nuc)

    return <char *> newnuc.c_str()


def alara_to_id(nuc):
    """Converts a nuclide directly from ALARA form ('am:242') to
    the canonical identifier form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide in ALARA form.

    Returns
    -------
    newnuc : int 
        Output nuclide in identifier form.

    """
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.alara_to_id(<char *> nuc)
    #elif isinstance(nuc, int) or isinstance(nuc, long):
    #    newnuc = cpp_nucname.alara_to_id(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return newnuc


def sza(nuc):
    """Converts a nuclide to its SZA form (SSSZZZAAA). 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    newnuc : int 
        Output nuclide in SZA form.

    """
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.sza(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        newnuc = cpp_nucname.sza(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return newnuc


def sza_to_id(nuc):
    """Converts a nuclide directly from SZA form (SSSZZZAAA) to
    the canonical identifier form. 

    Parameters
    ----------
    nuc : int or str 
        Input nuclide in SZA form.

    Returns
    -------
    newnuc : int 
        Output nuclide in identifier form.

    """
    if isinstance(nuc, basestring):
        newnuc = cpp_nucname.sza_to_id(<char *> nuc)
    elif isinstance(nuc, int) or isinstance(nuc, long):
        newnuc = cpp_nucname.sza_to_id(<int> nuc)
    else:
        raise NucTypeError(nuc)
    return newnuc


#
# C++ Helper Functions
#

cdef cpp_set[int] id_set(object nuc_sequence):
    cdef int nucid
    cdef cpp_set[int] nuc_set = cpp_set[int]()
    for nuc in nuc_sequence:
        nucid = id(nuc)
        nuc_set.insert(nucid)
    return nuc_set

cdef cpp_set[int] zzaaam_set(object nuc_sequence):
    cdef int nuc_zz
    cdef cpp_set[int] nuc_set = cpp_set[int]()
    for nuc in nuc_sequence:
        nuc_zz = zzaaam(nuc)
        nuc_set.insert(nuc_zz)
    return nuc_set
