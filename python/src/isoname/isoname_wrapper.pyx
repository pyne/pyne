"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std
cimport cpp_isoname
cimport stlconverters as conv

#
# Conversion dictionaries
#

LLzz = conv.map_to_dict_str_int(cpp_isoname.LLzz)
zzLL = conv.map_to_dict_int_str(cpp_isoname.zzLL)


#
# Elemental string sets
#

LAN = conv.cpp_to_py_set_str(cpp_isoname.LAN)
ACT = conv.cpp_to_py_set_str(cpp_isoname.ACT)
TRU = conv.cpp_to_py_set_str(cpp_isoname.TRU)
MA = conv.cpp_to_py_set_str(cpp_isoname.MA)
FP = conv.cpp_to_py_set_str(cpp_isoname.FP)


#
# Elemental integer sets
#

lan = conv.cpp_to_py_set_int(cpp_isoname.lan)
act = conv.cpp_to_py_set_int(cpp_isoname.act)
tru = conv.cpp_to_py_set_int(cpp_isoname.tru)
ma = conv.cpp_to_py_set_int(cpp_isoname.ma)
fp = conv.cpp_to_py_set_int(cpp_isoname.fp)


#
# Current Form Function
#

def CurrentForm(nuc):
    """Find the current form of a nuclide.

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * form_flag (str): The form identifier string from ["zzaaam", "LLAAAM", "MCNP"].
    """
    cdef std.string cpp_CurrentForm 

    if isinstance(nuc, basestring):
        cpp_CurrentForm = cpp_isoname.CurrentForm(std.string(nuc))
    elif isinstance(nuc, int):
        cpp_CurrentForm = cpp_isoname.CurrentForm(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return cpp_CurrentForm.c_str()


#
# LLAAAM_2_* Functions
#

def LLAAAM_2_zzaaam(char * nuc):
    """Converts a nuclide from LLAAAM form to its zzaaam form. 

    Args:
        * nuc (str): Input nuclide in LLAAAM form.

    Returns:
        * newnuc (int): Output nuclide in zzaaam form.
    """
    return cpp_isoname.LLAAAM_2_zzaaam(std.string(nuc))


def LLAAAM_2_MCNP(char * nuc):
    """Converts a nuclide from LLAAAM form to its MCNP form. 

    Args:
        * nuc (str): Input nuclide in LLAAAM form.

    Returns:
        * newnuc (int): Output nuclide in MCNP form.
    """
    return cpp_isoname.LLAAAM_2_MCNP(std.string(nuc))

#
# zzaaam_2_* Functions 
#

def zzaaam_2_LLAAAM(int nuc):
    """Converts a nuclide from zzaaam form to its LLAAAM form. 

    Args:
        * nuc (str): Input nuclide in zzaaam form.

    Returns:
        * newnuc (int): Output nuclide in LLAAAM form.
    """
    cdef std.string cpp_LLAAAM = cpp_isoname.zzaaam_2_LLAAAM(nuc)
    return cpp_LLAAAM.c_str()


def zzaaam_2_MCNP(int nuc):
    """Converts a nuclide from zzaaam form to its MCNP form. 

    Args:
        * nuc (str): Input nuclide in zzaaam form.

    Returns:
        * newnuc (int): Output nuclide in MCNP form.
    """
    return cpp_isoname.zzaaam_2_MCNP(nuc)


#
# MCNP_2_* Functions
#

def MCNP_2_zzaaam(int nuc):
    """Converts a nuclide from MCNP form to its zzaaam form. 

    Args:
        * nuc (str): Input nuclide in MCNP form.

    Returns:
        * newnuc (int): Output nuclide in zzaaam form.
    """
    return cpp_isoname.MCNP_2_zzaaam(nuc)


def MCNP_2_LLAAAM(int nuc):
    """Converts a nuclide from MCNP form to its LLAAAM form. 

    Args:
        * nuc (str): Input nuclide in MCNP form.

    Returns:
        * newnuc (int): Output nuclide in LLAAAM form.
    """
    cdef std.string cpp_LLAAAM = cpp_isoname.MCNP_2_LLAAAM(nuc)
    return cpp_LLAAAM.c_str()


#
# mixed_2_* Functions
#

def mixed_2_zzaaam(nuc):
    """Converts an arbitrary nuclide and its zzaaam form. 

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * newnuc (int): Output nuclide in zzaaam form.
    """

    if isinstance(nuc, basestring):
        newnuc = cpp_isoname.mixed_2_zzaaam(std.string(nuc))
    elif isinstance(nuc, int):
        newnuc = cpp_isoname.mixed_2_zzaaam(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return newnuc


def mixed_2_LLAAAM(nuc):
    """Converts an arbitrary nuclide and its LLAAAM form. 

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * newnuc (int): Output nuclide in LLAAAM form.
    """
    cdef std.string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_isoname.mixed_2_LLAAAM(std.string(nuc))
    elif isinstance(nuc, int):
        newnuc = cpp_isoname.mixed_2_LLAAAM(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return newnuc.c_str()


def mixed_2_MCNP(nuc):
    """Converts an arbitrary nuclide and its MCNP form. 

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * newnuc (int): Output nuclide in MCNP form.
    """

    if isinstance(nuc, basestring):
        newnuc = cpp_isoname.mixed_2_MCNP(std.string(nuc))
    elif isinstance(nuc, int):
        newnuc = cpp_isoname.mixed_2_MCNP(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return newnuc


#
# Helper Functions
#

def nuc_weight_zzaaam(int nuc):
    """Calculates the weight of a nuclide in [amu].

    Args:
        * nuc (int): Input nuclide.

    Returns:
        * weight (float): Atomic weight of this nuclide [amu].
    """
    return cpp_isoname.nuc_weight_zzaaam(nuc)


def nuc_weight(nuc):
    """Calculates the weight of a nuclide in [amu].

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * weight (float): Atomic weight of this nuclide [amu].
    """
    if isinstance(nuc, basestring):
        weight = cpp_isoname.nuc_weight(std.string(nuc))
    elif isinstance(nuc, int):
        weight = cpp_isoname.nuc_weight(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return weight


#
# (a*)_2_(b*)_List Functions
#

def LLAAAM_2_zzaaam_List(list nuclist):
    """Converts a list of LLAAAM form to a list of zzaaam form.

    Args:
        * `nuclist` (str list): List of LLAAAM nuclides.

    Returns:
        * `newnuclist` (int list): List of zzaaam nuclides.
    """
    return [LLAAAM_2_zzaaam(nuc) for nuc in nuclist]


def LLAAAM_2_MCNP_List(list nuclist):
    """Converts a list of LLAAAM form to a list of MCNP form.

    Args:
        * `nuclist` (str list): List of LLAAAM nuclides.

    Returns:
        * `newnuclist` (int list): List of MCNP nuclides.
    """
    return [LLAAAM_2_MCNP(nuc) for nuc in nuclist]


def zzaaam_2_LLAAAM_List(list nuclist):
    """Converts a list of zzaaam form to a list of LLAAAM form.

    Args:
        * `nuclist` (int list): List of zzaaam nuclides.

    Returns:
        * `newnuclist` (str list): List of LLAAAM nuclides.
    """
    return [zzaaam_2_LLAAAM(nuc) for nuc in nuclist]


def zzaaam_2_MCNP_List(list nuclist):
    """Converts a list of zzaaam form to a list of MCNP form.

    Args:
        * `nuclist` (int list): List of zzaaam nuclides.

    Returns:
        * `newnuclist` (int list): List of MCNP nuclides.
    """
    return [zzaaam_2_MCNP(nuc) for nuc in nuclist]


def MCNP_2_LLAAAM_List(list nuclist):
    """Converts a list of MCNP form to a list of LLAAAM form.

    Args:
        * `nuclist` (int list): List of MCNP nuclides.

    Returns:
        * `newnuclist` (str list): List of LLAAAM nuclides.
    """
    return [MCNP_2_LLAAAM(nuc) for nuc in nuclist]


def MCNP_2_zzaaam_List(list nuclist):
    """Converts a list of MCNP form to a list of zzaaam form.

    Args:
        * `nuclist` (int list): List of MCNP nuclides.

    Returns:
        * `newnuclist` (int list): List of zzaaam nuclides.
    """
    return [MCNP_2_zzaaam(nuc) for nuc in nuclist]


#
# mixed_2_*_List Functions
#

def RearRemoveDuplicates(list l):
    """Removes duplicate entries from list l, starting from the back. 
    Used internally in the [form]_2_[form]_List() functions.

    Args:
       * `l` (list): input list.

    Returns:
       * `l` (list): input with duplicates removed.
    """

    for n in range(len(l)-1, -1, -1):
        if 1 < l.count(l[n]):
            l.pop(n)
    return l
    

def mixed_2_zzaaam_List(list nuclist):
    """Converts a list of mixed form to a list of zzaaam form.

    Args:
        * `nuclist` (str or int list): List of nuclides of mixed form.

    Returns:
        * `newnuclist` (int list): List of zzaaam nuclides.
    """
    return RearRemoveDuplicates( [mixed_2_zzaaam(nuc) for nuc in nuclist] )


def mixed_2_LLAAAM_List(list nuclist):
    """Converts a list of mixed form to a list of LLAAAM form.

    Args:
        * `nuclist` (str or int list): List of nuclides of mixed form.

    Returns:
        * `newnuclist` (str list): List of LLAAAM nuclides.
    """
    return RearRemoveDuplicates( [mixed_2_LLAAAM(nuc) for nuc in nuclist] )


def mixed_2_MCNP_List(list nuclist):
    """Converts a list of mixed form to a list of MCNP form.

    Args:
        * `nuclist` (str or int list): List of nuclides of mixed form.

    Returns:
        * `newnuclist` (int list): List of MCNP nuclides.
    """
    return RearRemoveDuplicates( [mixed_2_MCNP(nuc) for nuc in nuclist] )

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
    newvec = {}

    for iso in isovec.keys():
        newvec[mixed_2_zzaaam(iso)] = isovec[iso]

    return newvec


def isovec_keys_2_LLAAAM(dict isovec):
    """Converts all keys of an isotopic vector dictionary to LLAAAM form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of LLAAAM (str) form.
    """
    newvec = {}

    for iso in isovec.keys():
        newvec[mixed_2_LLAAAM(iso)] = isovec[iso]

    return newvec


def isovec_keys_2_MCNP(dict isovec):
    """Converts all keys of an isotopic vector dictionary to MCNP form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of MCNP (int) form.
    """
    newvec = {}

    for iso in isovec.keys():
        newvec[mixed_2_MCNP(iso)] = isovec[iso]

    return newvec

