"""Python wrapper for rxname library.

The rxname module provides a string-based canonical form for reaction
channel names, an integer identifier, labeling and documentation, and
conversion to/from other conventions for these names.

    **Reaction Names:** The names themselves are strings chosen such
    that they are valid variable names in most programming languages
    (including Python and C/C++).  This strategy is known as *natural
    naming* and enables a *namespace* of reactions.  Therefore, all
    names must match the regular expression ``[A-Za-z_][A-Za-z0-9_]*``.
    For example, the elastic scattering cross section is simply
    "elastic" while the pair production reaction is given by
    "pair_prod".

    A number of patterns dictate how a reaction should be named.
    Foremost among these are particle names.  Where required, "z" is a
    variable for any incident particle.  The following table displays
    particle types and their names:

        =========  =========
        particle   name (z)
        =========  =========
        neutron    n
        proton     p
        deuterium  d
        tritium    t
        Helium-3   He3
        alpha      a
        gamma      gamma
        =========  =========

    From this we can see that a reaction which produces a neutron and a
    proton is called "np".  If multiple particles of the same type are
    produced, then the number precedes the particle type.  Thus, one
    neutron and two protons are given by "n2p".  However if this would
    result in the name starting with a number, then the name is
    prepended with ``z_`` to indicate any incident particle.  For
    example, the reaction which yields two neutrons is "z_2n" (because
    "2n" is not a valid variable name in most programming languages).

    Furthermore, if a reaction name ends in ``_[0-9]+`` (underscore plus
    digits), then this means that the nucleus is left in the nth excited
    state after the interaction.  For example, "n_0" produces a neutron
    and leaves the nucleus in the ground state, "n_1" produces a neutron
    and the nucleus is in the first excited state, and so on.  However,
    "_continuum" means that the nucleus in an energy state in the
    continuum.

    If a reaction name begins with ``erel_``, then this channel is for
    the energy release from the reaction by the name without ``erel_``.
    E.g. "erel_p" is the energy release from proton emission.

    **Reaction IDs:** While the reaction names are sufficient for
    defining all possible reactions in a partially physically meaningful
    way, they do so using a variable length format (strings).  It is
    often advantageous to have a fixed-width format, namely for storage.
    To this end, unique unsigned 32-bit integers are given to each name.
    These identifiers are computed based on a custom hash of the
    reaction name.  This hash function reserves space for MT numbers by
    not producing values below 1000.  It is recommended that the
    reaction identifiers be used for most performance critical tasks,
    rather than the names that they are calculated from.

    **Reaction Labels:** Reaction labels are short, human-readable
    strings.  These do not follow the naming convention restrictions
    that the names themselves are subject to.  Additionally, labels need
    not be unique.  The labels are provided as a convenience for
    building user interfaces with.

    **Reaction Docstrings:** Similar to labels, reactions also come with
    a documentation string which gives a description of the reaction in
    a sentence or two. These provide more help and information for a
    reaction.  This may be useful in a tool tip context for user
    interfaces.

    **Other Canonical Forms:** This module provides mappings between
    other reaction canonical forms and the naming conventions and IDs
    used here.  The most widespread of these are arguably the MT
    numbers.  MT numbers are a strict subset of the reactions used here.
    Further information may be found at [NNDC]_, [NEA]_, [T2]_, and
    [JAEA]_.

The rxname module implements a suite of functions for computing or
retrieving reaction names and their associated data described above.
These functions have a variety of interfaces. Lookup may occur either by
name, ID, MT number, a string of ID, or a string of MT number.

However, lookup may also occur via alternate names or abbreviations.
For example, "tot" and "abs" will give the names "total" and
"absorption".  Spelling out particles will also work; "alpha" and "duet"
will give "a" and "d".  For a listing of all alternative names see the
``altnames`` variable.

Furthermore, certain reactions may be inferred from the nuclide prior
and post reaction.  For example, if an incident neutron hit U-235 and
Th-232 was produced then an alpha production reaction is assumed to have
occurred. Thus most of the functions in rxname will take a from nuclide,
a to nuclide, and z -- the incident particle type (which defaults to "n"
neutron).  Note that z may also be "decay", indicating a radioactive
decay occurrence.

-------------------------
Adding New Reaction Names
-------------------------

If for some reason the existing reaction names are insufficient, the
following procedure will add a new reaction to the suite provided.

1. Increment ``NUM_RX_NAMES`` in "rxname.h".
2. Add the name to the ``_names`` array at the top of "rxname.cpp".
   Note the location in this array that the new name was added.
3. In the ``_fill_maps()`` function in "rxname.cpp", add the MT number
   to the ``_mts`` array at the same index the name was added in step 2.
   If the reaction does not have a corresponding MT number, add a zero
   (0) at this location instead.
4. In the ``_fill_maps()`` function in "rxname.cpp", add a short label
   to the ``_labels`` array at the same index the name was added in step
   2.
5. In the ``_fill_maps()`` function in "rxname.cpp", add a docstring to
   the ``_docs`` array at the same index the name was added in step 2.

Repeat this procedure as necessary.

--------------------------

.. [NNDC] http://www.nndc.bnl.gov/endfdocs/ENDF-102-2001.pdf
.. [NEA] http://www.oecd-nea.org/dbdata/data/manual-endf/endf102_MT.pdf
.. [T2] http://t2.lanl.gov/nis/endf/mts.html
.. [JAEA] http://wwwndc.jaea.go.jp/form/ENDF6/mt.html

--------------------------

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
cimport cpp_rxname
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv



QA_warn(__name__)

# names
cdef conv._SetStr names_proxy = conv.SetStr(False)
names_proxy.set_ptr = &cpp_rxname.names
names = names_proxy

# altnames
cdef conv._MapStrUInt altnames_proxy = conv.MapStrUInt(False)
altnames_proxy.map_ptr = &cpp_rxname.altnames
altnames = altnames_proxy

# id_name
cdef conv._MapUIntStr id_name_proxy = conv.MapUIntStr(False)
id_name_proxy.map_ptr = &cpp_rxname.id_name
id_name = id_name_proxy

# name_id
cdef conv._MapStrUInt name_id_proxy = conv.MapStrUInt(False)
name_id_proxy.map_ptr = &cpp_rxname.name_id
name_id = name_id_proxy

# id_mt
cdef conv._MapUIntUInt id_mt_proxy = conv.MapUIntUInt(False)
id_mt_proxy.map_ptr = &cpp_rxname.id_mt
id_mt = id_mt_proxy

# mt_id
cdef conv._MapUIntUInt mt_id_proxy = conv.MapUIntUInt(False)
mt_id_proxy.map_ptr = &cpp_rxname.mt_id
mt_id = mt_id_proxy

# labels
cdef conv._MapUIntStr labels_proxy = conv.MapUIntStr(False)
labels_proxy.map_ptr = &cpp_rxname.labels
labels = labels_proxy

# docs
cdef conv._MapUIntStr docs_proxy = conv.MapUIntStr(False)
docs_proxy.map_ptr = &cpp_rxname.docs
docs = docs_proxy



def hash(s):
    """hash(s)

    Hashes a string to be used as a reaction id.  This hash specifically reserves
    the h < 1000 for MT numbers.
    """
    s_bytes = s.encode()
    return int(cpp_rxname.hash(<const_char *> s_bytes))


def name(x, y=None, z="n"):
    """name(x, y=None, z="n")

    Gives the unique reaction name.  Note that this name follows the 'natural naming'
    convention and may be used as a variable in most languages.

    Parameters
    ----------
    x : str, int, or long
        name, abbreviation, id, MT number, or from nuclide.
    y : str, int, or None, optional
        to nuclide.
    z : str, optional
        incident particle type ("n", "p", ...) when x and y are nuclides.

    Returns
    -------
    n : str
        a unique reaction name.
    """
    cdef std_string cn
    cdef int from_nuc, to_nuc
    if y is None:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            cn = cpp_rxname.name(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            cn = cpp_rxname.name(<extra_types.uint32> long(x))
        elif isinstance(x, long):
            cn = cpp_rxname.name(<extra_types.uint32> x)
    else:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            from_nuc = cpp_nucname.id(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            from_nuc = cpp_nucname.id(<int> x)
        if isinstance(y, basestring):
            y_bytes = y.encode()
            to_nuc = cpp_nucname.id(std_string(<char *> y_bytes))
        elif isinstance(y, int):
            to_nuc = cpp_nucname.id(<int> y)
        z_bytes = z.encode()
        cn = cpp_rxname.name(from_nuc, to_nuc, std_string(<char *> z_bytes))
    n = bytes(<char *> cn.c_str()).decode()
    return n


def id(x, y=None, z="n"):
    """id(x, y=None, z="n")

    Gives the unique reaction id.  This is originally calculated as the hash of the
    name.  This id may not be less than 1000, since these integers are reserved for
    MT numbers.

    Parameters
    ----------
    x : str, int, or long
        name, abbreviation, id, MT number, or from nuclide.
    y : str, int, or None, optional
        to nuclide.
    z : str, optional
        incident particle type ("n", "p", ...) when x and y are nuclides.

    Returns
    -------
    rxid : int or long
        a unique reaction identifier.
    """
    cdef int from_nuc, to_nuc
    if y is None:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            rxid = cpp_rxname.id(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            rxid = cpp_rxname.id(<extra_types.uint32> long(x))
        elif isinstance(x, long):
            rxid = cpp_rxname.id(<extra_types.uint32> x)
    else:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            from_nuc = cpp_nucname.id(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            from_nuc = cpp_nucname.id(<int> x)
        if isinstance(y, basestring):
            y_bytes = y.encode()
            to_nuc = cpp_nucname.id(std_string(<char *> y_bytes))
        elif isinstance(y, int):
            to_nuc = cpp_nucname.id(<int> y)
        z_bytes = z.encode()
        rxid = cpp_rxname.id(from_nuc, to_nuc, std_string(<char *> z_bytes))
    return int(rxid)


def mt(x, y=None, z="n"):
    """mt(x, y=None, z="n")

    Gives the reaction MT number. This may not be greater than 1000.

    Parameters
    ----------
    x : str, int, or long
        name, abbreviation, id, MT number, or from nuclide.
    y : str, int, or None, optional
        to nuclide.
    z : str, optional
        incident particle type ("n", "p", ...) when x and y are nuclides.

    Returns
    -------
    mtnum : int or long
        a unique MT number.
    """
    cdef int from_nuc, to_nuc
    if y is None:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            mtnum = cpp_rxname.mt(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            mtnum = cpp_rxname.mt(<extra_types.uint32> long(x))
        elif isinstance(x, long):
            mtnum = cpp_rxname.mt(<extra_types.uint32> x)
    else:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            from_nuc = cpp_nucname.id(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            from_nuc = cpp_nucname.id(<int> x)
        if isinstance(y, basestring):
            y_bytes = y.encode()
            to_nuc = cpp_nucname.id(std_string(<char *> y_bytes))
        elif isinstance(y, int):
            to_nuc = cpp_nucname.id(<int> y)
        z_bytes = z.encode()
        mtnum = cpp_rxname.mt(from_nuc, to_nuc, std_string(<char *> z_bytes))
    return int(mtnum)


def label(x, y=None, z="n"):
    """label(x, y=None, z="n")

    Gives a short reaction label, useful for user interfaces.

    Parameters
    ----------
    x : str, int, or long
        name, abbreviation, id, MT number, or from nuclide.
    y : str, int, or None, optional
        to nuclide.
    z : str, optional
        incident particle type ("n", "p", ...) when x and y are nuclides.

    Returns
    -------
    lab : str
        a reaction label.
    """
    cdef std_string clab
    cdef int from_nuc, to_nuc
    if y is None:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            clab = cpp_rxname.label(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            clab = cpp_rxname.label(<extra_types.uint32> long(x))
        elif isinstance(x, long):
            clab = cpp_rxname.label(<extra_types.uint32> x)
    else:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            from_nuc = cpp_nucname.id(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            from_nuc = cpp_nucname.id(<int> x)
        if isinstance(y, basestring):
            y_bytes = y.encode()
            to_nuc = cpp_nucname.id(std_string(<char *> y_bytes))
        elif isinstance(y, int):
            to_nuc = cpp_nucname.id(<int> y)
        z_bytes = z.encode()
        clab = cpp_rxname.label(from_nuc, to_nuc, std_string(<char *> z_bytes))
    lab = bytes(<char *> clab.c_str()).decode()
    return lab


def doc(x, y=None, z="n"):
    """doc(x, y=None, z="n")

    Gives documentation string for the reaction.

    Parameters
    ----------
    x : str, int, or long
        name, abbreviation, id, MT number, or from nuclide.
    y : str, int, or None, optional
        to nuclide.
    z : str, optional
        incident particle type ("n", "p", ...) when x and y are nuclides.

    Returns
    -------
    d : str
        a reaction docstring.
    """
    cdef std_string cd
    cdef int from_nuc, to_nuc
    if y is None:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            cd = cpp_rxname.doc(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            cd = cpp_rxname.doc(<extra_types.uint32> long(x))
        elif isinstance(x, long):
            cd = cpp_rxname.doc(<extra_types.uint32> x)
    else:
        if isinstance(x, basestring):
            x_bytes = x.encode()
            from_nuc = cpp_nucname.id(std_string(<char *> x_bytes))
        elif isinstance(x, int):
            from_nuc = cpp_nucname.id(<int> x)
        if isinstance(y, basestring):
            y_bytes = y.encode()
            to_nuc = cpp_nucname.id(std_string(<char *> y_bytes))
        elif isinstance(y, int):
            to_nuc = cpp_nucname.id(<int> y)
        z_bytes = z.encode()
        cd = cpp_rxname.doc(from_nuc, to_nuc, std_string(<char *> z_bytes))
    d = bytes(<char *> cd.c_str()).decode()
    return d


def child(nuc, rx, z="n"):
    """child(nuc, rx, z="n")

    Gives the child nuclide that comes from a parent and a reaction.

    Parameters
    ----------
    nuc : str or int
        parent nuclide name or id.
    rx : str or int
        reaction name or id.
    z : str, optional
        incident particle type ("n", "p", ...).

    Returns
    -------
    to_nuc : int
        a nuclide identifier.
    """
    cdef std_string ptype
    cdef int to_nuc
    cdef bint nuc_is_str = isinstance(nuc, basestring)
    cdef bint z_is_str = isinstance(z, basestring)
    cdef bint rx_is_str = isinstance(rx, basestring)
    # convert particle to std::string
    z_bytes = z.encode() if z_is_str else z
    ptype = std_string(<char *> z_bytes)
    # call child function
    if nuc_is_str and rx_is_str:
        nuc_bytes = nuc.encode()
        rx_bytes = rx.encode()
        to_nuc = cpp_rxname.child(std_string(<char *> nuc_bytes),
                                  std_string(<char *> rx_bytes), ptype)
    elif not nuc_is_str and rx_is_str:
        rx_bytes = rx.encode()
        to_nuc = cpp_rxname.child(<int> nuc, std_string(<char *> rx_bytes), ptype)
    elif nuc_is_str and not rx_is_str:
        nuc_bytes = nuc.encode()
        to_nuc = cpp_rxname.child(std_string(<char *> nuc_bytes),
                                  <extra_types.uint32> long(rx), ptype)
    elif not nuc_is_str and not rx_is_str:
        to_nuc = cpp_rxname.child(<int> nuc, <extra_types.uint32> long(rx), ptype)
    return int(to_nuc)


def parent(nuc, rx, z="n"):
    """parent(nuc, rx, z="n")

    Gives the parent nuclide that produces a child from a reaction.

    Parameters
    ----------
    nuc : str or int
        child nuclide name or id.
    rx : str or int
        reaction name or id.
    z : str, optional
        incident particle type ("n", "p", ...).

    Returns
    -------
    from_nuc : int
        a nuclide identifier.
    """
    cdef std_string ptype
    cdef int from_nuc
    cdef bint nuc_is_str = isinstance(nuc, basestring)
    cdef bint z_is_str = isinstance(z, basestring)
    cdef bint rx_is_str = isinstance(rx, basestring)
    # convert particle to std::string
    z_bytes = z.encode() if z_is_str else z
    ptype = std_string(<char *> z_bytes)
    # call parent function
    if nuc_is_str and rx_is_str:
        nuc_bytes = nuc.encode()
        rx_bytes = rx.encode()
        from_nuc = cpp_rxname.parent(std_string(<char *> nuc_bytes),
                                    std_string(<char *> rx_bytes), ptype)
    elif not nuc_is_str and rx_is_str:
        rx_bytes = rx.encode()
        from_nuc = cpp_rxname.parent(<int> nuc, std_string(<char *> rx_bytes), ptype)
    elif nuc_is_str and not rx_is_str:
        nuc_bytes = nuc.encode()
        from_nuc = cpp_rxname.parent(std_string(<char *> nuc_bytes),
                                    <extra_types.uint32> long(rx), ptype)
    elif not nuc_is_str and not rx_is_str:
        from_nuc = cpp_rxname.parent(<int> nuc, <extra_types.uint32> long(rx), ptype)
    return int(from_nuc)

