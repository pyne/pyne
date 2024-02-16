"""Python wrapper for material library."""
from __future__ import division, unicode_literals

# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.set cimport set as cpp_set
#from cython.operator cimport reference as ref
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport malloc, free
from libcpp.string cimport string as std_string
from libcpp.map cimport map as cpp_map
from libcpp.vector cimport vector as cpp_vector
from libcpp cimport bool as cpp_bool

# Python imports
import collections
try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections
cimport numpy as np
import numpy as np
from pyne.utils import QA_warn
import os
import sys
if sys.version_info[0] >= 3:
    #Python2 basestring is now Python3 string
    basestring = str

import tables as tb

# local imports
from pyne cimport cpp_material
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv

cimport cpp_jsoncpp
cimport pyne.jsoncpp as jsoncpp
import pyne.jsoncpp as jsoncpp

cimport pyne.nucname as nucname
import pyne.nucname as nucname

cimport pyne.data as data
import pyne.data as data


QA_warn(__name__)

# Maximum 32-bit signed int
DEF INT_MAX = 2147483647

cdef cpp_map[int, double] dict_to_comp(dict nucvec):
    """Converts a dictionary with arbitraily-typed keys to component map."""
    cdef int key_zz
    cdef cpp_map[int, double] comp = cpp_map[int, double]()

    for key, value in nucvec.items():
        key_zz = nucname.id(key)
        comp[key_zz] = value

    return comp


cdef class _Material:

    def __cinit__(self, nucvec=None, double mass=-1.0, double density=-1.0,
                  double atoms_per_molecule=-1.0, metadata=None, bint free_mat=True,
                  *args, **kwargs):
        """Material C++ constructor."""
        cdef cpp_map[int, double] comp
        cdef jsoncpp.Value cmetadata = jsoncpp.Value({} if metadata is None else metadata)

        if isinstance(nucvec, _Material):
            # Material from Material
            self.mat_pointer = (<_Material> nucvec).mat_pointer
        elif isinstance(nucvec, dict):
            # Material from dict
            comp = dict_to_comp(nucvec)
            self.mat_pointer = new cpp_material.Material(
                    comp, mass, density, atoms_per_molecule, deref(cmetadata._inst))
        elif isinstance(nucvec, basestring):
            # Material from file
            nucvec = nucvec.encode()
            self.mat_pointer = new cpp_material.Material(
                    <char *> nucvec, mass, density, atoms_per_molecule,
                    deref(cmetadata._inst))
        elif (nucvec is None):
            if free_mat:
                # Make empty mass stream
                self.mat_pointer = new cpp_material.Material(comp, mass, density,
                                        atoms_per_molecule, deref(cmetadata._inst))
            else:
                self.mat_pointer = NULL
        else:
            # Bad Material
            raise TypeError("The mass stream nucvec must be a dict, str, "
                    "or None, but is a {0}".format(type(nucvec)))

        # Init some meta-data
        self._comp = None
        self._free_mat = free_mat

    def __dealloc__(self):
        """Material C++ destructor."""
        if self._free_mat:
            del self.mat_pointer


    #
    # Class Attributes
    #

    property comp:
        def __get__(self):
            cdef conv._MapIntDouble comp_proxy
            if self._comp is None:
                comp_proxy = conv.MapIntDouble(False, False)
                comp_proxy.map_ptr = &self.mat_pointer.comp
                self._comp = comp_proxy
            return self._comp

        def __set__(self, value):
            cdef cpp_pair[int, double] item
            cdef cpp_map[int, double]  m

            if isinstance(value, conv._MapIntDouble):
                self.mat_pointer.comp = deref(
                        (<conv._MapIntDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[int, double]()
                for k, v in value.items():
                    item = cpp_pair[int, double](k, v)
                    m.insert(item)
                self.mat_pointer.comp = m
            elif hasattr(value, '__len__'):
                m = cpp_map[int, double]()
                for i in value:
                    item = cpp_pair[int, double](i[0], i[1])
                    m.insert(item)
                self.mat_pointer.comp = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(
                        type(value)))

            self._comp = None

    property mass:
        def __get__(self):
            return self.mat_pointer.mass

        def __set__(self, double value):
            self.mat_pointer.mass = value

    property density:
        def __get__(self):
            return self.mat_pointer.density

        def __set__(self, double value):
            self.mat_pointer.density = value

    property atoms_per_molecule:
        def __get__(self):
            return self.mat_pointer.atoms_per_molecule

        def __set__(self, double value):
            self.mat_pointer.atoms_per_molecule = value

    property metadata:
        def __get__(self):
            cdef jsoncpp.Value val = jsoncpp.Value(view=True)
            val._inst = &self.mat_pointer.metadata
            return val

        def __set__(self, value):
            cdef jsoncpp.Value val = jsoncpp.Value(value)
            val._view = True
            self.mat_pointer.metadata = deref(val._inst)

    #
    # Class Methods
    #

    def norm_comp(self):
        """Normalizes the composition, preserving the mass of the nuclide
        vector as mass.

        """
        self.mat_pointer.norm_comp()


    def from_hdf5(self, filename, datapath, int row=-1,
                  int protocol=1):
        """from_hdf5(char * filename, char * datapath, int row=-1, int protocol=1)
        Initialize a Material object from an HDF5 file.

        Parameters
        ----------
        filename : str
            Path to HDF5 file that contains the data to read in.
        datapath : str
            Path to HDF5 table or group that represents the data.
            In the example below, datapath = "/mat_name".
        row : int, optional
            The index of the arrays from which to read the data.  This
            ranges from 0 to N-1.  Defaults to the last element of the array.
            Negative indexing is allowed (row[-N] = row[0]).
        protocol : int, optional
            Specifies the protocol to use to read in the data.  Different
            protocols are used to represent different internal structures in
            the HDF5 file.

        Notes
        -----
        There are currently two protocols which are implemented for how to
        store materials inside of an HDF5 file.  Protocol 0 is the older,
        deprecated method using a group of arrays.  Protocol 1 is the newer,
        prefered method which uses a table of materials plus a side array of
        nuclides.

        The Protocol 0 HDF5 representation of a Material is a group that holds
        several extendable array datasets.  One array is entitled "Mass" while
        the other datasets are nuclide names in name form ("U235", "NP237",
        *etc*).  For example::

            file.h5 (file)
                |-- material (group)
                    |-- Mass (array)
                    |-- H1 (array)
                    |-- O16 (array)
                    |-- U235 (array)
                    |-- PU239 (array)
                    |-- ...

        The arrays are all of length N, where each row typically represents a
        different fuel cycle pass.  The sum of all of the nuclide arrays should
        sum to one, like Material.comp. This method is deprecated.

        Protocol 1 is the newer, more efficient protocol for storing many
        materials.  It consists of a table which stores the material
        information and an array that stores the nuclides (id) which index
        the comp array::

            file.h5 (file)
                |-- material (table)
                    |-- mass (double col)
                    |-- density (double col)
                    |-- atoms_per_molecule (double col)
                    |-- comp (double array col, len of nuc_zz)
                |-- nuc_zz (int array)
                |-- material_attr (variable length char array)

        The material table has a string attribute called 'nucpath' which holds
        the path to the nuclide array inside this HDF5 file.  The same nucpath
        may be used for multiple material tables.  The length of the nucpath
        must match the length of the comp arrays.

        Examples
        --------
        This method loads data into a pre-existing :class:`Material`.
        Initialization is therefore a two-step process::

            mat = Material()
            mat.from_hdf5("afile.h5", "/foo/bar/mat", -3)

        """
        cdef char * c_filename
        if isinstance(filename, unicode):
            filename_bytes = filename.encode('UTF-8')
        else:
            filename_bytes = filename
        c_filename = filename_bytes
        cdef char * c_datapath
        if isinstance(datapath, unicode):
            datapath_bytes = datapath.encode('UTF-8')
        else:
            datapath_bytes = datapath
        c_datapath = datapath_bytes
        self.mat_pointer.from_hdf5(c_filename, c_datapath, row, protocol)


    def write_hdf5(self, filename, datapath="/mat_name", nucpath="",
                   row=-0.0, chunksize=100):
        """write_hdf5(filename, datapath="/mat_name", nucpath="", row=-0.0, chunksize=100)
        Writes the material to an HDF5 file, using Protocol 1 (see the
        from_hdf5() method).

        Parameters
        ----------
        filename : str
            Path to HDF5 file to write the data out to.  If the file does not
            exist, it will be created.
        datapath : str, optional
            Path to HDF5 table that represents the data.  If the table does not
            exist, it will be created.
        nucpath : str, optional
            Path to id array of nuclides to write out.  If this array does
            not exist, it is created with the nuclides present in this
            material. Nuclides present in this material but not in nucpath will
            not be written out.
        row : float, optional
            The row index of the HDF5 table to write this material to.  This
            ranges from 0 to N.  Negative indexing is allowed (row[-N] =
            row[0]).  Defaults to the appending this material to the table
            (row[N] = row[-0.0]).  This value must be a float since in integer
            repesentation 0 **is** -0, but in float representation 0.0 **is
            not** -0.0.
        chunksize : int, optional
            In Protocol 1, materials are stored in an HDF5 table which is an
            extensible data type. The chunksize determines the number of rows
            per chunk.  For better performance, this number should be as close
            as possible to the final table size.  This parameter is only
            relevant if a new table is being created.

        Examples
        --------
        The following writes out ten low-enriched uranium materials to a new
        table::

            leu = Material({'U235': 0.04, 'U238': 0.96}, 4.2, "LEU", 1.0)
            leu.write_hdf5('proto1.h5', chunksize=10)

            for i in range(2, 11):
                leu = Material({'U235': 0.04, 'U238': 0.96}, i*4.2, "LEU",
                               1.0*i)
                leu.write_hdf5('proto1.h5')

        """
        cdef char * c_filename
        filename_bytes = filename.encode('UTF-8')
        c_filename = filename_bytes
        cdef char * c_datapath
        datapath_bytes = datapath.encode('UTF-8')
        c_datapath = datapath_bytes
        cdef char * c_nucpath

        if nucpath != "":
            nucpath_bytes = nucpath.encode('UTF-8')
            c_nucpath = nucpath_bytes
            self.mat_pointer.deprecated_write_hdf5(c_filename, c_datapath, c_nucpath, row, chunksize)
        else:
            self.mat_pointer.write_hdf5(c_filename, c_datapath, row, chunksize)


    def phits(self, frac_type='mass', mult_den=True):
        """phits(frac_type='mass', mult_den=True)
        Return an phits card
        Parameters
        ----------
        frac_type : str, optional
            Either 'mass' or 'atom'. Speficies whether mass or atom fractions
            are used to describe material composition. (default 'mass')
        mult_den : bool, optional
            Flag for whether material cards are written in mass density if True,
            or mass fraction if False. (default True)
        """
        cdef std_string card
        card = self.mat_pointer.phits(frac_type.encode(), mult_den)
        return card.decode()

    def mcnp(self, frac_type='mass', mult_den=True):
        """mcnp(frac_type='mass', mult_den=True)
        Return an mcnp card
        Parameters
        ----------
        frac_type : str, optional
            Either 'mass' or 'atom'. Speficies whether mass or atom fractions
            are used to describe material composition. (default 'mass')
        mult_den : bool, optional
            Flag for whether material cards are written in mass density if True,
            or mass fraction if False. (default True)
        """
        cdef std_string card
        card = self.mat_pointer.mcnp(frac_type.encode(), mult_den)
        return card.decode()

    def get_uwuw_name(self):
        """get_uwuw_name()
        Return a uwuw material name
        """
        cdef std_string uwuw_name
        uwuw_name = self.mat_pointer.get_uwuw_name()
        return uwuw_name.decode()

    def openmc(self, frac_type='mass', indent_lvl=1):
        """openmc(frac_type)
        Return an openmc xml element for the material
        """
        cdef std_string mat
        mat = self.mat_pointer.openmc(frac_type.encode(), indent_lvl)
        return mat.decode()

    def fluka(self, fid, frac_type='mass'):
        """fluka()
        Return a fluka material record if there is only one component,
        otherwise return the compound material record and the fluka
        compound record
        Parameters
        ----------
        The sequential material id starting from 26 unless predefined
        """
        cdef std_string card
        card = self.mat_pointer.fluka(fid, frac_type.encode())
        return card.decode()

    def not_fluka_builtin(self, fluka_name):
        """not_fluka_builtin()
        Return whether a string is in the fluka built-in list
        Parameters
        ----------
        A string representing a FLUKA material name
        """
        cdef cpp_bool card
        card = self.mat_pointer.not_fluka_builtin(fluka_name)
        return card

    def fluka_material_str(self, id):
        """fluka_material_str()
        Return the FLUKA MATERIAL record with the given id.
        A single-component material is expected
        Parameters
        ----------
        The sequential material id starting from 26 unless predefined
        """
        cdef std_string card
        card = self.mat_pointer.fluka_material_str(id)
        return card

    def fluka_material_component(self, id, nucid, fluka_name):
        """fluka_material_component()
        Return the FLUKA MATERIAL record with the given id, nucid and name.
        Parameters
        ----------
        The sequential material id, the (single) nucid, and the fluka name
        """
        cdef std_string card
        card = self.mat_pointer.fluka_material_component(id, nucid, fluka_name)
        return card

    def fluka_material_line(self, znum, mass, id, name):
        """fluka_material_line()
        Return the FLUKA MATERIAL record with the given znum, atomic mass, id,
        and fluka name
        Parameters
        ----------
        The znum, atomic mass, material id, and the fluka name
        """
        cdef std_string card
        card = self.mat_pointer.fluka_material_line(znum, mass, id, name)
        return card

    def fluka_format_field(self, field):
        """fluka_format_field()
        Return a string for a single field in the FLUKA MATERIAL record
        Parameters
        ----------
        The field value
        """
        cdef std_string card
        card = self.mat_pointer.fluka_format_field(field)
        return card

    def fluka_compound_str(self, id, frac_type='mass'):
        """fluka_compound_str()
        Return the FLUKA MATERIAL record for the compound, and the
    FLUKA COMPOUND record for the components
        Parameters
        ----------
        The sequential compound id starting from 26 unless predefined
        """
        cdef std_string card
        card = self.mat_pointer.fluka_compound_str(id, frac_type)
        return card

    def from_text(self, filename):
        """from_text(char * filename)
        Initialize a Material object from a simple text file.

        Parameters
        ----------
        filename : str
            Path to text file that contains the data to read in.

        Notes
        -----
        The text representation of Materials are nuclide identifiers in the
        first column and mass or weight values in the second column.  For
        example, for natural uranium::

            922340  0.000055
            U235    0.00720
            92238   0.992745

        Data in this file must be whitespace separated.  Any valid nuclide
        naming scheme may be used for the nuclide identifiers.  Moreover,
        material metadata may be optionally supplied::

            Name    NatU
            Mass    42.0
            APerM   1
            922340  0.000055
            U235    0.00720
            92238   0.992745

        Examples
        --------
        This method loads data into a pre-existing Material.
        Initialization is therefore a two-step process::

            mat = Material()
            mat.from_text("natu.txt")

        This method is most often called implicitly by the Material constructor.

        """
        cdef char * c_filename
        if isinstance(filename, unicode):
            filename_bytes = filename.encode('UTF-8')
        else:
            filename_bytes = filename
        c_filename = filename_bytes
        self.mat_pointer.from_text(c_filename)


    def write_text(self, filename):
        """write_text(filename)
        Writes the material to a plain text file.

        Parameters
        ----------
        filename : str
            Path to text file to write the data to.  If the file already
            exists, it will be overwritten.

        Examples
        --------
        The following writes out a low-enriched uranium material to a new file::

            leu = Material({'U235': 0.04, 'U238': 0.96}, 42.0, "LEU", 1.0)
            leu.write_text('leu.txt')

        """
        cdef char * c_filename
        if isinstance(filename, unicode):
            filename_bytes = filename.encode('UTF-8')
        else:
            filename_bytes = filename
        c_filename = filename_bytes
        self.mat_pointer.write_text(c_filename)

    def load_json(self, json):
        """load_json(json)
        Loads a JSON instance into this Material.

        Parameters
        ----------
        json : jsoncpp.Value
            An object-type JSON value.

        """
        self.mat_pointer.load_json(deref((<jsoncpp.Value> json)._inst))

    def dump_json(self):
        """dump_json()
        Dumps the material to a JSON object.

        Returns
        -------
        val : jsoncpp.Value
            An object-type JSON value.

        """
        cdef jsoncpp.Value val = jsoncpp.Value(view=False)
        val._inst[0] = self.mat_pointer.dump_json()
        return val

    def from_json(self, filename):
        """from_json(char * filename)
        Initialize a Material object from a JSON file.

        Parameters
        ----------
        filename : str
            Path to text file that contains the data to read in.

        """
        cdef char * c_filename
        filename_bytes = filename.encode('UTF-8')
        c_filename = filename_bytes
        self.mat_pointer.from_json(c_filename)

    def write_json(self, filename):
        """write_json(filename)
        Writes the material to a JSON file.

        Parameters
        ----------
        filename : str
            Path to text file to write the data to.  If the file already
            exists, it will be overwritten.

        Examples
        --------
        The following writes out a low-enriched uranium material to a new file::

            leu = Material({'U235': 0.04, 'U238': 0.96}, 42.0, "LEU", 1.0)
            leu.write_json('leu.json')

        """
        filename = filename.encode()
        self.mat_pointer.write_json(filename)

    def normalize(self):
        """This convenience method normalizes the mass stream by setting its
        mass = 1.0.

        """
        self.mat_pointer.normalize()


    def mult_by_mass(self):
        """This multiplies comp by mass and returns the resultant
        nuctopic vector.

        Returns
        -------
        nucvec : dict
            For a Material mat,

            .. math:: \\mbox{nucvec[nuc]} = \\mbox{mat.comp[nuc]} \\times \\mbox{mat.mass}

        """
        cdef conv._MapIntDouble nucvec_proxy = conv.MapIntDouble()
        nucvec_proxy.map_ptr = new cpp_map[int, double](
                self.mat_pointer.mult_by_mass())
        return nucvec_proxy


    def activity(self):
        """This provides the activity of the comp of the material. It assumes
        that the mass of the material is given in units of [grams] and returns
        activities in units of [Bq].

        Returns
    -------
    nucvec : dict
        For a Material mat

        """
        cdef conv._MapIntDouble nucvec_proxy = conv.MapIntDouble()
        nucvec_proxy.map_ptr = new cpp_map[int, double](
                self.mat_pointer.activity())
        return nucvec_proxy


    def decay_heat(self):
        """This provides the decay heat using the comp of the the Material. It
        assumes that the composition of material is given in units of [grams]
        and returns decay heat in units of [MW].

        Returns
        -------
        nucvec : dict
            For a Material mat
        """
        cdef conv._MapIntDouble nucvec_proxy = conv.MapIntDouble()
        nucvec_proxy.map_ptr = new cpp_map[int, double](
                self.mat_pointer.decay_heat())
        return nucvec_proxy


    def dose_per_g(self, dose_type, source=0):
        """This provides the dose per gram using the comp of the the Material.

        Parameters
        ----------
        dose_type : string
            One of: ext_air, ext_soil, ingest, inhale
        source : int
            optional; default is EPA
            0 for EPA, 1 for DOE, 2 for GENII

        Returns
        -------
        nucvec : dict
            For a Material mat:
            ext_air_dose returns mrem/h per g per m^3
            ext_soil_dose returns mrem/h per g per m^2
            ingest_dose returns mrem per g
            inhale_dose returns mrem per g
        """
        cdef conv._MapIntDouble nucvec_proxy = conv.MapIntDouble()
        cdef std_string dosetype
        if not isinstance(dose_type, bytes):
            dose_type = dose_type.encode()
        dosetype = std_string(<char *> dose_type)
        nucvec_proxy.map_ptr = new cpp_map[int, double](
                self.mat_pointer.dose_per_g(dosetype, source))
        return nucvec_proxy


    def molecular_mass(self, atoms_per_molecule=-1.0):
        """molecular_mass(atoms_per_molecule=-1.0)
        This method returns the molecular mass of the comp of this
        material.

        Parameters
        ----------
        atoms_per_molecule : double, optional
            Number of atoms to per molecule of material.  Needed to obtain
            proper scaling.  For example, this value for water is 3.0.

        Returns
        -------
        mol_mass : float
            Molecular mass in [amu].

        """
        return self.mat_pointer.molecular_mass(atoms_per_molecule)

    def expand_elements(self, nucset=set()):
        """expand_elements(self)
        Expands the elements ('U', 'C', etc) in the material by
        replacing them with their natural isotopic distributions with
        the exception of the ids in nucset. This function returns a
        copy.

        Parameters
        ----------
        nucset : set, optional
            A set of integers representing nucids which should not
            be expanded.

        Returns
        -------
        newmat : Material
            A copied and expanded material.

        """
        cdef _Material newmat = Material()
        newmat.mat_pointer[0] = self.mat_pointer.expand_elements(nucset)
        return newmat

    def collapse_elements(self, nucset):
        """collapse_elements(self, nucset)
        Collapses the elements in the material, excluding the nucids in
        the set nucset. This function returns a copy of the material.

        Parameters
        ----------
        nucset : set, optional
            A set of integers representing nucids which should not
            be collapsed.

        Returns
        -------
        newmat : Material
            A copied and collapsed material.

        """
        cdef _Material newmat = Material()
        newmat.mat_pointer[0] = self.mat_pointer.collapse_elements(nucset)
        return newmat

    def mass_density(self, double num_dens=-1.0, double atoms_per_molecule=-1.0):
        """mass_density(self, num_dens=-1.0, atoms_per_molecule=-1.0)
        Computes, sets, and returns the mass density when num_dens is greater
        than or equal zero.  If num_dens is negative, this simply returns the
        current value of the density attribute.

        Parameters
        ----------
        num_dens : float, optional
            The number density from which to compute the mass density in units
            of [1/cc].
        atoms_per_molecule : float, optional
            Number of atoms to per molecule of material. For example, this value
            for water is 3.0.

        Returns
        -------
        density : float
            The density attr [g/cc].

        """
        return self.mat_pointer.mass_density(num_dens, atoms_per_molecule)

    def number_density(self, double mass_dens=-1.0, double atoms_per_molecule=-1.0):
        """number_density(self, mass_dens=-1.0, atoms_per_molecule=-1.0)
        Computes and returns the number density from the mass_dens argument if this
        is greater than or equal zero.  If mass_dens is negative, then the number
        density is computed using the current value of the density attribute.

        Parameters
        ----------
        mass_dens : float, optional
            The mass density from which to compute the number density in units
            of [g/cc].
        atoms_per_molecule : float, optional
            Number of atoms to per molecule of material. For example, this value
            for water is 3.0.

        Returns
        -------
        num_dens : float
            The number density [1/cc] of the material.

        """
        return self.mat_pointer.number_density(mass_dens, atoms_per_molecule)



    #
    # submaterial Methods
    #

    def sub_mat(self, nuc_sequence):
        """sub_mat(nuc_sequence)
        Grabs a subset of the material and returns a new material comprised
        of only the specified nuclides.

        Parameters
        ----------
        nuc_sequence : sequence
            Nuctopes --OR-- elements to be taken from current stream.
            Members of this list must be integers.  For example, [922350, 942390]
            would take U-235 and Pu-239.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only
            has the members given in nuc_sequence.  The mass of the submaterial
            is calculated based on the mass fraction composition and mass
            of the original mass stream.

        Notes
        -----
        The input here is seen as a suggestion and so no error is raised if a
        nuclide is asked for via nuc_sequence that is not present in the
        original material.

        """
        # Make an nuctopic set
        cdef cpp_set[int] nuc_set = nucname.id_set(nuc_sequence)

        # Make new python version of this material
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_mat(nuc_set)
        return pymat


    def set_mat(self, nuc_sequence, value):
        """set_mat(nuc_sequence, value)
        Sets a subset of the material to a new value and returns a new
        material.

        Parameters
        ----------
        nuc_sequence : sequence
            Nuctopes --OR-- elements to be taken from current stream.
            Members of this list must be integers.  For example, [922350, 942390]
            would take U-235 and Pu-239.
        value : float
            Mass value to set all nuclides in sequence to on the material.

        Returns
        -------
        submaterial : Material
            A new material object whose members in nuc_sequence have the
            cooresponding mass value.  The mass of the submaterial is
            calculated based on the mass fraction composition and mass of the
            original material.

        """
        # Make an nuctopic set
        cdef cpp_set[int] nuc_set = nucname.id_set(nuc_sequence)

        # Make new python version of this material
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.set_mat(nuc_set, <double> value)
        return pymat


    def del_mat(self, nuc_sequence):
        """del_mat(nuc_sequence)
        Removes a subset of the material and returns a new material
        comprised of only the non-specified nuclides.

        Parameters
        ----------
        nuc_sequence : sequence
            Nuclides to be taken out of the current material.

        Returns
        -------
        submaterial : Material
            A new material object that only has the members not given in
            nuc_sequence.  The mass of the submaterial is calculated based on
            the mass fraction composition and mass of the original material.

        Notes
        -----
        The input here is seen as a suggestion and so no error is raised if a
        nuclide is asked for via nuc_sequence that is not present in the
        original material.

        """
        # Make an nuctopic set
        cdef cpp_set[int] nuc_set = nucname.id_set(nuc_sequence)

        # Make new python version of this material
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.del_mat(nuc_set)
        return pymat


    def sub_range(self, lower=0, upper=INT_MAX):
        """sub_range(lower=0, upper=INT_MAX)
        Grabs a sub-material from this mat based on a range [lower, upper)
        of values.

        Parameters
        ----------
        lower : nuclide-name, optional
            Lower bound on nuclide range.
        upper : nuclide-name, optional
            Upper bound on nuclide range.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has nuclides on the given range.

        """
        cdef int clower, cupper

        if isinstance(lower, int):
            clower = lower
        else:
            clower = nucname.id(lower)

        if isinstance(upper, int):
            cupper = upper
        else:
            cupper = nucname.id(upper)

        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_range(clower, cupper)
        return pymat


    def set_range(self, lower=0, upper=INT_MAX, value=0.0):
        """set_range(lower=0, upper=INT_MAX, value=0.0)
        Sets a sub-material from this mat based on a range [lower, upper) to
        a new mass weight value.

        Parameters
        ----------
        lower : nuclide-name, optional
            Lower bound on nuclide range.
        upper : nuclide-name, optional
            Upper bound on nuclide range.
        value : float
            Mass value to set all nuclides on the range to on the material.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has nuclides on the given range.

        """
        cdef int clower, cupper

        if isinstance(lower, int):
            clower = lower
        else:
            clower = nucname.id(lower)

        if isinstance(upper, int):
            cupper = upper
        else:
            cupper = nucname.id(upper)

        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.set_range(clower, cupper, <double> value)
        return pymat


    def del_range(self, lower=0, upper=INT_MAX):
        """del_range(lower=0, upper=INT_MAX)
        Remove a range [lower, upper) of nuclides from this material and
        returns a submaterial.

        Parameters
        ----------
        lower : nuclide-name, optional
            Lower bound on nuclide range.
        upper : nuclide-name, optional
            Upper bound on nuclide range.

        Returns
        -------
        submaterial : Material
            A new mass stream object that does not contain nuclides on the
            given range.

        """
        cdef int clower, cupper

        if isinstance(lower, int):
            clower = lower
        else:
            clower = nucname.id(lower)

        if isinstance(upper, int):
            cupper = upper
        else:
            cupper = nucname.id(upper)

        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.del_range(clower, cupper)
        return pymat


    def sub_elem(self, element):
        """sub_elem(element)
        Grabs a subset of the material and returns a new material comprised of
        only the nuclides of the specified element.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has members of the given element.

        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_elem(nucname.id(element))
        return pymat


    def sub_u(self):
        """sub_u()
        Convenience method that gets the Uranium portion of a mass stream.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has Uranium members.

        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_elem(nucname.id('U'))
        return pymat


    def sub_pu(self):
        """sub_pu()
        Convenience method that gets the Plutonium portion of a mass stream.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has Plutonium members.

        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_elem(nucname.id('Pu'))
        return pymat


    def sub_lan(self):
        """sub_lan()
        Convenience method that gets the Lanthanide portion of a mass stream.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has Lanthanide members.

        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_lan()
        return pymat


    def sub_act(self):
        """sub_act()
        Convenience method that gets the Actinide portion of a mass stream.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has Actinide members.

        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_act()
        return pymat


    def sub_tru(self):
        """sub_tru()
        Convenience method that gets the Transuranic portion of a mass
        stream.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has Transuranic members.

        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_tru()
        return pymat


    def sub_ma(self):
        """sub_ma()
        Convenience method that gets the Minor Actinide portion of a mass
        stream.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has Minor Actinide members.

        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_ma()
        return pymat


    def sub_fp(self):
        """sub_fp()
        Convenience method that gets the Fission Product portion of a mass
        stream.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has Fission Product members.

        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_fp()
        return pymat


    #
    # Atom Fraction Methods
    #
    def to_atom_frac(self):
        """Converts the material to a map of nuclides to atom fractions.

        Returns
        -------
        atom_fracs : mapping
            Dictionary-like object that maps nuclides to atom fractions in the
            material.

        """
        cdef conv._MapIntDouble comp_proxy = conv.MapIntDouble()
        comp_proxy.map_ptr = new cpp_map[int, double](self.mat_pointer.to_atom_frac())
        return comp_proxy


    def from_atom_frac(self, atom_fracs):
        """from_atom_frac(atom_fracs)
        Loads the material composition based on a mapping of atom fractions.

        Parameters
        ----------
        atom_fracs : dict
            Dictionary that maps nuclides or materials to atom fractions for
            the material.  The keys may be intergers, strings, or materials.
            The values must be castable to floats.

        Examples
        --------
        To get a material from water, based on atom fractions::

            h2o = {10010: 2.0, 'O16': 1.0}
            mat = Material()
            mat.from_atom_frac(h2o)

        Or for Uranium-Oxide, based on an initial fuel vector::

            # Define initial heavy metal
            ihm = Material()
            ihm.from_atom_frac({'U235': 0.05, 'U238': 0.95})

            # Define Uranium-Oxide
            uox = {ihm: 1.0, 80160: 2.0}
            mat = Material()
            mat.from_atom_frac(uox)

        Note that the initial heavy metal was used as a key in a dictionary.
        This is possible because Materials are hashable.

        """
        cdef int key_zz
        cdef double val
        cdef cpp_map[int, double] key_af
        cdef cpp_map[int, double].iterator keyiter, keyend
        cdef cpp_map[int, double] af = cpp_map[int, double]()

        # Convert atom_fracs to something usable in C++
        for key, value in atom_fracs.items():
            val = <double> value
            if isinstance(key, int):
                key_zz = <int> nucname.id(key)
                if 0 == af.count(key_zz):
                    af[key_zz] = 0.0
                af[key_zz] = af[key_zz] + val
            elif isinstance(key, basestring):
                key_zz = nucname.id(key)
                if 0 == af.count(key_zz):
                    af[key_zz] = 0.0
                af[key_zz] = af[key_zz] + val
            elif isinstance(key, _Material):
                key_af = deref((<_Material> key).mat_pointer).to_atom_frac()
                keyiter = key_af.begin()
                keyend = key_af.end()
                while keyiter != keyend:
                    key_zz = deref(keyiter).first
                    if 0 == af.count(key_zz):
                        af[key_zz] = 0.0
                    af[key_zz] = af[key_zz] + (val * deref(keyiter).second)
                    inc(keyiter)
            else:
                raise TypeError("Atom fraction keys must be integers, "
                        "strings, or Materials.")

        self.mat_pointer.from_atom_frac(af)

    def from_activity(self, activities):
        """from_activity(activities)
        Loads the material composition based on a mapping of radionuclide
        activities. It assumes that activities are supplied in units of [Bq]
        and sets the material mass to units of [grams].

        Parameters
        ----------
        activities : dict
            Dictionary that maps radionuclides to activities for the material.
            The keys may be intergers or strings. The values must be castable
            to floats.

        Examples
        --------
        To get a material of natural uranium, based on activities::

            natu = {'U234': 12223.2, 'U235': 568.648, 'U238': 12347.1}
            mat = Material()
            mat.from_activity(natu)

        """
        cdef int key_zz
        cdef double val
        cdef cpp_map[int, double] key_act
        cdef cpp_map[int, double].iterator keyiter, keyend
        cdef cpp_map[int, double] act = cpp_map[int, double]()

        # Convert atom_fracs to something usable in C++
        for key, value in activities.items():
            val = <double> value
            if isinstance(key, int):
                key_zz = <int> nucname.id(key)
            elif isinstance(key, basestring):
                key_zz = nucname.id(key)
            else:
                raise TypeError("Activity keys must be integers, "
                        "or strings.")

            if 0 == act.count(key_zz):
                act[key_zz] = 0.0
            act[key_zz] = act[key_zz] + val

        self.mat_pointer.from_activity(act)

    def to_atom_dens(self):
        """Converts the material to a map of nuclides to atom densities.

        Returns
        -------
        atom_dens : mapping
            Dictionary-like object that maps nuclides to atom densites in the
            material.

        """
        cdef conv._MapIntDouble comp_proxy = conv.MapIntDouble()
        comp_proxy.map_ptr = new cpp_map[int, double](self.mat_pointer.to_atom_dens())
        return comp_proxy


    #
    # Radioactive Properties
    #

    def gammas(self, norm=False):
        """
        Returns a vector of gamma rays and intensities in decays/s/atom material

        Returns
        -------
        gammas : a vector of pairs of gamma-rays and intensities. The
            intensities are in decays/s/atom material
        """
        return self.mat_pointer.gammas()

    def xrays(self, norm=False):
        """
        Returns a vector of X rays and intensities in decays/s/atom material.
        Includes only X rays from internal conversion and electron capture

        Returns
        -------
        x-rays : a vector of pairs of X-rays and intensities. The
            intensities are in decays/s/atom material
        """
        return self.mat_pointer.xrays()

    def photons(self, norm=False):
        """
        Returns a vector of photons and intensities in decays/s/atom material.
        This vector is the combination of X-rays and gamma-rays produced in the
        decay of the material.


        Parameters
        ----------
        norm : boolean
            Whether or not to normalize the returned data if True then
            intensities

        Returns
        -------
        photons : a vector of pairs of photon energies and intensities. The
            intensities are in decays/s/atom material
        """
        return self.mat_pointer.photons(<cpp_bool> norm)

    def decay(self, double t):
        """decay(double t)
        Decays a material for a time t, in seconds. Returns a new material.
        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.decay(t)
        return pymat

    def cram(self, A, int order=14):
        """Transmutes the material via the CRAM method.

        Parameters
        ----------
        A : 1D array-like
            The transmutation matrix [unitless]
        order : int, optional
            The CRAM approximation order (default 14).

        Returns
        -------
        A new material which has been transmuted.
        """
        A = np.asarray(A, dtype=np.float64)
        cdef int Alen = len(A)
        cdef double* Aptr = <double*> np.PyArray_DATA(A)
        cdef cpp_vector[double] cpp_A = cpp_vector[double]()
        cpp_A.reserve(Alen)
        cpp_A.assign(Aptr, Aptr + Alen)
        # transmute and return
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.cram(cpp_A, order)
        return pymat


    #
    # Operator Overloads
    #

    # Addition

    def __add_float__(x, double y):
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = x.mat_pointer[0] + y
        return pymat


    def __add_material__(x, y):
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = (<_Material> x).mat_pointer[0] + (<_Material> y).mat_pointer[0]
        return pymat


    def __add__(x, y):
        if isinstance(x, _Material) and isinstance(y, _Material):
            return x.__add_material__(y)
        elif isinstance(y, float):
            return x.__add_float__(y)
        elif isinstance(x, float):
            return y.__add_float__(x)
        elif isinstance(y, int):
            return x.__add_float__(float(y))
        elif isinstance(x, int):
            return y.__add_float__(float(x))
        else:
            return NotImplemented


    # Multiplication

    def __mul_float__(x, double y):
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = x.mat_pointer[0] * y
        return pymat


    def __mul__(x, y):
        if isinstance(y, float):
            return x.__mul_float__(y)
        elif isinstance(x, float):
            return y.__mul_float__(x)
        elif isinstance(y, int):
            return x.__mul_float__(float(y))
        elif isinstance(x, int):
            return y.__mul_float__(float(x))
        else:
            return NotImplemented


    # Division
    def __div_float__(self, double y):
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer[0] * (1 / y)
        return pymat

    def __div__(self, y):
        if isinstance(y, float):
            return self.__div_float__(y)
        elif isinstance(y, int):
            return self.__div_float__(float(y))
        else:
            return NotImplemented

    def __rdiv__(self, y):
        return self.__div__(y)

    def __truediv__(self, y):
        if isinstance(y, float):
            return self.__div_float__(y)
        elif isinstance(y, int):
            return self.__div_float__(float(y))
        else:
            return NotImplemented

    #
    # Mapping interface
    #
    def __len__(self):
        return self.mat_pointer.comp.size()


    def __contains__(self, int key):
        if 0 < self.mat_pointer.comp.count(key):
            return True
        else:
            return False


    def __getitem__(self, key):
        cdef double key_mass

        # Get single integer-key
        if isinstance(key, int):
            if 0 < self.mat_pointer.comp.count(key):
                key_mass = self.mat_pointer.comp[key] * self.mat_pointer.mass
                return key_mass
            else:
                raise KeyError("key {0} not found".format(repr(key)))

        # Get single string-key
        elif isinstance(key, basestring):
            key_zz = nucname.id(key)
            return self[key_zz]

        # Get slice-based sub-material
        elif isinstance(key, slice):
            lower = key.start
            if lower is None:
                lower = 0

            upper = key.stop
            if upper is None:
                upper = INT_MAX

            return self.sub_range(lower, upper)

        # Get sequance-based sub-material
        elif hasattr(key, '__len__'):
            return self.sub_mat(key)

        # Fail-Yurt
        else:
            raise TypeError("key {0} is of unsupported type {1}".format(
                    repr(key), type(key)))


    def __setitem__(self, key, double value):
        cdef _Material new_mat
        cdef matp new_matp
        cdef conv._MapIntDouble mbm
        cdef int key_zz
        cdef cpp_map[int, double].iterator mbmiter, mbmend

        # Set single integer-key
        if isinstance(key, int):
            mbm = self.mult_by_mass()
            mbm.map_ptr[0][key] = value
            new_matp = new cpp_material.Material(mbm.map_ptr[0], -1.0, -1.0)
            self.mat_pointer = new_matp
            self._comp = None

        # Set single string-key
        elif isinstance(key, basestring):
            key_zz = nucname.id(key)
            self[key_zz] = value

        # Set slice-based sub-material
        elif isinstance(key, slice):
            lower = key.start
            if lower is None:
                lower = 0

            upper = key.stop
            if upper is None:
                upper = INT_MAX

            # set values back on instance
            new_mat = self.set_range(lower, upper, value)
            self.mat_pointer[0] = new_mat.mat_pointer[0]
            self._comp = None

        # Set sequance-based sub-material
        elif hasattr(key, '__len__'):
            new_mat = self.set_mat(key, value)
            self.mat_pointer[0] = new_mat.mat_pointer[0]
            self._comp = None

        # Fail-Yurt
        else:
            msg = "key {0} is of unsupported type {1}".format(repr(key), type(key))
            raise TypeError(msg)


    def __delitem__(self, key):
        cdef _Material new_mat
        cdef matp new_matp
        cdef conv._MapIntDouble mbm
        cdef int key_zz
        cdef cpp_map[int, double].iterator mbmiter, mbmend

        # Remove single key
        if isinstance(key, int):
            if 0 == self.mat_pointer.comp.count(key):
                return
            mbm = self.mult_by_mass()
            mbm.map_ptr.erase(<int> key)
            new_matp = new cpp_material.Material(mbm.map_ptr[0], -1.0)
            new_matp = new cpp_material.Material(mbm.map_ptr[0], -1.0, -1.0)
            self.mat_pointer = new_matp
            self._comp = None

        # Remove single string-key
        elif isinstance(key, basestring):
            key_zz = nucname.id(key)
            del self[key_zz]

        # Remove slice-based sub-material
        elif isinstance(key, slice):
            lower = key.start
            if lower is None:
                lower = 0

            upper = key.stop
            if upper is None:
                upper = INT_MAX

            # set values back on instance
            new_mat = self.del_range(lower, upper)
            self.mat_pointer[0] = new_mat.mat_pointer[0]
            self._comp = None

        # Remove sequance-based sub-material
        elif hasattr(key, '__len__'):
            new_mat = self.del_mat(key)
            self.mat_pointer[0] = new_mat.mat_pointer[0]
            self._comp = None

        # Fail-Yurt
        else:
            raise TypeError("key {0} is of unsupported type {1}".format(
                    repr(key), type(key)))


    def __iter__(self):
        mbm = self.mult_by_mass()
        self._iter_mbm = mbm
        mbm_iter = iter(mbm)
        return mbm_iter


    #
    # Make materials hasable so that they may be used as keys in a dictionary
    #
    def __hash__(self):
        # not the most sofisticated hash...
        return id(self)


class Material(_Material, collectionsAbc.MutableMapping):
    """Material composed of nuclides.

    Parameters
    ----------
    comp : dict or str
        This is the input nuclide component dictionary.  This dictionary need
        not be normalized; Material initialization will automatically
        renormalize the stream.  Thus the comp simply is a dictionary of
        relative mass.  The keys of comp must be integers representing
        nuclides in id-form.  The values are floats for each nuclide's
        mass fraction. If a string is provided instead of a dictionary, then
        Material will read in the comp vector from a file at the string's
        location.  This either plaintext or hdf5 files. If no comp is provided,
        an empty Material object is constructed.
    mass : float, optional
        This is the mass of the new stream. If the mass provided is negative
        (default -1.0) then the mass of the new stream is calculated from the
        sum of compdict's components before normalization.  If the mass here is
        positive or zero, then this mass overrides the calculated one.
    density : float, optional
        This is the density of the material.
    atoms_per_molecule : float, optional
        Number of atoms to per molecule of material.  Needed to obtain proper
        scaling of molecular mass.  For example, this value for water is
        3.0.
    metadata : JSON-convertable Python object, optional
        Initial attributes to build the material with.  At the top-level this is
        usually a dictionary with string keys.  This container is used to store
        arbitrary metadata about the material.
    free_mat : bool, optional
        Flag for whether this wrapper 'owns' this underlying C++ pyne::Material
        object, and thus determines whether or not to deallocate it on wrapper
        destruction.

    """
    def __str__(self):
        header = ["Material:"]
        header += ["mass = {0}".format(self.mass)]
        header += ["density = {0}".format(self.density)]
        header += ["atoms per molecule = {0}".format(self.atoms_per_molecule)]
        if self.metadata.isobject():
            for key, value in self.metadata.items():
                header += ["{0} = {1}".format(key, value)]
        header += ['-' * max([len(h) for h in header])]
        header = "\n".join(header) + "\n"

        s = header + "\n".join(["{0:<7}{1}".format(
                nucname.name(key), value) for key, value in self.comp.items()])
        return s

    def __repr__(self):
        return "pyne.material.Material({0}, {1}, {2}, {3}, {4})".format(
                repr(self.comp), self.mass, self.density, self.atoms_per_molecule, repr(self.metadata))

    def __deepcopy__(self, memo):
        cdef _Material other = Material(free_mat=False)
        cdef cpp_material.Material * self_ptr = (<_Material> self).mat_pointer
        cdef cpp_material.Material * other_ptr = new cpp_material.Material()
        other_ptr.comp = self_ptr.comp
        other_ptr.mass = self_ptr.mass
        other_ptr.density = self_ptr.density
        other_ptr.atoms_per_molecule = self_ptr.atoms_per_molecule
        other_ptr.metadata = self_ptr.metadata
        other.mat_pointer = other_ptr
        other._free_mat = True
        return other


    def write_mcnp(self, filename, frac_type='mass'):
        """write_mcnp(self, filename, frac_type='mass')
        The method appends an MCNP mass fraction definition, with
        attributes to the file with the supplied filename.

        Parameters
        ----------
        filename : str
            The file to append the material definition to.
        frac_type : str, optional
            Either 'mass' or 'atom'. Speficies whether mass or atom fractions
            are used to describe material composition.
        """
        with open(filename, 'a') as f:
            f.write(self.mcnp(frac_type))

    def write_openmc(self, filename, frac_type='mass', indent_lvl=1):
        """write_openmc(self, filename, frac_type='mass')
        The method appends an OpenMC mass fraction definition, with
        attributes to the file with the supplied filename.

        Parameters
        ----------
        filename : str
            The file to append the material definition to.
        frac_type : str, optional
            Either 'mass' or 'atom'. Speficies whether mass or atom fractions
            are used to describe material composition.
        """
        with open(filename, 'a') as f:
            f.write(self.openmc(frac_type, indent_lvl))

    def alara(self):
        """alara(self)
        This method returns an ALARA material in string form, with relevant
        attributes as ALARA valid comments.

        Returns
        -------
        s : str
            The MCNP material card.
        """
        s = ''

        if 'mat_number' in self.metadata:
            s += '# mat number: {0}\n'.format(self.metadata['mat_number'])
            mat_num = self.metadata['mat_number']  # for use in mat_name
        else:
            mat_num = '<mat_num>'

        if 'source' in self.metadata:
            s += '# source: {0}\n'.format(self.metadata['source'])

        if 'comments' in self.metadata:
            comment_string= 'comments: ' + self.metadata['comments']
            # split up lines so comments are less than 80 characters
            for n in range(0, int(np.ceil(float(len(comment_string))/77))):
                s += '# {0}\n'.format(comment_string[n*77:(n + 1)*77])

        # set density. If not present, set it equal to "<rho>"
        if str(self.density) != '-1.0':
            density = self.density
        else:
            density = '<rho>'

        # if a name is present, use it. Otherwise the name is is:
        # mat<mat_number>_rho<rho>
        if 'name' in self.metadata:
            mat_name = self.metadata['name']
        else:
            mat_name = 'mat{0}_rho-{1}'.format(mat_num, density)

        s += '{0} {1} {2}\n'\
                    .format(mat_name, density, len(self.comp))

        # Multiply frac by 100 because ALARA uses mass percent for matlib files
        for iso, frac in self.comp.items():
            s += '     {0} {1:.4E} {2}\n'.format(nucname.alara(iso),
                                                 frac*100, str(nucname.znum(iso)))

        return s

    def write_alara(self, filename):
        """write_alara(self, filename)
        The method appends an ALARA material d$efinition, with attributes
        to the file with the supplied filename.

        Parameters
        ----------
        filename : str
            The file to append the material definition to.
        """
        with open(filename, 'a') as f:
            f.write(self.alara())



#####################################
### Material generation functions ###
#####################################

def from_atom_frac(atom_fracs, double mass=-1.0, double density=-1.0,
                   double atoms_per_molecule=-1.0, metadata=None):
    """from_atom_frac(atom_fracs, double mass=-1.0, double atoms_per_molecule=-1.0)
    Create a Material from a mapping of atom fractions.

    Parameters
    ----------
    atom_fracs : dict
        Dictionary that maps nuclides or materials to atom fractions for the
        material.  The keys may be intergers, strings, or materials. The values
        must be castable to floats.
    mass : float, optional
        This is the mass of the new stream. If the mass provided is negative
        (default -1.0) then the mass of the new stream is calculated from the
        sum of compdict's components before normalization.  If the mass here is
        positive or zero, then this mass overrides the calculated one.
    density : float, optional
        This is the density of the material.
    atoms_per_molecule : float, optional
        Number of atoms per molecule of material.  Needed to obtain proper
        scaling of molecular mass.  For example, this value for water is
        3.0.
    metadata : JSON-convertable Python object, optional
        Initial attributes to build the material with.  At the top-level this is
        usually a dictionary with string keys.  This container is used to store
        arbitrary metadata about the material.

    Returns
    -------
    mat : Material
        A material generated from atom fractions.

    Examples
    --------
    To get a material from water, based on atom fractions::

        h2o = {10010: 2.0, 'O16': 1.0}
        mat = from_atom_frac(h2o)

    Or for Uranium-Oxide, based on an initial fuel vector::

        # Define initial heavy metal
        ihm = from_atom_frac({'U235': 0.05, 'U238': 0.95})

        # Define Uranium-Oxide
        uox = {ihm: 1.0, 80160: 2.0}
        mat = from_atom_frac(uox)

    Note that the initial heavy metal was used as a key in a dictionary.
    This is possible because Materials are hashable.

    See Also
    --------
    Material.from_atom_frac : Underlying method class method.

    """
    mat = Material(metadata=metadata)
    mat.from_atom_frac(atom_fracs)

    if 0.0 <= mass:
        mat.mass = mass

    if 0.0 <= density:
        mat.density = density

    if 0.0 <= atoms_per_molecule:
        mat.atoms_per_molecule = atoms_per_molecule

    return mat



def from_activity(activities, double mass=-1.0, double density=-1.0,
                   double atoms_per_molecule=-1.0, metadata=None):
    """from_activity(activities, double mass=-1.0, double atoms_per_molecule=-1.0)
    Create a Material from a mapping of radionuclide activities. If mass < 0.0,
    it assumes the activities are supplied in units of [Bq] and sets the
    material mass to units of [grams]. Otherewise when mass >= 0.0, it
    treats the supplied activities as relative activities.

    Parameters
    ----------
    activities : dict
        Dictionary that maps radionuclides to activities for the material. The
        keys may be intergers or strings. The values must be castable to
        floats.
    mass : float, optional
        This is the mass of the new stream. If the mass provided is negative
        (default -1.0) then the mass of the new stream is calculated from the
        sum of compdict's components before normalization.  If the mass here is
        positive or zero, then this mass overrides the calculated one.
    density : float, optional
        This is the density of the material.
    atoms_per_molecule : float, optional
        Number of atoms per molecule of material.  Needed to obtain proper
        scaling of molecular mass.  For example, this value for water is
        3.0.
    metadata : JSON-convertable Python object, optional
        Initial attributes to build the material with.  At the top-level this is
        usually a dictionary with string keys.  This container is used to store
        arbitrary metadata about the material.

    Returns
    -------
    mat : Material
        A material generated from radionuclide activities.

    Examples
    --------
    To get a material of natural uranium, based on activities::

            natu = {'U234': 12223.2, 'U235': 568.648, 'U238': 12347.1}
            mat = from_activity(natu)

    See Also
    --------
    Material.from_activity : Underlying method class method.

    """
    mat = Material(metadata=metadata)
    mat.from_activity(activities)

    if 0.0 <= mass:
        mat.mass = mass

    if 0.0 <= density:
        mat.density = density

    if 0.0 <= atoms_per_molecule:
        mat.atoms_per_molecule = atoms_per_molecule

    return mat



def from_hdf5(filename, datapath, int row=-1, int protocol=1):
    """from_hdf5(char * filename, char * datapath, int row=-1, int protocol=1)
    Create a Material object from an HDF5 file.

    Parameters
    ----------
    filename : str
        Path to HDF5 file that contains the data to read in.
    datapath : str
        Path to HDF5 table or group that represents the data.
    row : int, optional
        The index of the arrays from which to read the data.  This
        ranges from 0 to N-1.  Defaults to the last element of the array.
        Negative indexing is allowed (row[-N] = row[0]).
    protocol : int, optional
        Specifies the protocol to use to read in the data.  Different
        protocols are used to represent different internal structures in
        the HDF5 file.

    Returns
    -------
    mat : Material
        A material found in the HDF5 file.

    Examples
    --------
    This method loads data into a new material::

        mat = from_hdf5("afile.h5", "/foo/bar/mat", -3)

    See Also
    --------
    Material.from_hdf5 : Underlying method class method.

    """
    cdef char * c_filename
    filename_bytes = filename.encode('UTF-8')
    c_filename = filename_bytes
    cdef char * c_datapath
    datapath_bytes = datapath.encode('UTF-8')
    c_datapath = datapath_bytes
    mat = Material()
    mat.from_hdf5(c_filename, c_datapath, row, protocol)
    return mat



def from_text(filename, double mass=-1.0, double atoms_per_molecule=-1.0, metadata=None):
    """from_text(char * filename, double mass=-1.0, double atoms_per_molecule=-1.0)
    Create a Material object from a simple text file.

    Parameters
    ----------
    filename : str
        Path to text file that contains the data to read in.
    mass : float, optional
        This is the mass of the new stream. If the mass provided is negative
        (default -1.0) then the mass of the new stream is calculated from the
        sum of compdict's components before normalization.  If the mass here is
        positive or zero, then this mass overrides the calculated one.
    atoms_per_molecule : float, optional
        Number of atoms to per molecule of material.  Needed to obtain proper
        scaling of molecular mass.  For example, this value for water is
        3.0.
    metadata : JSON-convertable Python object, optional
        Initial attributes to build the material with.  At the top-level this is
        usually a dictionary with string keys.  This container is used to store
        arbitrary metadata about the material.

    Returns
    -------
    mat : Material
        A material found in the HDF5 file.

    Examples
    --------
    This method loads data into a new Material::

        mat = from_text("natu.txt")

    See Also
    --------
    Material.from_text : Underlying method class method.

    """
    cdef char * c_filename
    filename_bytes = filename.encode('UTF-8')
    c_filename = filename_bytes
    mat = Material(metadata=metadata)

    if 0.0 <= mass:
        mat.mass = mass

    if 0.0 <= atoms_per_molecule:
        mat.atoms_per_molecule = atoms_per_molecule

    mat.from_text(c_filename)
    return mat


###########################
### Material Converters ###
###########################




# (Str, Material)
cdef class MapIterStrMaterial(object):
    cdef void init(self, cpp_map[std_string, matp] * map_ptr):
        cdef cpp_map[std_string, matp].iterator * itn = <cpp_map[std_string,
                matp].iterator *> malloc(sizeof(map_ptr.begin()))
        itn[0] = map_ptr.begin()
        self.iter_now = itn

        cdef cpp_map[std_string, matp].iterator * ite = <cpp_map[std_string,
                matp].iterator *> malloc(sizeof(map_ptr.end()))
        ite[0] = map_ptr.end()
        self.iter_end = ite

    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_map[std_string, matp].iterator inow = deref(self.iter_now)
        cdef cpp_map[std_string, matp].iterator iend = deref(self.iter_end)

        if inow != iend:
            pyval = str(<char *> deref(inow).first.c_str())
        else:
            raise StopIteration

        inc(deref(self.iter_now))
        return pyval


cdef class _MapStrMaterial:
    def __cinit__(self, new_map=True, bint free_map=True):
        cdef std_string s
        cdef cpp_pair[std_string, matp] item

        # Cache needed to prevent Python from deref'ing
        # pointers before their time.
        self._cache = {}

        # Decide how to init map, if at all
        if isinstance(new_map, _MapStrMaterial):
            self.map_ptr = (<_MapStrMaterial> new_map).map_ptr
            self._cache = (<_MapStrMaterial> new_map)._cache
        elif hasattr(new_map, 'items'):
            self.map_ptr = new cpp_map[std_string, matp]()
            for key, value in new_map.items():
                #s = std_string(key)
                #item = cpp_pair[std_string, matp](s, (<_Material>
                # value).mat_pointer)
                #self.map_ptr.insert(item)
                self[key] = value
        elif hasattr(new_map, '__len__'):
            self.map_ptr = new cpp_map[std_string, matp]()
            for i in new_map:
                #s = std_string(i[0])
                #item = cpp_pair[std_string, matp](s, (<_Material>
                # i[1]).mat_pointer)
                #self.map_ptr.insert(item)
                self[i[0]] = i[1]
        elif bool(new_map):
            self.map_ptr = new cpp_map[std_string, matp]()

        # Store free_map
        self._free_map = free_map

    def __dealloc__(self):
        if self._free_map:
            del self.map_ptr

    def __contains__(self, key):
        cdef std_string s
        if isinstance(key, str):
            s = std_string(<char *> key)
        else:
            return False

        if 0 < self.map_ptr.count(s):
            return True
        else:
            return False

    def __len__(self):
        return self.map_ptr.size()

    def __iter__(self):
        cdef MapIterStrMaterial mi = MapIterStrMaterial()
        mi.init(self.map_ptr)
        return mi

    def __getitem__(self, key):
        cdef std_string s
        cdef _Material pymat

        if isinstance(key, basestring):
            key = key.encode()
            s = std_string(<char *> key)
        else:
            raise TypeError("Only string keys are valid.")

        if 0 < self.map_ptr.count(s):
            if key not in self._cache:
                pymat = Material(nucvec=None, free_mat=False)
                pymat.mat_pointer = deref(self.map_ptr)[s]
                self._cache[key] = pymat
            return self._cache[key]
        else:
            raise KeyError(repr(key) + " not in map.")

    def __setitem__(self, key, value):

        cdef char * c_key
        key_bytes = key.encode('UTF-8')
        c_key = key_bytes
        cdef std_string s = std_string(c_key)
        if not isinstance(value, _Material):
            raise TypeError("may only set materials into this mapping.")
        cdef cpp_pair[std_string, matp] item = cpp_pair[std_string, matp](s,
                (<_Material> value).mat_pointer)
        self.map_ptr.insert(item)
        self._cache[c_key] = value

    def __delitem__(self, char * key):
        cdef std_string s
        if key in self:
            s = std_string(key)
            self.map_ptr.erase(s)
            del self._cache[key]


class MapStrMaterial(_MapStrMaterial, collectionsAbc.MutableMapping):
    """Wrapper class for C++ standard library maps of type <string, Material
    \*>.  Provides dictionary like interface on the Python level.

    Parameters
    ----------
    new_map : bool or dict-like
        Boolean on whether to make a new map or not, or dict-like object
        with keys and values which are castable to the appropriate type.
    free_map : bool
        Flag for whether the pointer to the C++ map should be deallocated
        when the wrapper is dereferenced.

    """
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "{" + ", ".join(["{0}: {1}".format(repr(key), value) for key, value in self.items()]) + "}"



class MultiMaterial(collectionsAbc.MutableMapping):
    """ This class is serves as a way of storing a collection of materials.
    There sole argument of this function is a dictionary with material
    objects and keys and vol/mass fractions as values. There are two
    main uses cases. A collection of materials can be mixed together
    by volume of mass. Alternatively, a collection of materials can be
    used to describe a materials with multiple densities. In this latter
    case, the dict values are irrevelant"""
    def __init__(self, mats):
        # This function reads a dict of materials and either mass or volume factions,
        # then normalizes the fraction and assigns the dict as an attribute of self.
        # Normalize Mixture fractions
        # First calculate the normalization factor:
        norm = 0
        for mat, mix_frac in mats.items():
            norm = norm + mix_frac
        # Now divide all mixture fraction by the normalization factor
        # and assign them as MultiMaterial attributes
        for mat, mix_frac in mats.items():
            mats[mat] = mix_frac / norm
            mat.mass = 1  # set all mass to 1 for mixing
        self._mats = mats

    def __getitem__(self, key):
        return self._mats[key]

    def __setitem__(self, key, value):
        self._mats[key] = value

    def __add__(self, other):
        pass

    def __delitem__(self, key):
        pass

    def __iter__(self):
        pass

    def  __len__(self):
        pass

    def mix_by_mass(self):
        """This function reads in a python dict of materials and mass fractions
        then mixes the material by mass fractions and returns a material of mass=1.
        """
        total = 0
        total_mass_frac = 0
        mix = Material()
        for mat, mat_frac in self._mats.items():
            mix = mix + mat*mat_frac
            total_mass_frac += mat_frac
            total += (mat_frac/mat.density)
        mix.mass = 1
        mix.density = total_mass_frac/total

        return mix

    def mix_by_volume(self):
        """This function reads in a python dict of materials and volume fractions
        then mixes the material by volume fractions and returns a material of mass=1.
        """
        total = 0
        mix = Material()
        for mat, mat_frac in self._mats.items():
            mix = mix + mat*mat_frac*mat.density
            total += (mat.density*mat_frac)
        mix.mass = 1
        mix.density = total
        return mix


def mats_latex_table(mats, labels=None, align=None, format=".5g"):
    if align is None:
        align = '|l|' + 'c|'*len(mats)
    if labels is None:
        labels = []
    if isinstance(format, basestring):
        format = [format] * len(mats)

    nucs = set()
    colnames = ["Nuclide"]
    for i, mat in enumerate(mats):
        nucs |= set(mat.comp.keys())
        name = mat.metadata['name'] if 'name' in mat.metadata else "mat{0}".format(i)
        colnames.append(name)
    nucs = sorted(nucs)
    colnames = labels + colnames[len(labels):]
    colnames = " & ".join([r"\bf{"+n+'}' for n in colnames]) + r" \\ " + "\n"

    tab = r"\begin{tabular}{" + align + "}\n\\hline\n" + colnames + "\\hline\n"
    for nuc in nucs:
        tab += nucname.name(nuc)
        for mat, f in zip(mats, format):
            val = mat[nuc] if nuc in mat else 0.0
            tab += " & {v:{f}}".format(v=val, f=f)
        tab += r" \\ " + "\n\\hline\n"
    tab += "\\end{tabular}\n"
    return tab


ensure_material = lambda m: m if isinstance(m, Material) else Material(m)

