"""Python wrapper for material library."""
# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
#from cython.operator cimport reference as ref
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport malloc, free

# Python imports
import collections

# local imports
cimport std
cimport cpp_material
cimport pyne.stlconverters as conv
import pyne.stlconverters as conv

cimport cpp_jsoncpp
cimport jsoncpp
import jsoncpp

cimport pyne.nucname as nucname
import pyne.nucname as nucname
import os


cdef cpp_map[int, double] dict_to_comp(dict nucvec):
    """Converts a dictionary with arbitraily-typed keys to component map."""
    cdef int key_zz
    cdef cpp_map[int, double] comp = cpp_map[int, double]()

    for key, value in nucvec.items():
        if isinstance(key, int):
            comp[key] = value
        else:
            key_zz = nucname.zzaaam(key)
            comp[key_zz] = value

    return comp


cdef class _Material:

    def __cinit__(self, nucvec=None, double mass=-1.0, double density=-1.0,
                  double atoms_per_mol=-1.0, attrs=None, bint free_mat=True,
                  *args, **kwargs):
        """Material C++ constuctor."""
        cdef cpp_map[int, double] comp
        cdef jsoncpp.Value cattrs = jsoncpp.Value({} if attrs is None else attrs)

        if isinstance(nucvec, dict):
            # Material from dict
            comp = dict_to_comp(nucvec)
            self.mat_pointer = new cpp_material.Material(
                    comp, mass, density, atoms_per_mol, deref(cattrs._inst))

        elif isinstance(nucvec, basestring):
            # Material from file
            self.mat_pointer = new cpp_material.Material(
                    <char *> nucvec, mass, density, atoms_per_mol, deref(cattrs._inst))

        elif (nucvec is None):
            if free_mat:
                # Make empty mass stream
                self.mat_pointer = new cpp_material.Material(comp, 
                                        mass, atoms_per_mol, deref(cattrs._inst))
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

    property atoms_per_mol:
        def __get__(self):
            return self.mat_pointer.atoms_per_mol

        def __set__(self, double value):
            self.mat_pointer.atoms_per_mol = value

    property attrs:
        def __get__(self):
            cdef jsoncpp.Value val = jsoncpp.Value(view=True)
            val._inst = &self.mat_pointer.attrs
            return val

        def __set__(self, value):
            cdef jsoncpp.Value val = jsoncpp.Value(value)
            val._view = True
            self.mat_pointer.attrs = deref(val._inst)

    #
    # Class Methods
    #

    def norm_comp(self):
        """Normalizes the composition, preserving the mass of the nuclide
        vector as mass.

        """
        self.mat_pointer.norm_comp()


    def from_hdf5(self, char * filename, char * datapath, int row=-1,
                  int protocol=1):
        """from_hdf5(char * filename, char * datapath, int row=-1, int protocol=1)
        Initialize a Material object from an HDF5 file.

        Parameters
        ----------
        filename : str
            Path to HDF5 file that contains the data to read in.
        datapath : str
            Path to HDF5 table or group that represents the data.
            In the example below, datapath = "/material".
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
        information and an array that stores the nuclides (zzaaam) which index
        the comp array::

            file.h5 (file)
                |-- material (table)
                    |-- name (string col, len 20)
                    |-- mass (double col)
                    |-- atoms_per_mol (double col)
                    |-- comp (double array col, len of nuc_zz)
                |-- nuc_zz (int array)

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
        self.mat_pointer.from_hdf5(filename, datapath, row, protocol)


    def write_hdf5(self, filename, datapath="/material", nucpath="/nuc_zz",
                   row=-0.0, chunksize=100):
        """write_hdf5(filename, datapath="/material", nucpath="/nuc_zz", row=-0.0, chunksize=100)
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
            Path to zzaaam array of nuclides to write out.  If this array does
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
        self.mat_pointer.write_hdf5(filename, datapath, nucpath, row, chunksize)


    def from_text(self, char * filename):
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
        self.mat_pointer.from_text(filename)


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
        self.mat_pointer.write_text(filename)




    def normalize(self):
        """This convenience method normalizes the mass stream by setting its
        mass = 1.0.

        """
        self.mat_pointer.normalize()


    def mult_by_mass(self):
        """This multiplies multiplies comp by mass and returns the resultant
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


    def molecular_weight(self, atoms_per_mol=-1.0):
        """molecular_weight(atoms_per_mol=-1.0)
        This method returns the molecular weight of the comp of this
        material.

        Parameters
        ----------
        atoms_per_mol : double, optional
            Number of atoms to per molecule of material.  Needed to obtain
            proper scaling.  For example, this value for water is 3.0.

        Returns
        -------
        mol_weight : float
            Molecular weight in [amu].

        """
        return self.mat_pointer.molecular_weight(atoms_per_mol)


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
            Elements and nuctopes to be taken from current stream.
            Members of this list must be integers.  For example, [92, 942390]
            would take all uranium atoms and Pu-239.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only
            has the members given in nuc_sequence.  The mass of the submaterial
            is calculated based on the weight fraction composition and mass
            of the original mass stream.

        Notes
        -----
        The input here is seen as a suggestion and so no error is raised if a
        nuclide is asked for via nuc_sequence that is not present in the
        original material.

        """
        # Make an nuctopic set
        cdef cpp_set[int] nuc_set = nucname.zzaaam_set(nuc_sequence)

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
            Elements and nuctopes to be taken from current stream.
            Members of this list must be integers.  For example, [92, 942390]
            would take all uranium atoms and Pu-239.
        value : float
            Mass value to set all nuclides in sequence to on the material.

        Returns
        -------
        submaterial : Material
            A new material object whose members in nuc_sequence have the
            cooresponding mass value.  The mass of the submaterial is
            calculated based on the weight fraction composition and mass of the
            original material.

        """
        # Make an nuctopic set
        cdef cpp_set[int] nuc_set = nucname.zzaaam_set(nuc_sequence)

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
            the weight fraction composition and mass of the original material.

        Notes
        -----
        The input here is seen as a suggestion and so no error is raised if a
        nuclide is asked for via nuc_sequence that is not present in the
        original material.

        """
        # Make an nuctopic set
        cdef cpp_set[int] nuc_set = nucname.zzaaam_set(nuc_sequence)

        # Make new python version of this material
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.del_mat(nuc_set)
        return pymat


    def sub_range(self, lower=0, upper=10000000):
        """sub_range(lower=0, upper=10000000)
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
            clower = nucname.zzaaam(lower)

        if isinstance(upper, int):
            cupper = upper
        else:
            cupper = nucname.zzaaam(upper)

        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_range(clower, cupper)
        return pymat


    def set_range(self, lower=0, upper=10000000, value=0.0):
        """set_range(lower=0, upper=10000000, value=0.0)
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
            clower = nucname.zzaaam(lower)

        if isinstance(upper, int):
            cupper = upper
        else:
            cupper = nucname.zzaaam(upper)

        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.set_range(clower, cupper, <double> value)
        return pymat


    def del_range(self, lower=0, upper=10000000):
        """del_range(lower=0, upper=10000000)
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
            clower = nucname.zzaaam(lower)

        if isinstance(upper, int):
            cupper = upper
        else:
            cupper = nucname.zzaaam(upper)

        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.del_range(clower, cupper)
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
        pymat.mat_pointer[0] = self.mat_pointer.sub_u()
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
        pymat.mat_pointer[0] = self.mat_pointer.sub_pu()
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
                key_zz = <int> key
                if 0 == af.count(key_zz):
                    af[key_zz] = 0.0
                af[key_zz] = af[key_zz] + val
            elif isinstance(key, basestring):
                key_zz = nucname.zzaaam(key)
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
        pymat.mat_pointer[0] = self.mat_pointer[0] / y
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
        return self.__div__(y)


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
            key_zz = nucname.zzaaam(key)
            return self[key_zz]

        # Get slice-based sub-material
        elif isinstance(key, slice):
            lower = key.start
            if lower is None:
                lower = 0

            upper = key.stop
            if upper is None:
                upper = 10000000

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
            key_zz = nucname.zzaaam(key)
            self[key_zz] = value

        # Set slice-based sub-material
        elif isinstance(key, slice):
            lower = key.start
            if lower is None:
                lower = 0

            upper = key.stop
            if upper is None:
                upper = 10000000

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
            new_matp = new cpp_material.Material(mbm.map_ptr[0], -1.0, -1.0)
            self.mat_pointer = new_matp
            self._comp = None

        # Remove single string-key
        elif isinstance(key, basestring):
            key_zz = nucname.zzaaam(key)
            del self[key_zz]

        # Remove slice-based sub-material
        elif isinstance(key, slice):
            lower = key.start
            if lower is None:
                lower = 0

            upper = key.stop
            if upper is None:
                upper = 10000000

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



class Material(_Material, collections.MutableMapping):
    """Material composed of nuclides.

    Parameters
    ----------
    comp : dict or str
        This is the input nuclide component dictionary.  This dictionary need
        not be normalized; Material initialization will automatically
        renormalize the stream.  Thus the comp simply is a dictionary of
        relative weights.  The keys of comp must be integers representing
        nuclides in zzaaam-form.  The values are floats for each nuclide's
        weight fraction. If a string is provided instead of a dictionary, then
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
    atoms_per_mol : float, optional
        Number of atoms to per molecule of material.  Needed to obtain proper
        scaling of molecular weights.  For example, this value for water is
        3.0.
    attrs : JSON-convertable Python object, optional
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
        header += ["density= {0}".format(self.density)]
        header += ["atoms per molecule = {0}".format(self.atoms_per_mol)]
        if self.attrs.isobject():
            for key, value in self.attrs.items():
                header += ["{0} = {1}".format(key, value)]
        header += ['-' * max([len(h) for h in header])]
        header = "\n".join(header) + "\n"

        s = header + "\n".join(["{0:<7}{1}".format(
                nucname.name(key), value) for key, value in self.comp.items()])
        return s

    def __repr__(self):
        return "pyne.material.Material({0}, {1}, {2}, {3}, {4})".format(
                repr(self.comp), self.mass, self.density, self.atoms_per_mol, repr(self.attrs))


#####################################
### Material generation functions ###
#####################################

def from_atom_frac(atom_fracs, double mass=-1.0, double
                   atoms_per_mol=-1.0, attrs=None):
    """from_atom_frac(atom_fracs, double mass=-1.0, double atoms_per_mol=-1.0)
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
    atoms_per_mol : float, optional
        Number of atoms per molecule of material.  Needed to obtain proper
        scaling of molecular weights.  For example, this value for water is
        3.0.
    attrs : JSON-convertable Python object, optional
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
    mat = Material(attrs=attrs)
    mat.from_atom_frac(atom_fracs)

    if 0.0 <= mass:
        mat.mass = mass

    if 0.0 <= atoms_per_mol:
        mat.atoms_per_mol = atoms_per_mol

    return mat



def from_hdf5(char * filename, char * datapath, int row=-1, int protocol=1):
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
    mat = Material()
    mat.from_hdf5(filename, datapath, row, protocol)
    return mat



def from_text(char * filename, double mass=-1.0, double atoms_per_mol=-1.0, attrs=None):
    """from_text(char * filename, double mass=-1.0, double atoms_per_mol=-1.0)
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
    atoms_per_mol : float, optional
        Number of atoms to per molecule of material.  Needed to obtain proper
        scaling of molecular weights.  For example, this value for water is
        3.0.
    attrs : JSON-convertable Python object, optional
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
    mat = Material(attrs=attrs)

    if 0.0 <= mass:
        mat.mass = mass

    if 0.0 <= atoms_per_mol:
        mat.atoms_per_mol = atoms_per_mol

    mat.from_text(filename)
    return mat





###########################
### Material Converters ###
###########################

# <string, Material *>

cdef cpp_map[std.string, matp] dict_to_map_str_matp(dict pydict):
    cdef _Material pymat
    cdef cpp_material.Material * cpp_matp
    cdef cpp_map[std.string, matp] cppmap = cpp_map[std.string, matp]()
    cdef cpp_pair[std.string, matp] item

    for key, value in pydict.items():
        pymat = value
        cpp_matp = pymat.mat_pointer
        #cppmap[std.string(key)] = cpp_matp
        item = cpp_pair[std.string, matp](std.string(key), cpp_matp)
        cppmap.insert(item)

    return cppmap


cdef dict map_to_dict_str_matp(cpp_map[std.string, matp] cppmap):
    pydict = {}
    cdef _Material pymat
    cdef cpp_map[std.string, matp].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pymat = Material()
        pymat.mat_pointer[0] = deref(deref(mapiter).second)
        pydict[deref(mapiter).first.c_str()] = pymat
        inc(mapiter)

    return pydict



# (Str, Material)
cdef class MapIterStrMaterial(object):
    cdef void init(self, cpp_map[std.string, matp] * map_ptr):
        cdef cpp_map[std.string, matp].iterator * itn = <cpp_map[std.string,
                matp].iterator *> malloc(sizeof(map_ptr.begin()))
        itn[0] = map_ptr.begin()
        self.iter_now = itn

        cdef cpp_map[std.string, matp].iterator * ite = <cpp_map[std.string,
                matp].iterator *> malloc(sizeof(map_ptr.end()))
        ite[0] = map_ptr.end()
        self.iter_end = ite

    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_map[std.string, matp].iterator inow = deref(self.iter_now)
        cdef cpp_map[std.string, matp].iterator iend = deref(self.iter_end)

        if inow != iend:
            pyval = str(deref(inow).first.c_str())
        else:
            raise StopIteration

        inc(deref(self.iter_now))
        return pyval


cdef class _MapStrMaterial:
    def __cinit__(self, new_map=True, bint free_map=True):
        cdef std.string s
        cdef cpp_pair[std.string, matp] item

        # Cache needed to prevent Python from deref'ing
        # pointers before their time.
        self._cache = {}

        # Decide how to init map, if at all
        if isinstance(new_map, _MapStrMaterial):
            self.map_ptr = (<_MapStrMaterial> new_map).map_ptr
            self._cache = (<_MapStrMaterial> new_map)._cache
        elif hasattr(new_map, 'items'):
            self.map_ptr = new cpp_map[std.string, matp]()
            for key, value in new_map.items():
                #s = std.string(key)
                #item = cpp_pair[std.string, matp](s, (<_Material>
                # value).mat_pointer)
                #self.map_ptr.insert(item)
                self[key] = value
        elif hasattr(new_map, '__len__'):
            self.map_ptr = new cpp_map[std.string, matp]()
            for i in new_map:
                #s = std.string(i[0])
                #item = cpp_pair[std.string, matp](s, (<_Material>
                # i[1]).mat_pointer)
                #self.map_ptr.insert(item)
                self[i[0]] = i[1]
        elif bool(new_map):
            self.map_ptr = new cpp_map[std.string, matp]()

        # Store free_map
        self._free_map = free_map

    def __dealloc__(self):
        if self._free_map:
            del self.map_ptr

    def __contains__(self, key):
        cdef std.string s
        if isinstance(key, str):
            s = std.string(key)
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
        cdef std.string s
        cdef _Material pymat

        if isinstance(key, basestring):
            s = std.string(key)
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

    def __setitem__(self, char * key, value):
        cdef std.string s = std.string(key)
        if not isinstance(value, _Material):
            raise TypeError("may only set materials into this mapping.")
        cdef cpp_pair[std.string, matp] item = cpp_pair[std.string, matp](s,
                (<_Material> value).mat_pointer)
        self.map_ptr.insert(item)
        self._cache[key] = value

    def __delitem__(self, char * key):
        cdef std.string s
        if key in self:
            s = std.string(key)
            self.map_ptr.erase(s)
            del self._cache[key]


class MapStrMaterial(_MapStrMaterial, collections.MutableMapping):
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



class MultiMaterial(collections.MutableMapping):

	def __init__(self, mats):
		self._mats = mats

	def __getitem__(self, key):
		return self._mats[key]

	def __setitem__(self, key, value):
		self._mats[key] = value

	def write_hdf5(self):
		pass

	def write_iter(self):
		pass

	def __add__(self, other):
		pass
