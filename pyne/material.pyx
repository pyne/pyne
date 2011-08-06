"""Python wrapper for material library."""
# Cython imports
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

cimport pyne.nucname as nucname
import pyne.nucname as nucname
import os



cdef class _Material:

    def __cinit__(self, nucvec=None, float mass=-1.0, char * name=''):
        """Material C++ constuctor."""
        cdef cpp_map[int, double] comp

        if isinstance(nucvec, dict):
            # Material from dict
            comp = conv.dict_to_map_int_dbl(nucvec)
            self.mat_pointer = new cpp_material.Material(comp, mass, std.string(name))

        elif isinstance(nucvec, basestring):
            # Material from file
            self.mat_pointer = new cpp_material.Material(<char *> nucvec, mass, std.string(name))

        elif nucvec is None:
            # Empty mass stream
            self.mat_pointer = new cpp_material.Material()

        else:
            # Bad Material 
            raise TypeError("The mass stream nucvec must be a dict, str, or None.")


    def __dealloc__(self):
        """Material C++ destructor."""
        del self.mat_pointer


    #
    # Class Attributes
    #

    property comp:
        def __get__(self):
            cdef conv._MapProxyIntDouble comp_proxy = conv.MapProxyIntDouble()
            comp_proxy.map_ptr = new cpp_map[int, double](self.mat_pointer.comp)
            return comp_proxy

        def __set__(self, dict comp):
            self.mat_pointer.comp = conv.dict_to_map_int_dbl(comp)


    property mass:
        def __get__(self):
            return self.mat_pointer.mass

        def __set__(self, double mass):
            self.mat_pointer.mass = mass


    property name:
        def __get__(self):
            cdef std.string mat_name = self.mat_pointer.name
            return mat_name.c_str()

        def __set__(self, char * name):
            self.mat_pointer.name = std.string(name)

    #
    # Class Methods
    #

    def norm_comp(self):
        """Normalizes the composition, preserving the mass of the nuclide vector as mass."""
        self.mat_pointer.norm_comp()


    def load_from_hdf5(self, char * filename, char * groupname, int row=-1):
        """Initialize a Material object from an HDF5 file.

        Parameters
        ----------
        filename : str
            Path to HDF5 file that contains the data to read in.    
        groupname : str 
            Path to HDF5 group that represents the data. 
            In the above example, groupname = "/Material".    
        row : int, optional 
            The index of the arrays from which to read the data.  This 
            ranges from 0 to N-1.  Defaults to the last element of the array.
            Negative indexing is allowed (row[-N] = row[0]).

        Notes
        -----
        The HDF5 representation of a Material is a group that holds several 
        extendable array datasets.  One array is entitled "Mass" while the other datasets
        are nuclide names in LLAAAM form ("U235", "NP237", *etc*).  For example::

            File.h5 (file)
                |-- Material (group)
                    |-- Mass (array)
                    |-- H1 (array)
                    |-- O16 (array)
                    |-- U235 (array)
                    |-- PU239 (array)
                    |-- ...

        The arrays are all of length N, where each row typically represents a different 
        fuel cycle pass.  The sum of all of the nuclide arrays should sum to one, like 
        Material.comp. 

        Examples
        --------
        This method loads data into a pre-existing :class:`Material`.  
        Initialization is therefore a two-step process::

            mat = Material()
            mat.load_from_hdf5("afile.h5", "/foo/bar/mat", -3)
        """
        self.mat_pointer.load_from_hdf5(filename, groupname, row)


    def load_from_text(self, char * filename):
        """Initialize a Material object from a simple text file.

        Parameters
        ----------
        filename : str 
            Path to HDF5 file that contains the data to read in.    


        Notes
        -----
        The text representation of Materials are nuclide identifiers in the 
        first column and mass or weight values in the second column.  For example, 
        for natural uranium::

            922340  0.000055
            U235    0.00720
            92238   0.992745

        Data in this file must be whitespace separated.  Any valid nuclide naming
        scheme may be used for any nuctope.

        Examples
        --------
        This method loads data into a pre-existing Material.  
        Initialization is therefore a two-step process::

            mat = Material()
            mat.load_from_text("natu.h5")

        This method is most often called implicitly by the Material constructor.
        """
        self.mat_pointer.load_from_text(filename)



    def normalize(self):
        """This convenience method normalizes the mass stream by setting its mass = 1.0."""
        self.mat_pointer.normalize()


    def mult_by_mass(self):
        """This multiplies multiplies comp by mass and returns the resultant nuctopic vector.

        Returns
        -------
        nucvec : dict
            For a Material mat, 

            .. math:: \\mbox{nucvec[nuc]} = \\mbox{mat.comp[nuc]} \\times \\mbox{mat.mass}
        """
        cdef conv._MapProxyIntDouble nucvec_proxy = conv.MapProxyIntDouble()
        nucvec_proxy.map_ptr = new cpp_map[int, double](self.mat_pointer.mult_by_mass())
        return nucvec_proxy


    def molecular_weight(self):
        """This method returns the molecular weight of the comp of this material.  Note that this is 
        only a rough estimate since this function is not yet coupled with measured atomic weights.

        Returns
        -------
        mol_weight : float 
            Molecular weight in [amu]."""
        return self.mat_pointer.molecular_weight()


    #
    # submaterial Methods
    #

    def sub_mat(self, nuc_sequence, char * name=""):
        """Grabs a subset of the material and returns a new material comprised of only
        the specified nuclides.  

        Parameters
        ----------
        nuc_sequence : sequence
            Elements and nuctopes to be taken from current stream.
            Members of this list must be integers.  For example, [92, 942390]
            would take all uranium atoms and Pu-239.  
        name : str, optional
            The name of the submaterial.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only 
            has the members given in nuc_sequence.  The mass of the submaterial
            is calculated based on the weight fraction composition and mass
            of the original mass stream.

        Notes
        -----
        The input here is seen as a suggestion and so no error is raised if a nuclide 
        is asked for via nuc_sequence that is not present in the original material.
        """
        # Make an nuctopic set 
        cdef cpp_set[int] nuc_set = nucname.zzaaam_set(nuc_sequence)

        # Make new python version of this material
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_mat(nuc_set, std.string(name))
        return pymat


    def set_mat(self, nuc_sequence, value, char * name=""):
        """Sets a subset of the material to a new value and returns a new material.

        Parameters
        ----------
        nuc_sequence : sequence
            Elements and nuctopes to be taken from current stream.
            Members of this list must be integers.  For example, [92, 942390]
            would take all uranium atoms and Pu-239.
        value : float
            Mass value to set all nuclides in sequence to on the material. 
        name : str, optional
            The name of the submaterial.

        Returns
        -------
        submaterial : Material
            A new material object whose members in nuc_sequence have the cooresponding 
            mass value.  The mass of the submaterial is calculated based on the weight 
            fraction composition and mass of the original material.
        """
        # Make an nuctopic set 
        cdef cpp_set[int] nuc_set = nucname.zzaaam_set(nuc_sequence)

        # Make new python version of this material
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.set_mat(nuc_set, <double> value, std.string(name))
        return pymat


    def del_mat(self, nuc_sequence, char * name=""):
        """Removes a subset of the material and returns a new material comprised of only
        the non-specified nuclides.  

        Parameters
        ----------
        nuc_sequence : sequence
            Nuclides to be taken out of the current material.
        name : str, optional
            The name of the submaterial.

        Returns
        -------
        submaterial : Material
            A new material object that only has the members not given in nuc_sequence. 
            The mass of the submaterial is calculated based on the weight fraction 
            composition and mass of the original material.

        Notes
        -----
        The input here is seen as a suggestion and so no error is raised if a nuclide 
        is asked for via nuc_sequence that is not present in the original material.
        """
        # Make an nuctopic set 
        cdef cpp_set[int] nuc_set = nucname.zzaaam_set(nuc_sequence)

        # Make new python version of this material
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.del_mat(nuc_set, std.string(name))
        return pymat


    def sub_range(self, lower=0, upper=10000000, char * name=""):
        """Grabs a sub-material from this mat based on a range [lower, upper) of values.

        Parameters
        ----------
        lower : nuclide-name, optional
            Lower bound on nuclide range.
        upper : nuclide-name, optional
            Upper bound on nuclide range.
        name : str, optional 
            The name of the submaterial.

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
        pymat.mat_pointer[0] = self.mat_pointer.sub_range(clower, cupper, std.string(name))
        return pymat
        

    def set_range(self, lower=0, upper=10000000, value=0.0, char * name=""):
        """Sets a sub-material from this mat based on a range [lower, upper) to a new 
        mass weight value.

        Parameters
        ----------
        lower : nuclide-name, optional
            Lower bound on nuclide range.
        upper : nuclide-name, optional
            Upper bound on nuclide range.
        value : float
            Mass value to set all nuclides on the range to on the material. 
        name : str, optional 
            The name of the submaterial.

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
        pymat.mat_pointer[0] = self.mat_pointer.set_range(clower, cupper, <double> value, std.string(name))
        return pymat
        

    def del_range(self, lower=0, upper=10000000, char * name=""):
        """Remove a range [lower, upper) of nuclides from this material and returns a submaterial.

        Parameters
        ----------
        lower : nuclide-name, optional
            Lower bound on nuclide range.
        upper : nuclide-name, optional
            Upper bound on nuclide range.
        name : str, optional 
            The name of the submaterial.

        Returns
        -------
        submaterial : Material 
            A new mass stream object that does not contain nuclides on the given range.
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
        pymat.mat_pointer[0] = self.mat_pointer.del_range(clower, cupper, std.string(name))
        return pymat


    def sub_u(self, char * name=""):
        """Convenience method that gets the Uranium portion of a mass stream.

        Parameters
        ----------
        name : str, optional 
            The name of the submaterial.

        Returns
        -------
        submaterial : Material 
            A new mass stream object that only has Uranium members. 
        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_u(std.string(name))
        return pymat
        

    def sub_pu(self, char * name=""):
        """Convenience method that gets the Plutonium portion of a mass stream.

        Parameters
        ----------
        name : str, optional 
            The name of the submaterial.

        Returns
        -------
        submaterial : Material 
            A new mass stream object that only has Plutonium members. 
        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_pu(std.string(name))
        return pymat
        

    def sub_lan(self, char * name=""):
        """Convenience method that gets the Lanthanide portion of a mass stream.

        Parameters
        ----------
        name : str, optional 
            The name of the submaterial.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has Lanthanide members. 
        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_lan(std.string(name))
        return pymat
        

    def sub_act(self, char * name=""):
        """Convenience method that gets the Actinide portion of a mass stream.

        Parameters
        ----------
        name : str, optional
            The name of the submaterial.

        Returns
        -------
        submaterial : Material
            A new mass stream object that only has Actinide members. 
        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_act(std.string(name))
        return pymat
        

    def sub_tru(self, char * name=""):
        """Convenience method that gets the Transuranic portion of a mass stream.

        Parameters
        ----------
        name : str, optional 
            The name of the submaterial.

        Returns
        -------
        submaterial : Material 
            A new mass stream object that only has Transuranic members. 
        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_tru(std.string(name))
        return pymat
        

    def sub_ma(self, char * name=""):
        """Convenience method that gets the Minor Actinide portion of a mass stream.

        Parameters
        ----------
        name : str, optional 
            The name of the submaterial.

        Returns
        -------
        submaterial : Material 
            A new mass stream object that only has Minor Actinide members. 
        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_ma(std.string(name))
        return pymat
        

    def sub_fp(self, char * name=""):
        """Convenience method that gets the Fission Product portion of a mass stream.

        Parameters
        ----------
        name : str, optional 
            The name of the submaterial.

        Returns
        -------
        submaterial : Material 
            A new mass stream object that only has Fission Product members. 
        """
        cdef _Material pymat = Material()
        pymat.mat_pointer[0] = self.mat_pointer.sub_fp(std.string(name))
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
        pymat.mat_pointer[0] =  (<_Material> x).mat_pointer[0] + (<_Material> y).mat_pointer[0]
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
            raise TypeError("key {0} is of unsupported type {1}".format(repr(key), type(key)))


    def __setitem__(self, key, double value):
        cdef _Material new_mat
        cdef matp new_matp 
        cdef conv._MapProxyIntDouble mbm 
        cdef int key_zz 
        cdef cpp_map[int, double].iterator mbmiter, mbmend

        # Set single integer-key
        if isinstance(key, int):
            mbm = self.mult_by_mass()
            mbm.map_ptr[0][key] = value
            new_matp = new cpp_material.Material(mbm.map_ptr[0], -1.0, self.mat_pointer.name)
            self.mat_pointer = new_matp

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
            new_mat = self.set_range(lower, upper, value, self.name)
            self.mat_pointer[0] = new_mat.mat_pointer[0]

        # Set sequance-based sub-material
        elif hasattr(key, '__len__'):
            new_mat = self.set_mat(key, value, self.name)
            self.mat_pointer[0] = new_mat.mat_pointer[0]

        # Fail-Yurt
        else:
            raise TypeError("key {0} is of unsupported type {1}".format(repr(key), type(key)))


    def __delitem__(self, key):
        cdef _Material new_mat
        cdef matp new_matp
        cdef conv._MapProxyIntDouble mbm
        cdef int key_zz
        cdef cpp_map[int, double].iterator mbmiter, mbmend

        # Remove single key
        if isinstance(key, int):
            if 0 == self.mat_pointer.comp.count(key):
                return
            mbm = self.mult_by_mass()
            mbm.map_ptr.erase(<int> key)
            new_matp = new cpp_material.Material(mbm.map_ptr[0], -1.0, self.mat_pointer.name)
            self.mat_pointer = new_matp

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
            new_mat = self.del_range(lower, upper, self.name)
            self.mat_pointer[0] = new_mat.mat_pointer[0]

        # Remove sequance-based sub-material
        elif hasattr(key, '__len__'):
            new_mat = self.del_mat(key, self.name)
            self.mat_pointer[0] = new_mat.mat_pointer[0]

        # Fail-Yurt
        else:
            raise TypeError("key {0} is of unsupported type {1}".format(repr(key), type(key)))


    def __iter__(self):
        mbm = self.mult_by_mass()
        mbm_iter = iter(mbm)
        return mbm_iter





class Material(_Material, collections.MutableMapping):
    """Material composed of nuclides.

    Parameters
    ----------
    comp : dict or str 
        This is the input nuclide component dictionary.
        This dictionary need not be normalized; Material initialization will
        automatically renormalize the stream.  Thus the comp simply is a dictionary
        of relative weights.  The keys of comp must be integers representing
        nuclides in zzaaam-form.  The values are floats for each nuclide's weight 
        fraction. If a string is provided instead of a dictionary, then Material will
        read in the comp vector from a file at the string's location.  This  
        either plaintext or hdf5 files. If no comp is provided, an empty Material 
        object is constructed.
    mass : float, optional
        This is the mass of the new stream. If the mass provided
        is negative (default -1.0) then the mass of the new stream is calculated from 
        the sum of compdict's components before normalization.  If the mass here
        is positive or zero, then this mass overrides the calculated one.
    name : str, optional
        A string label for the material.  Helpful for large numbers of 
        streams. Default ''.

    """
    def __str__(self):
        title = "Material: {0}".format(self.name)
        underline = '-' * len(title)
        title += "\nmass = {0}\n{1}\n".format(self.mass, underline)
        s = title + "\n".join(["{0:<7}{1}".format(nucname.LLAAAM(key), value) for key, value in self.comp.items()])
        return s

    def __repr__(self):
        return "pyne.material.Material({0}, {1}, {2})".format(repr(self.comp), self.mass, repr(self.name))


###########################
### Material Converters ###
###########################

# <string, Material *>

cdef cpp_map[std.string, matp] dict_to_map_str_matp(dict pydict):
    cdef _Material pymat 
    cdef cpp_material.Material * cpp_matp
    cdef cpp_map[std.string, matp] cppmap = cpp_map[std.string, matp]()

    for key, value in pydict.items():
        pymat = value
        cpp_matp = pymat.mat_pointer
        cppmap[std.string(key)] = cpp_matp

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

