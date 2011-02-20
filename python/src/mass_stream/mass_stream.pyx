"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
#from cython.operator cimport reference as ref
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std
cimport cpp_mass_stream
cimport stlconverters as conv

import isoname
import os

cdef class MassStream:
    """MassStream fuel cycle flow object.

    Args:
        * compdict (dict or str): This is the input isotopic component dictionary.
          This dictionary need not be normalized; MassStream initialization will
          automatically renormalize the stream.  Thus compdict simply is a dictionary
          of relative weights.  The keys of compdict must be integers representing
          isotopes in zzaaam <isoform> form.  The values are floats for 
          each isotope's weight fraction.

          If a string is provided instead of a dictionary, then MassStream will
          read in the compdict vector from a file at the string's location using 
          the load_from_text() method.  

          If no compdict is provided, an empty MassStream object is constructed.

    Keyword Args:
        * mass (float): This is the mass of the new stream. If the mass provided
          is negative (default -1.0) then the mass of the new stream is calculated from 
          the sum of compdict's components before normalization.  If the mass here
          is positive or zero, then this mass overrides the calculated one.
        * name (str):  A string label for the stream.  Helpful for large numbers of 
          streams. Default ''.
    """

    #cdef cpp_mass_stream.MassStream * ms_pointer

    def __cinit__(self, isovec=None, float mass=-1.0, char * name=''):
        """MassStream C++ constuctor."""
        cdef cpp_map[int, double] comp_dict

        if isinstance(isovec, dict):
            # Mass Stream from dict
            comp_dict = conv.dict_to_map_int_dbl(isovec)
            self.ms_pointer = new cpp_mass_stream.MassStream(comp_dict, mass, std.string(name))

        elif isinstance(isovec, basestring):
            # Mass Stream from file
            if os.path.exists(isovec):
                print "Here"
                self.ms_pointer = new cpp_mass_stream.MassStream(<char *> isovec, mass, std.string(name))
            else:
                raise IOError("The file {0} does not exist!".format(isovec))

        elif isovec is None:
            # Empty mass stream
            self.ms_pointer = new cpp_mass_stream.MassStream()

        else:
            # Bad MassStream 
            raise TypeError("The mass stream isotopic vector must be a dict, str, or None.")


    def __dealloc__(self):
        """MassStream C++ destructor."""
        del self.ms_pointer


    #
    # Class Attributes
    #

    property comp:
        def __get__(self):
            comp_dict = conv.map_to_dict_int_dbl(self.ms_pointer.comp)
            return comp_dict

        def __set__(self, dict comp):
            self.ms_pointer.comp = conv.dict_to_map_int_dbl(comp)


    property mass:
        def __get__(self):
            return self.ms_pointer.mass

        def __set__(self, double mass):
            self.ms_pointer.mass = mass


    property name:
        def __get__(self):
            cdef std.string ms_name = self.ms_pointer.name
            return ms_name.c_str()

        def __set__(self, char * name):
            self.ms_pointer.name = std.string(name)

    #
    # Class Methods
    #

    def norm_comp_dict(self):
        """Normalizes the comp(osition), preserving the mass of the isotopic vector as mass."""
        self.ms_pointer.norm_comp_dict()


    def load_from_hdf5(self, char * filename, char * groupname, int row=-1):
        """A MassStream object may be initialized from an HDF5 file.
        The HDF5 representation of a MassStream is a group that holds several 
        extendable array datasets.  One array is entitled "Mass" while the other datasets
        are nuclide names in LLAAAM form ("U235", "NP237", *etc*).  For example::

            File.h5 (file)
                |-- MassStream (group)
                    |-- Mass (array)
                    |-- H1 (array)
                    |-- O16 (array)
                    |-- U235 (array)
                    |-- PU239 (array)
                    |-- ...

        The arrays are all of length N, where each row typically represents a different 
        fuel cycle pass.  The sum of all of the nuclide arrays should sum to one, like 
        MassStream.comp. 

        Args:
            * filename  (str): Path to HDF5 file that contains the data to read in.    
            * groupname (str): Path to HDF5 group that represents the data. 
              In the above example, groupname = "/MassStream".    

        Keyword Args:
            * row (int): The index of the arrays from which to read the data.  This 
              ranges from 0 to N-1.  Defaults to the last element of the array.
              Negative indexing is allowed (row[-N] = row[0]).

        Usage:
            This function loads data into a pre-existing :class:`MassStream`.  
            Initialization is therefore a two-step process::

                ms = MassStream()
                ms.load_from_hdf5("afile.h5", "/foo/bar/ms", -3)
        """
        self.ms_pointer.load_from_hdf5(std.string(filename), std.string(groupname), row)


    def load_from_text(self, char * filename):
        """A MassStream object may be initialized from a simple text file.
        The text representation of MassStreams are nuclide identifiers in the 
        first column and mass or weight values in the second column.  For example, 
        for natural uranium::

            922340  0.000055
            U235    0.00720
            92238   0.992745

        Data in this file must be whitespace separated.  Any valid nuclide naming
        scheme may be used for any isotope.

        Args:
            * filename (str): Path to HDF5 file that contains the data to read in.    

        Usage:
            This function loads data into a pre-existing MassStream.  
            Initialization is therefore a two-step process::

            ms = MassStream()
            ms.load_from_text("natu.h5")

        This function is most often called implicitly by the MassStream constructor.
        """
        self.ms_pointer.load_from_text(filename)


    def Print(self):
        """This prints a string representation of the MassStream to stdout.  Print is 
        particularly useful in C++.  In Python, this method simply duplicates 
        the functionality you would get from the built-in str() function.
        """
        self.ms_pointer.Print()


    def Normalize(self):
        """This convenience function normalizes the mass stream by setting its mass = 1.0."""
        self.ms_pointer.Normalize()


    def multByMass(self):
        """This function multiplies comp by mass and returns the resultant isotopic vector.

        Returns:
            * isovec(dict): For a MassStream ms, 

              .. math:: \mbox{isovec[iso]} = \mbox{ms.comp[iso]} \times \mbox{ms.mass}
        """
        cdef cpp_map[int, double] cpp_isovec = self.ms_pointer.multByMass()
        isovec = conv.map_to_dict_int_dbl(cpp_isovec)
        return isovec


    def atomic_weight(self):
        """This method returns the atomic weight of the comp of this MassStream.  Note that this is 
        only a rough estimate since this function is not yet coupled with measured atomic weights.

        Returns:
            * atomic_weight (float): Atomic weight in [amu]."""
        return self.ms_pointer.atomic_weight()


    #
    # Substream Methods
    #

    def getSubStream(self, iso_sequence, char * name=""):
        """Grabs a subset of the mass streams and returns a new stream comprised of only
        the specified nuclides.  The elements or isotopes included in the new substream
        are determined by iso_sequence. 

        The input here is seen as a suggestion and so no error is raised if a nuclide 
        is asked for via iso_sequence that is not present in the original mass stream.

        Args:
            * isoname (sequence): Elements and isotopes to be taken from current stream.
              Members of this list must be integers.  For example, [92, 942390]
              would take all uranium atoms and Pu-239.  
            * name (str): The name of the substream.

        Returns:
            * substream (MassStream): A new mass stream object that only 
              has the members given in iso_sequence.  The mass of the substream
              is calculated based on the weight fraction composition and mass
              of the original mass stream.
        """
        # Make an isotopic set 
        cdef int iso_zz
        cdef cpp_set[int] iso_set = cpp_set[int]()
        for iso in iso_sequence:
            if isinstance(iso, int):
                if (iso in isoname.zzLL):
                    iso_zz = iso
                else:
                    iso_zz = isoname.mixed_2_zzaaam(iso)

            elif isinstance(iso, basestring):
                iso_str = iso.upper()
                if (iso_str in isoname.LLzz):
                    iso_zz = isoname.LLzz[iso_str]
                else:
                    iso_zz = isoname.mixed_2_zzaaam(iso)

            else:
                raise TypeError("Isotopes must be strings or integers.")

            iso_set.insert(iso_zz)

        # Make new python version of this mass stream
        py_ms = MassStream()
        py_ms.ms_pointer[0] = self.ms_pointer.getSubStream(iso_set, std.string(name))
        return py_ms


    def getU(self, char * name=""):
        """Convenience method that gets the Uranium portion of a mass stream.

        Args:
            * name (str): The name of the substream.

        Returns:
            * substream (MassStream): A new mass stream object that only 
              has Uranium members. 
        """
        py_ms = MassStream()
        py_ms.ms_pointer[0] = self.ms_pointer.getU(std.string(name))
        return py_ms
        

    def getPU(self, char * name=""):
        """Convenience method that gets the Plutonium portion of a mass stream.

        Args:
            * name (str): The name of the substream.

        Returns:
            * substream (MassStream): A new mass stream object that only 
              has Plutonium members. 
        """
        py_ms = MassStream()
        py_ms.ms_pointer[0] = self.ms_pointer.getPU(std.string(name))
        return py_ms
        

    def getLAN(self, char * name=""):
        """Convenience method that gets the Lanthanide portion of a mass stream.

        Args:
            * name (str): The name of the substream.

        Returns:
            * substream (MassStream): A new mass stream object that only 
              has Lanthanide members. 
        """
        py_ms = MassStream()
        py_ms.ms_pointer[0] = self.ms_pointer.getLAN(std.string(name))
        return py_ms
        

    def getACT(self, char * name=""):
        """Convenience method that gets the Actinide portion of a mass stream.

        Args:
            * name (str): The name of the substream.

        Returns:
            * substream (MassStream): A new mass stream object that only 
              has Actinide members. 
        """
        py_ms = MassStream()
        py_ms.ms_pointer[0] = self.ms_pointer.getACT(std.string(name))
        return py_ms
        

    def getTRU(self, char * name=""):
        """Convenience method that gets the Transuranic portion of a mass stream.

        Args:
            * name (str): The name of the substream.

        Returns:
            * substream (MassStream): A new mass stream object that only 
              has Transuranic members. 
        """
        py_ms = MassStream()
        py_ms.ms_pointer[0] = self.ms_pointer.getTRU(std.string(name))
        return py_ms
        

    def getMA(self, char * name=""):
        """Convenience method that gets the Minor Actinide portion of a mass stream.

        Args:
            * name (str): The name of the substream.

        Returns:
            * substream (MassStream): A new mass stream object that only 
              has Minor Actinide members. 
        """
        py_ms = MassStream()
        py_ms.ms_pointer[0] = self.ms_pointer.getMA(std.string(name))
        return py_ms
        

    def getFP(self, char * name=""):
        """Convenience method that gets the Fission Product portion of a mass stream.

        Args:
            * name (str): The name of the substream.

        Returns:
            * substream (MassStream): A new mass stream object that only 
              has Fission Product members. 
        """
        py_ms = MassStream()
        py_ms.ms_pointer[0] = self.ms_pointer.getFP(std.string(name))
        return py_ms
        

    #
    # Operator Overloads
    #

    # Addition

    def __add_float__(MassStream x, double y):
        py_ms = MassStream()
        py_ms.ms_pointer[0] = x.ms_pointer[0] + y
        return py_ms
        

    def __add_mass_stream__(MassStream x, MassStream y):
        py_ms = MassStream()
        py_ms.ms_pointer[0] = x.ms_pointer[0] + y.ms_pointer[0]
        return py_ms


    def __add__(x, y): 
        if isinstance(x, MassStream) and isinstance(y, MassStream):
            return x.__add_mass_stream__(y)
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

    def __mul_float__(MassStream x, double y):
        py_ms = MassStream()
        py_ms.ms_pointer[0] = x.ms_pointer[0] * y
        return py_ms


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

    def __div_float__(MassStream self, double y):
        py_ms = MassStream()
        py_ms.ms_pointer[0] = self.ms_pointer[0] / y
        return py_ms


    def __div__(MassStream self, y):
        if isinstance(y, float):
            return self.__div_float__(y)
        elif isinstance(y, int):
            return self.__div_float__(float(y))
        else:
            return NotImplemented


    def __rdiv__(MassStream self, y):
        return self.__div__(y)

    
    def __truediv__(MassStream self, y):
        return self.__div__(y)




##############################
### Mass Stream Converters ###
##############################

# <string, MassStream *>

cdef cpp_map[std.string, msp] dict_to_map_str_msp(dict pydict):
    cdef MassStream pyms 
    cdef cpp_mass_stream.MassStream * cpp_msp
    cdef cpp_map[std.string, msp] cppmap = cpp_map[std.string, msp]()

    for key, value in pydict.items():
        pyms = value
        cpp_msp = pyms.ms_pointer
        cppmap[std.string(key)] = cpp_msp

    return cppmap


cdef dict map_to_dict_str_msp(cpp_map[std.string, msp] cppmap):
    pydict = {}
    cdef MassStream pyms 
    cdef cpp_map[std.string, msp].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pyms = MassStream()
        pyms.ms_pointer[0] = deref(deref(mapiter).second)
        pydict[deref(mapiter).first.c_str()] = pyms
        inc(mapiter)

    return pydict

