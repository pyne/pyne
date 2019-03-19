#!/usr/bin/python
"""Module for parsing FLUKA output data. FLUKA is a fully integrated particle
physics MonteCarlo simulation package. It has many applications in high
energy experimental physics and engineering, shielding, detector and telescope
design, cosmic ray studies, dosimetry, medical physics and radio-biology.
Further information on FLUKA can be obtained from
http://www.fluka.org/fluka.php

Currently, only usrbin output files can be read.

If PyMOAB is not installed, then Usrbin and UsrbinTally will not be
available to use.

"""

from warnings import warn
from pyne.utils import QAWarning

# Mesh specific imports
from pyne.mesh import Mesh, StatMesh, MeshError, HAVE_PYMOAB

if HAVE_PYMOAB:
    from pyne.mesh import NativeMeshTag
else:
    warn("The PyMOAB optional dependency could not be imported. "
         "Some aspects of the fluka module may be incomplete.",
         QAWarning)


class Usrbin(object):
    """This class is the wrapper class for UsrbinTally. This class stores
    all information for a single file that contains one or more usrbin
    tallies. The "tally" attribute provides key/value access to individual
    UsrbinTally objects.

    Attributes
    ----------
    filename : string
        Path to Fluka usrbin file
    tally : dict
        A dictionary with user-specified tally names as keys and UsrbinTally
        objects as values.
    """

    def __init__(self, filename):
        """Parameters
        ----------
        filename : string
            FLUKA USRBIN file
        """

        if not HAVE_PYMOAB:
            raise RuntimeError("PyMOAB is not available, "
                               "unable to create Meshtal.")

        self.tally = {}

        with open(filename, 'r') as fh:
            self._read_tallies(fh)

    def _read_tallies(self, fh):
        """Read in all of the USRBIN tallies from the USRBIN file.
        """
        line = fh.readline()

        while (line != "" and line[0] == '1'):
            new_tally = UsrbinTally(fh)
            self.tally[new_tally.name] = new_tally
            line = fh.readline()


class UsrbinTally(Mesh):
    """This class reads a single FLUKA USRBIN tally from a USRBIN file.

    Attributes
    ----------
    coord_sys : string
        The coordinate system used. Either "Cartesian", "R-Z", "R-Phi-Z", or
        user-defined. Only "Cartesian" is supported.
    name : string
        The user-defined name for the tally
    particle : string
        The number code corresponding to the particle tracked in tally.
        For complete list visit http://www.fluka.org/fluka.php?id=man_onl&sub=7
    x_bounds : list of floats
        The locations of mesh vertices in the x direction
    y_bounds : list of floats
        The locations of mesh vertices in the y direction
    z_bounds : list of floats
        The locations of mesh vertices in the z direction
    part_data_tag : string
        The name of the tag for the track-length tally data.
        Follows form "part_data_X" where X is the number of the particle
    error_data_tag : string
        The name of the tag for the error data.
        Follows form "error_data_X" where X is the number of the particle
    """

    def __init__(self, fh):
        """Creates a UsrbinTally object by reading through the file

        Parameters
        ----------
        fh : filehandle
            An open usrbin file
        """

        if not HAVE_PYMOAB:
            raise RuntimeError("PyMOAB is not available, "
                               "unable to create Meshtal.")

        part_data = []
        error_data = []

        line = fh.readline()

        # Read the header for the tally.
        # Information obtained: coordinate system used, user-defined tally
        # name, particle, and x, y, and z dimension information.
        [self.coord_sys, self.name, self.particle] = line.split('"')
        self.name = self.name.strip()
        self.coord_sys = self.coord_sys.split()[0]
        self.particle = self.particle.split()[-1]

        if self.coord_sys != 'Cartesian':
            raise ValueError(
                "Only cartesian coordinate system currently supported")

        [x_info, y_info, z_info] = self._read_usrbin_head(fh)

        # Advance to start of tally data skipping blank and/or text lines.
        line = fh.readline()
        line = fh.readline()
        if "accurate deposition" in line:
            line = fh.readline()
        if "track-length binning" in line:
            line = fh.readline()

        # Read the track-length binning data (part_data) and percentage error
        # data (error_data).
        num_volume_element = x_info[2]*y_info[2]*z_info[2]
        part_data += [float(x) for x in line.split()]
        while (len(part_data) < num_volume_element):
            line = fh.readline()
            part_data += [float(x) for x in line.split()]
        for count in range(0, 3):
            line = fh.readline()
        while (len(error_data) < num_volume_element):
            line = fh.readline()
            error_data += [float(x) for x in line.split()]

        # create mesh object
        self.x_bounds = self._generate_bounds(x_info)
        self.y_bounds = self._generate_bounds(y_info)
        self.z_bounds = self._generate_bounds(z_info)
        self._create_mesh(part_data, error_data)

    def _read_usrbin_head(self, fh):
        """Get the minimum bound, maximum bound, number of bins, and bin width
        for each of the x, y, and z dimensions contained within the header.
        """
        line = fh.readline()
        # assume next line is x coord info
        x_info = self._parse_dimensions(line)
        line = fh.readline()
        # assume next line is y coord info
        y_info = self._parse_dimensions(line)
        line = fh.readline()
        # assume next line is z coord info
        z_info = self._parse_dimensions(line)

        line = fh.readline()

        # return lists of info for each dimension:
        # [min, max, number of bins, width]
        return x_info, y_info, z_info

    def _parse_dimensions(self, line):
        """This retrieves the specific dimensions and binning information for
        the x, y, and z dimensions. Information retrieved is the minimum and
        maximum value for each dimension, the number of bins in each direction,
        and the width of each evenly spaced bin.
        """
        tokens = line.split()
        return float(tokens[3]), float(tokens[5]), int(tokens[7]), \
            float(tokens[10])

    def _generate_bounds(self, dim_info):
        """This takes in the dimension information (min, max, bins, and width)
        and returns a list of bound values for that given dimension.
        """
        [dim_min, dim_max, bins, width] = dim_info
        bound_data = []
        for i in range(0, bins + 1):
            bound_data.append(dim_min+(i*width))
        return bound_data

    def _create_mesh(self, part_data, error_data):
        """This will create the mesh object with the name of the tally
        specified by the user. One mesh object contains both the part_data and
        the error_data.
        """
        super(UsrbinTally, self).__init__(structured_coords=[self.x_bounds,
                                                             self.y_bounds, self.z_bounds],
                                          structured=True,
                                          structured_ordering='zyx',
                                          mats=None)
        self.part_data_tag = NativeMeshTag(size=1, dtype=float, mesh=self,
                                           name="part_data_{0}".format(self.particle))
        self.error_data_tag = NativeMeshTag(size=1, dtype=float, mesh=self,
                                            name="error_data_{0}".format(self.particle))
        self.part_data_tag[:] = part_data
        self.error_data_tag[:] = error_data
