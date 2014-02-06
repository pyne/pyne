#!/usr/bin/python

#from itaps import iMesh
#from pyne.mesh import Mesh, StatMesh

import math
#from pyne.mcnp import MeshTally
#from pyne.mesh import Mesh, StatMesh, MeshError
class Meshtal():
    """Wrapper class for Usrbin() class. This class stores information for the
       complete file given.
    """
    def __init__(self):
        self.track_tally = {}
        with open(filename) as fh:
            self.read_tally_sets(fh)
        pass

    def read_tally_sets(self, fh):
        """finds next start of new set of track-length binning data and indexes it
        """ 
        track_tally = 1
#        for line in fh:
#            if "Cart" in line:
#               track_tally = track_tally + 1
#               self.track_tally[track_tally] = Usrbin(f, track_tally)
#        print track_tally
	line = fh.readline()
	while ("Cart" not in line):
	    tally_number = track_tally
	    self.tally[tally_number] = Meshtally(fh, tally_number)

	track_tally = track_tally + 1

class Meshtally():
    """This class stores all the information from Fluka USRBIN output
    file with a single **or multiple** track length binnings and their 
    associated error binnings. Currently this class only supports a 
    cartesian coordinate system.

    Attributes
    ----------
    assigned_name : str
        the name given by user to the particular binning
    coord_max : float
        maximum value of one dimension
    coord_min : float
        minimum value of one dimension
    coordinate_system : str
        type of coordinate system being used
    data : list of lists
        list of track length data for a particular binning
    error_data : list of lists
        list of error data for each binning
    mesh_tally : int
        tallies the number of track length binnings in file
    number_of_bins : int
        the number of evenly spaced bins corresponding to one dimension
        (x, y, or z)
    particle : 
        particle type given by its corresponding numerical identifier
    total_bins :
        total number of bins in the three dimensional system
    Notes
    -----
    Attribute names are identical to names speficied in a USRBIN file
    description in the Fluka Manual, section 7.77 USRBIN input commands.

    """
    def __init__(self):
#(self, fh, track_tally, line):
#        self.track_tally = track_tally
#        self.get_name_and_system(line)
#        self.get_coordinate_system(line)
#        self.matrix_organization(line)
        pass

    def create_mesh(self,x_bounds,y_bounds,z_bounds):
        """Instantiate a Mesh object and tag the iMesh instance
        with results and relative errors.
        """    
        my_mesh = structured_coords=[self.x_bounds,self.y_bounds, self.z_bounds],structured=True
        return my_mesh

    def get_name_and_system(self, line): # retrieves user defined binning name, coordinate system, and particle type used
        # readlines until finds "cartesian" -- then parse as have already
        name = line.split('"')[1]
        coord_sys = line.split()[0]
        line_split = line.split()
        particle_type = line_split[len(line_split)-1]
        return name, coord_sys, particle_type

    def get_coordinate_system(self, line): # retrieves specific dimensions and binning information for each x, y, and z dimensions
        # automatically call after get_name_...
#        print line 
        number_of_bins = line.split()[7]
        coord_min = line.split()[3]
        coord_max = line.split()[5]
        bin_width = line.split()[10]
#        print number_of_bins
        return int(number_of_bins), float(coord_min), float(coord_max), float(bin_width)

    def matrix_organization(self, line): # retrieves information about how the matrix of information is organized for each binning
        matrix_organization = line.split()[5]
        # for a cartesian coordinate system:
        i = matrix_organization[3] # changes last
        j = matrix_organization[6] # changes second
        k = matrix_organization[9] # changes first
        column_format = line.split()[7]
        column_number = column_format.split(',')[2].split('(')[0] # number of columns in table
        print column_format, column_number
        return i, j, k, int(column_number)

    def read_data(self, line): # reads data
        """
        reads data from string, splits into space delimited chunks and returns a list of items
        """
        data_line = line.split()
        return data_line

    def generate_bounds(self,dir_min,dir_max,bin_width,bounds): 
        """ takes arguments of the the minimum and max values in a given direction and the number
        of splits, returns list of boundaries
        """
        bound_data=[]
        bound_data.append(dir_min) #append the first boundary
        for i in range(1,bounds+1):
            bound_data.append(dir_min+(i*bin_width))
        return bound_data

    def read_header(self, fh, line):
        # call get_name_and_system
        # if Cartesian.. then: (assume that x, y and z are next lines)
        #   x_bounds = 'get coord syst' (for x, y and z)
        # call matrix_org (Elliott's read_column_ordr)
        
        print line 
        name = self.get_name_and_system(line)[0]
        coord_sys = self.get_name_and_system(line)[1]
        particle = self.get_name_and_system(line)[2]
        if coord_sys == "Cartesian":
            line = fh.readline()
            print line
            # assume next line is X coord info
            x_bins = self.get_coordinate_system(line)[0]
            x_min = self.get_coordinate_system(line)[1]
            x_max = self.get_coordinate_system(line)[2]
            x_width = self.get_coordinate_system(line)[3]
            print x_bins, x_min, x_max, x_width
            x_info = [x_bins, x_min, x_max, x_width]
            line = fh.readline()
            # assume next line is y coord info
            y_bins = self.get_coordinate_system(line)[0]
            y_min = self.get_coordinate_system(line)[1]
            y_max = self.get_coordinate_system(line)[2]
            y_width = self.get_coordinate_system(line)[3]
            y_info = [y_bins, y_min, y_max, y_width]
            line = fh.readline()
            # assume next line is z coord info
            z_bins = self.get_coordinate_system(line)[0]
            z_min = self.get_coordinate_system(line)[1]
            z_max = self.get_coordinate_system(line)[2]
            z_width = self.get_coordinate_system(line)[3]
            z_info = [z_bins, z_min, z_max, z_width]
        else:
            print "Coordinate sytem is not Cartesian"
        line = fh.readline()
        # collect how data is arranged (number of columns in matrix)
        columns = self.matrix_organization(line)[3]
#        print x_info, y_info, z_info, columns
        return x_info, y_info, z_info, columns


    def read_usrbin(self, filename): # combines all above functions to place data in a list
        fh = open(filename)
        mesh_tally = 0
        cart_tf = False
        track_length_tf = False
        track_error_tf = False
        part_data = []
        error_data = []
        while (not cart_tf and not track_length_tf):
            line = fh.readline()
            if ("Cartesian" in line):
               [x_info, y_info, z_info, columns] = self.read_header(fh, line)
               print columns
               mesh_tally = mesh_tally + 1
               cart_tf = True
               track_error_tf =  False

        while (cart_tf and not track_length_tf):
            line = fh.readline()
            if ("this is a track-length binning" in line): # initializes the collection of track length data into list
                track_length_tf = True

        while (cart_tf and track_length_tf):
            line = fh.readline()
            data = self.read_data(line)
            part_data.append(data)
            if ("Percentage errors" in line): 
            # stops collection of track length data and initializes collection of error data
                track_length_tf = False
                track_error_tf = True
                cart_tf = False

        while (not cart_tf and track_error_tf and not track_length_tf):
            line = fh.readline()
            data = self.read_data(line)
            error_data.append(data)
            if ("Cartesian" in line or line == ""):
                cart_tf = True
	print part_data, error_data





my_file = Meshtally()
my_file.read_usrbin("fng_dose_usrbin_23.lis")
