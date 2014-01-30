#!/usr/bin/python

#from itaps import iMesh
#from pyne.mesh import Mesh, StatMesh

import math
#from pyne.mcnp import MeshTally
#from pyne.mesh import Mesh, StatMesh, MeshError
class USR():
    """Wrapper class for Usrbin() class. This class stores information for the
       complete file given.
    """
    def __init__(self):
#        self.track_tally = {}
#        with open(filename) as fh:
#            self.read_number_of_bins(fh)
        pass

    def read_number_of_bins(self, fh):
        """reads how many track-length binning and their associated errors are in one file
        """ 
        track_tally = 0
        for line in fh:
            if "Cart" in line:
               track_tally = track_tally + 1
               self.track_tally[track_tally] = Usrbin(f, track_tally)
        print track_tally
        return int(track_tally)

class Usrbin():
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
        # print my_mesh
        return my_mesh

    def get_name_and_system(self, line): # retrieves user defined binning name, coordinate system, and particle type used
        # readlines until finds "cartesian" -- then parse as have already
        name = line.split('"')[1]
        coord_sys = line.split()[0]
        particle_type = line[len(line)-1]
        return name, coord_sys, particle_type

    def get_coordinate_system(self, line): # retrieves specific dimensions and binning information for each x, y, and z dimensions
        # automatically call after get_name_... 
        number_of_bins = line.split()[7]
        coord_min = line.split()[3]
        coord_max = line.split()[5]
        bin_width = line.split()[10]
        return int(number_of_bins), float(coord_min), float(coord_max), float(bin_width)

    def matrix_organization(self, line): # retrieves information about how the matrix of information is organized for each binning
        matrix_organization = line.split()[5]
        # for a cartesian coordinate system:
        i = matrix_organization[3] # changes last
        j = matrix_organization[6] # changes second
        k = matrix_organization[9] # changes first
        column_format = line.split()[7]
        column_number = column_format.split(',')[2].split('(')[0] # number of columns in table
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
         
        name = self.get_name_and_system(line)[0]
        coord_sys = self.get_name_and_system(line)[1]
        print coord_sys
        particle = self.get_name_and_system(line)[2]
        if coord_sys == "Cartesian":
            fh.readline()
            # assume next line is X coord info
            x_bins = self.get_coordinate_system(line)[0]
            x_min = self.get_coordinate_system(line)[1]
            x_max = self.get_coordinate_system(line)[2]
            x_width = self.get_coordinate_system(line)[3]
#            x_info = [x_bins, x_min, x_max, x_width]
            fh.readline()
            # assume next line is y coord info
            y_bins = self.get_coordinate_system(line)[0]
            y_min = self.get_coordinate_system(line)[1]
            y_max = self.get_coordinate_system(line)[2]
            y_width = self.get_coordinate_system(line)[3]
#            y_info = {'bins': y_bins, 'min': y_min, 'max': y_max, 'width': y_width}
            fh.readline()
            # assume next line is z coord info
            z_bins = self.get_coordinate_system(line)[0]
            z_min = self.get_coordinate_system(line)[1]
            z_max = self.get_coordinate_system(line)[2]
            z_width = self.get_coordinate_system(line)[3]
#            z_info = {'bins': z_bins, 'min': z_min, 'max': z_max, 'width': z_width}
        else:
            print "Coordinate sytem is not Cartesian"
        fh.readline()
        # collect how data is arranged (number of columns in matrix)
        columns = self.matrix_organization(line)[2]
#        print x_info, y_info, z_info, columns
        return x_bins, x_min, x_max, x_width, y_bins, y_min, y_max, y_width, z_bins, z_min, z_max, z_width, columns


    def read_usrbin(self, filename): # combines all above functions to place data in a list
        fh = open(filename)
        mesh_tally = 0
        cart_tf = False
#        x_bin_tf = False
#        y_bin_tf = False
#        z_bin_tf = False
        track_length_tf = False
        track_error_tf = False
        for line in fh:
            part_data = []
            error_data = []
            if (cart_tf and track_length_tf):
#               x_bin_tf and y_bin_tf and z_bin_tf): # place track length data into list
               data = self.read_data(line)
               part_data.append(data)
#               mesh = self.create_mesh(x_bounds,y_bounds,z_bounds)
#               print part_data
#               print mesh
            if (cart_tf  and not track_length_tf and track_error_tf):
#               x_bin_tf and y_bin_tf and z_bin_tf): # place error data into list
               data = self.read_data(line)
               error_data.append(data)
#               print error_data
            if ("Cartesian" in line): # looks for start of binning and retrieves appropriate information
               x_bins = self.read_header(fh, line)[0]
               x_min = self.read_header(fh, line)[1]
               x_max = self.read_header(fh, line)[2]
               x_width = self.read_header(fh, line)[3] 
               y_bins = self.read_header(fh, line)[4]
               y_min = self.read_header(fh, line)[5]
               y_max = self.read_header(fh, line)[6]
               y_width = self.read_header(fh, line)[7]
               z_bins = self.read_header(fh, line)[8]
               z_min = self.read_header(fh, line)[9]
               z_max = self.read_header(fh, line)[10]
               z_width = self.read_header(fh, line)[11]
               columns = self.read_header(fh, line)[12]
#               assigned_name = self.get_name_and_system(line)[0]
#               coordinate_system = self.get_name_and_system(line)[1]
#               particle = self.get_name_and_system(line)[2]
#               print assigned_name, coordinate_system, particle
               mesh_tally = mesh_tally + 1
               print mesh_tally
               print x_bins, x_min, x_max, x_width, y_bins, y_min, y_max, y_width, z_bins, z_min, z_max, z_width, columns
               cart_tf = True
#            if ("X coordinate" in line): # retrieves all information about the x dimension
#               print x_bins, x_min, x_max, x_width
#              x_bounds=self.generate_bounds(x_min,x_max,x_width,x_bins)
#               print x_bounds
#               x_bin_tf = True
#            if ("Y coordinate" in line): # retrieves all information about the y dimension
#               print y_bins, y_min, y_max, y_width
#               y_bounds=self.generate_bounds(y_min,y_max,y_width,y_bins)
#               print y_bounds
#               y_bin_tf = True
#            if ("Z coordinate" in line): # retrieves all information about the z dimension
#               print z_bins, z_min, z_max, z_width
#               z_bounds=self.generate_bounds(z_min,z_max,z_width,z_bins)
#               print z_bounds
#               z_bin_tf = True
#            if ("Data" in line): # retrieves information about organization of data in table
#               col_n = self.matrix_organization(line)[3]
#               print col_n
            if ("this is a track-length binning" in line): # initializes the collection of track length data into list
               track_length_tf = True
            if ("Percentage errors" in line): # stops collection of track length data and initializes collection of error data
               track_length_tf = False
               track_error_tf = True
                      
        # the following will the total number of data points and the total number of rows in file of data for each binning
#        total_bins = x_bins*y_bins*z_bins
#        total_rows = int(math.ceil(float(total_bins/col_n)))
        # mh=open("new_mesh.h5m")
        
#        print total_bins, total_rows, my_mesh

#my_file = USR()
#my_file.read_number_of_bins("fng_dose_usrbin_25.lis")
my_file = Usrbin()
my_file.read_usrbin("fng_dose_usrbin_25.lis")
