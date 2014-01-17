#!/usr/bin/python

#from itaps import iMesh
#from pyne.mesh import Mesh, StatMesh

import math

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
        pass

    
    def get_name_and_system(self, line): # retrieves user defined binning name, coordinate system, and particle type used
        name = line.split('"')[1]
        coord_sys = line.split()[0]
        particle_type = line[len(line)-1]
        return name, coord_sys, particle_type

    def get_coordinate_system(self, line): # retrieves specific dimensions and binning information for each x, y, and z dimensions
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
        data_line = line.split()
        return data_line

    def read_usrbin(self, filename): # combines all above functions to place data in a list
        fh = open(filename)
        mesh_tally = 0
        cart_tf = False
        x_bin_tf = False
        y_bin_tf = False
        z_bin_tf = False
        track_length_tf = False
        track_error_tf = False
        for line in fh:
            part_data = []
            error_data = []
            if (cart_tf and x_bin_tf and y_bin_tf and z_bin_tf and track_length_tf): # place track length data into list
               data = self.read_data(line)
               part_data.append(data)
#               print part_data
            if (cart_tf  and x_bin_tf and y_bin_tf and z_bin_tf and not track_length_tf and track_error_tf): # place error data into list
               data = self.read_data(line)
               error_data.append(data)
#               print error_data
            if ("Cartesian" in line): # looks for start of binning and retrieves appropriate information
               assigned_name = self.get_name_and_system(line)[0]
               coordinate_system = self.get_name_and_system(line)[1]
               particle = self.get_name_and_system(line)[2]
               print assigned_name, coordinate_system, particle
               mesh_tally = mesh_tally + 1
               print mesh_tally
               cart_tf = True
            if ("X coordinate" in line): # retrieves all information about the x dimension
               x_bins = self.get_coordinate_system(line)[0]
               x_min = self.get_coordinate_system(line)[1]
               x_max = self.get_coordinate_system(line)[2]
               x_width = self.get_coordinate_system(line)[3]
               print x_bins, x_min, x_max, x_width
               x_bin_tf = True
            if ("Y coordinate" in line): # retrieves all information about the y dimension
               y_bins = self.get_coordinate_system(line)[0]
               y_min = self.get_coordinate_system(line)[1]
               y_max = self.get_coordinate_system(line)[2]
               y_width = self.get_coordinate_system(line)[3]
               print y_bins, y_min, y_max, y_width
               y_bin_tf = True
            if ("Z coordinate" in line): # retrieves all information about the z dimension
               z_bins = self.get_coordinate_system(line)[0]
               z_min = self.get_coordinate_system(line)[1]
               z_max = self.get_coordinate_system(line)[2]
               z_width = self.get_coordinate_system(line)[3]
               print z_bins, z_min, z_max, z_width
               z_bin_tf = True
            if ("Data" in line): # retrieves information about organization of data in table
               col_n = self.matrix_organization(line)[3]
               print col_n
            if ("this is a track-length binning" in line): # initializes the collection of track length data into list
               track_length_tf = True
            if ("Percentage errors" in line): # stops collection of track length data and initializes collection of error data
               track_length_tf = False
               track_error_tf = True
        # the following will the total number of data points and the total number of rows in file of data for each binning
        total_bins = x_bins*y_bins*z_bins
        total_rows = int(math.ceil(float(total_bins/col_n)))
        print total_bins, total_rows



my_file = Usrbin()
my_file.read_usrbin("fng_dose_usrbin_23.lis")
