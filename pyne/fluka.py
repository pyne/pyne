#!/usr/bin/python

#from itaps import iMesh
#from pyne.mesh import Mesh, StatMesh

import math

class Usrbin():
    """Module for parsing Fluka USRBIN output data.
    * insert description of fluka here *
    class Usrbin(Mesh):
    This Usrbin class will store all the information 
    from a Fluka USRBIN output file.

    Attributes
    ----------
    what1 : int
        type of binning selected, default = 0.0 
    what2 : int
        particle type
    what3 : int?
        logical output unit?, default = 11.0
    what4 :
    what5 :
    what6 :
    sdum  :

    Notes
    -----
    Attribute names are identical to names speficied in a USRBIN file
    description in the Fluka Manual, section 7.77 USRBIN input commands.

    """
    def __init__(self):
        pass

    def read_line1(self, line):
    	coord_sys = line[1]
        particle_type = line[len(line)] 
    
    def get_name_and_system(self, line):
        name = line.split('"')[1]
        coord_sys = line.split()[0]
        return name, coord_sys

    def get_coordinate_system(self, line):
        number_of_bins = line.split()[7]
        coord_min = line.split()[3]
        coord_max = line.split()[5]
        bin_width = line.split()[10]
        return int(number_of_bins), float(coord_min), float(coord_max), float(bin_width)

    def matrix_organization(self, line):
        matrix_organization = line.split()[5]
        # for a cartesian coordinate system:
        i = matrix_organization[3] # changes last
        j = matrix_organization[6] # changes second
        k = matrix_organization[9] # changes first
        column_format = line.split()[7]
        column_number = column_format.split(',')[2].split('(')[0]
        return i, j, k, int(column_number)

    def read_data(self, line):
        data_line = line.split()
        return data_line

    def read_usrbin(self, filename):
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
            if (cart_tf and x_bin_tf and y_bin_tf and z_bin_tf and track_length_tf):
               data = self.read_data(line)
               part_data.append(data)
               print part_data
            if (cart_tf  and x_bin_tf and y_bin_tf and z_bin_tf and not track_length_tf and track_error_tf):
               data = self.read_data(line)
               error_data.append(data)
               print error_data
            if ("Cartesian" in line):
#               print line
               assigned_name = self.get_name_and_system(line)[0]
               coordinate_system = self.get_name_and_system(line)[1]
               print assigned_name, coordinate_system
               mesh_tally = mesh_tally + 1
               print mesh_tally
               cart_tf = True
            if ("X coordinate" in line):
#               print line
               x_bins = self.get_coordinate_system(line)[0]
               x_min = self.get_coordinate_system(line)[1]
               x_max = self.get_coordinate_system(line)[2]
               x_width = self.get_coordinate_system(line)[3]
               print x_bins, x_min, x_max, x_width
               x_bin_tf = True
            if ("Y coordinate" in line):
#               print line
               y_bins = self.get_coordinate_system(line)[0]
               y_min = self.get_coordinate_system(line)[1]
               y_max = self.get_coordinate_system(line)[2]
               y_width = self.get_coordinate_system(line)[3]
               print y_bins, y_min, y_max, y_width
               y_bin_tf = True
            if ("Z coordinate" in line):
#               print line
               z_bins = self.get_coordinate_system(line)[0]
               z_min = self.get_coordinate_system(line)[1]
               z_max = self.get_coordinate_system(line)[2]
               z_width = self.get_coordinate_system(line)[3]
               print z_bins, z_min, z_max, z_width
               z_bin_tf = True
            if ("Data" in line):     
               col_n = self.matrix_organization(line)[3]
               print col_n
            if ("this is a track-length binning" in line):
               track_length_tf = True
            if ("Percentage errors" in line):
#               fh.readline()
               track_length_tf = False
               track_error_tf = True
        total_bins = x_bins*y_bins*z_bins
        total_rows = int(math.ceil(float(total_bins/col_n)))
        print total_bins, total_rows



my_file = Usrbin()
my_file.read_usrbin("fng_dose_usrbin_23.lis")
