#!/usr/bin/python

from itaps import iMesh
from pyne.mesh import Mesh, StatMesh, IMeshTag

import math
#from pyne.mesh import Mesh, StatMesh, MeshError, IMeshTag

class USRBIN(Mesh):
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
    def __init__(self, filename):
        fh = open(filename)
	self.read_usrbin(fh)

    def _create_mesh(self, file_handle):
        print "creating mesh"
	super(USRBIN, self).__init__(structured_coords=[self.x_bounds,
                                        self.y_bounds, self.z_bounds],
                                        structured=True,structured_ordering='zyx')

        print self.x_bounds,self.y_bounds,self.z_bounds

	print "done creating mesh"
	self.part_data_tag = IMeshTag(size=1, dtype=float, mesh=self, name="part_data_{0}".format(self.particle))
	self.error_data_tag = IMeshTag(size=1, dtype=float, mesh=self, name="error_data_{0}".format(self.particle))

        print len(self.part_data)
        print len(self.error_data)
	
	self.part_data_tag[:] = self.part_data
	self.error_data_tag[:] = self.error_data


 #       self.write_hdf5(self.name.strip()+".h5m")
#        from itaps import iMesh
        self.mesh.save(self.name.strip()+".h5m")
	

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
        for i, ve in enumerate(data_line):
	    data_line[i] = float(ve)
#        for item in data_line:
#	    item = float(item)
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
        name = self.get_name_and_system(line)[0]
        coord_sys = self.get_name_and_system(line)[1]
        particle = self.get_name_and_system(line)[2]
        if coord_sys == "Cartesian":
            line = fh.readline()
 #           print line
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
            print y_bins, y_min, y_max, y_width
            # assume next line is z coord info
            z_bins = self.get_coordinate_system(line)[0]
            z_min = self.get_coordinate_system(line)[1]
            z_max = self.get_coordinate_system(line)[2]
            z_width = self.get_coordinate_system(line)[3]
            z_info = [z_bins, z_min, z_max, z_width]
            print z_bins, z_min, z_max, z_width
        else:
            print "Coordinate sytem is not Cartesian"
        line = fh.readline()
        # collect how data is arranged (number of columns in matrix)
        columns = self.matrix_organization(line)[3]
#        print x_info, y_info, z_info, columns
        return x_info, y_info, z_info, columns


    def read_usrbin(self, file_handle): # combines all above functions to place data in a list
#        fh = open(filename)
        fh = file_handle
        self.mesh_tally = 0
        cart_tf = False
        track_length_tf = False
        track_error_tf = False
        EOF_tf = False
        self.part_data = []
        self.error_data = []
        line = True
	line = fh.readline()
        while line:
	    self.mesh_tally = self.mesh_tally + 1
	    print self.mesh_tally

#	    line = fh.readline()
	    print line
            if "1" not in line:
		print "error not a usrbin file"
                line = False
            line = fh.readline()    
#	    print line

	    self.name = self.get_name_and_system(line)[0]
	    self.particle = self.get_name_and_system(line)[2]

#	    print self.name

	    [x_info, y_info, z_info, columns] = self.read_header(fh, line)
#	    print line
	    for count in range (0,2):
		line = fh.readline()


	    if "accurate deposition" in line:
		print "!!!!!!!!!"
		print line
		print "!!!!!!!!!"
		line = fh.readline()

#	    print "!!!!!"
#	    print line
#	    print "!!!!!"
#            if "track-length binning" not in line:
#		print "not a track length binning?"
#		line = False
	    if "track-length binning" in line:
		print "!!!!!!!!!"
		print line
		print "!!!!!!!!!"
		fh.readline()
#		print line

	    # now reading track length data
            #int(math.ceil(x_info[0]*y_info[0]*z_info[0]/float(columns)))
            num_volume_element = x_info[0]*y_info[0]*z_info[0]
            while ( len(self.part_data) < num_volume_element ):
		self.part_data += [float(x) for x in line.split()]
                line = fh.readline()
	    line = fh.readline()
#	    print line
            for count in range (0,2):
		line = fh.readline()
#		print line
            while ( len(self.error_data) < num_volume_element ):
		self.error_data += [float(x) for x in line.split()]
                line = fh.readline()
                                

	    # create mesh object

            self.x_bounds = self.generate_bounds(x_info[1],x_info[2],x_info[3],x_info[0])
            self.y_bounds = self.generate_bounds(y_info[1],y_info[2],y_info[3],y_info[0])
            self.z_bounds = self.generate_bounds(z_info[1],z_info[2],z_info[3],z_info[0])

	    print "mesh #", self.mesh_tally
	    self._create_mesh(fh)

	    print "****"
	    print line
	    print "****"
#	    if "1" in line:
#		line = False
#                self.mesh_tally = self.mesh_tally + 1
#	    print line
#	    print self.mesh_tally
#	print self.part_data
   



#test1=USR("fng_dose_usrbin_22.lis")
# test2=read_tally_sets()
my_file = USRBIN("/data/prod/nasa/atic/big_run/fludag/atic_usrbin_76.lis")
#my_file.read_usrbin("fng_dose_usrbin_22.lis")
