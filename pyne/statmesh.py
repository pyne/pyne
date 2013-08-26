import math

from itaps import iMesh
from pyne.scdmesh import ScdMesh

class StatMesh(ScdMesh):
    """This is a class to store and manipulate information on a MOAB mesh
    containing statistical quantities (with an average and statistical error).
    """

    def __init__(self, mesh):
        self.mesh = mesh

    def add(self, mesh_2, tag_name, out_mesh):
        """For all volumes elements add the <tag_name> values from <mesh_2> to
        the values on self.
        """

        _verify_valid(self.mesh, mesh_2, tag_name):

        instance_mesh_iterator = self.mesh.iterHex("xyz")
        argument_mesh_iterator = mesh_2.iterHex("xyz")
        output_mesh_iterator = out_mesh.iterHex("xyz")

        instance_mesh_error_iterator = self.mesh.iterHex("xyz")
        argument_mesh_error_iterator = mesh_2.iterHex("xyz")
        output_mesh_error_iterator = out_mesh.iterHex("xyz")

        instance_tag = self.mesh.imesh.getTagHandle(tag_name)
        argument_tag = mesh2.imesh.getTagHandle(tag_name)
        output_tag = self.out_mesh.imesh.getTagHandle(tag_name)

        instance_error_tag = self.mesh.imesh.getTagHandle(tag_name + "error")
        argument_error_tag = mesh2.imesh.getTagHandle(tag_name + "error")
        output_error_tag = self.out_mesh.imesh.getTagHandle(tag_name + "error")

        for ve_1, ve_2, ve_out, ve_1_error, ve_2_error, ve_out_error in 
                zip(instance_mesh_iterator, argument_mesh_iterator, output_mesh_iterator):

            output_tag[ve_out] = instance_tag[ve_1] + argument_tag[ve_2]
            output_error_tag[ve_out_error] = math.sqrt(math.exp(instance_tag[ve_1], 2) + math.exp(argument_tag[ve_2], 2))       
   

    def _verify_valid(self, mesh_2, tag_name, out_mesh)

       #Ensure mesh division are the same
       dims = ["x", "y", "z"]
       for dim in dims:
           if self.mesh.getDivisions(dim) != mesh_2.getDivisions(dim):
               raise ValueError("Meshes are not the same size.")

        #Check that tag exists on the instance mesh.
        try:
            self.mesh.imesh.getTagHandle(tag_name)   
        except iBase.TagNotFoundError as e:
            print "Tag not found on instance mesh"
            raise e

        #Check that error tag exists on the instance mesh.
        try:
            self.mesh.imesh.getTagHandle(tag_name + "_error")
        except iBase.TagNotFoundError as e:
            print "Error tag not found on instance mesh"
            raise e

        #Check that tag exists on the argument mesh.
        try:
            mesh2.imesh.getTagHandle(tag_name)
        except iBase.TagNotFoundError as e:
            print "Tag not found on argument mesh"
            raise e

        #Check that error tag exists on the argument mesh.
        try:
            mesh2.imesh.getTagHandle(tag_name + "_error")
        except iBase.TagNotFoundError as e:
            print "Error tag not found on argument mesh"
            raise e

        #If tag does not already exist on output mesh, create it
        try:
            out_mesh.imesh.getTagHandle(tag_name)
        except iBase.TagNotFoundError:
            out_mesh.imesh.createTag(tag_name)

        #If error tag does not already exist on output mesh, create it
        try:
            out_mesh.imesh.getTagHandle(tag_name + "_error")
        except iBase.TagNotFoundError:
            out_mesh.imesh.createTag(tag_name + "_error")

     
    def subtract(self, mesh_2, out_mesh):


    def multiply(self, mesh_2, out_mesh):


    def divide(self, mesh_2, out_mesh):


    def minimum(tag_name):

    def maximum(tag_name):

    def scalar_multiply():

    def scalar_divide(): 
