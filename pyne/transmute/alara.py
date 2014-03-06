from pyne import alara
from pyne.mesh import Mesh
import uuid

class Transmuter(Transmuter):
    """A class for transmuting materials (possibly on meshes) using ALARA."""

    def __init__(self, element_lib, data_lib, t=0.0, phi=0.0,
                 tol=1e-7, irr_blocks=None, *args, **kwargs):
        """Parameters
        ----------
        element_lib : str
            Path to element library.
        data_library : str
            The data_library card (see ALARA user's guide).
        t : float
            Transmutations time [sec].
        phi : float
            Scalar neutron flux [n/cm^2/sec]. Required if transmuting a 
            material, not used if transmuting a mesh.
        tol : float
            Tolerance level for chain truncation.
        irr_blocks : str, optional
            Irradition-related ALARA input blocks.
        args : tuple, optional
            Other arguments ignored for compatibility with other Transmuters.
        kwargs : dict, optional
            Other keyword arguments ignored for compatibility with other 
            Transmuters.
        """
        self.element_lib = element_lib
        self.data_lib = data_lib
        self.t = t
        self.phi = phi
        self.tol = tol
        self.irr_blocks = irr_blocks

    def transmute(self, x=None, flux_tag=None, t=None, phi=None, tol=None, 
                  element_lib=None, data_lib=None, irr_blocks=None, *args, 
                  **kwargs):
        """Transmutes a material into its daughters or transmutes all materials
        in a mesh.

        Parameters
        ----------
        x : Material, PyNE Mesh object or similar
            Input material or mesh for transmutation.
        flux_tag : str, required if x is a mesh
            The name of the tag related to flux on the mesh, required if 
            transmuting a mesh, not used if transmuting a material.
        t : float, optional
            Transmutations time [sec].
        phi : float or array of floats, optional
            Neutron flux vector [n/cm^2/sec].  Currently this must either be 
            a scalar or match the group structure of EAF.
        tol : float, optional
            Tolerance level for chain truncation.
        element_lib : str, optional
            Path to element library.
        data_library : str, optional
            The data_library card (see ALARA user's guide).
        irr_blocks : str, optional
            Irradition-related ALARA input blocks.
        args : tuple, optional
            Other arguments ignored for compatibility with other Transmuters.
        kwargs : dict, optional
            Other keyword arguments ignored for compatibility with other 
            Transmuters.

        Returns
        -------
        y : Material if x is a Material
            The output material post-transmutation.
        """
        if isinstance(x, Mesh): 
            self._transmute_mesh(mesh, flux_tag, t, tol, element_lib, data_lib, 
                                 irr_blocks, args, kwargs)
        else:
            mesh = Mesh(structured_coords=[[0,1], [0,1], [0,1]], structured=True, 
                        structured_ordering='zyx', mats={0: x})
            flux_tag = "flux"
            tag = mesh.mesh.createTag(flux_tag, 1, float)
            ve = list(mesh.structured_iterate_hex())[0]
            tag[ve] = phi
            self._transmute_mesh(mesh, flux_tag, t, tol, element_lib, data_lib, args, 
                                 kwargs)
            y = mesh.mats[0]
            return y        

    def _transmute_mesh(self, mesh, flux_tag, t=None, tol=None, element_lib=None, 
                        data_lib=None, irr_blocks=None, *args, **kwargs):
        """Transmutes a mesh of material into their daughters.

        Parameters
        ----------
        t : float
            Transmutations time [sec].
        tol : float
            Tolerance level for chain truncation.
        mesh : PyNE mesh object
            The mesh to be transmuted
        flux_tag : str
            The name of the tag related to flux on the mesh
        element_lib : str
            Path to element library.
        data_library : str
            The data_library card (see ALARA user's guide).
        irr_blocks : str, optional
            Irradition-related ALARA input blocks.
        args : tuple, optional
            Other arguments ignored for compatibility with other Transmuters.
        kwargs : dict, optional
            Other keyword arguments ignored for compatibility with other 
            Transmuters.
        """
        if t is not None:
            self.t = t
        if tol is not None:
            self.tol = tol
        if element_lib is not None:
            self.element_lib = element_lib
        if data_lib is not None:
            self.data_lib = data_lib
        if irr_blocks is not None:
            self.irr_blocks = irr_blocks

        fluxin = str(uuid.uuid4())
        geom = str(uuid.uuid4())
        matlib = str(uuid.uuid4())

        # build alara input
        alara.mesh_to_fluxin(mesh, flux_tag, fluxin, reverse=False)
        alara.mesh_to_geom(mesh, geom, matlib)
        if self.irr_blocks is None:
            irr_blocks = alara.irradiation_blocks(matlib, element_lib, data_lib,
                                                  ["0 s"], fluxin, self.t, 
                                                  output="number_density",
                                                  truncation=self.tol)
        else:
            irr_blocks = self.irr_blocks

        with open(geom, 'a') as f:
            f.write(irr_blocks)
        # run alara
        p = subprocess.Popen(["alara", geom], stdout=subprocess.PIPE)
        alara_out, err = p.communicate()
        # tag transmuted materials back to mesh
        alara.num_density_to_mesh(alara_out.split("\n"), 'shutdown', mesh)

        if os.path.isfile(fluxin):
            os.remove(fluxin)
        if os.path.isfile(geom):
            os.remove(geom)
        if os.path.isfile(matlib):
            os.remove(matlib)
