from pyne import alara
from pyne.transmute.transmuter import Transmuter
from pyne.mesh import Mesh
import uuid

class AlaraTransmuter(Transmuter):
    """A class for transmuting materials using a ALARA."""

    def __init__(self, element_lib, data_lib, t=0.0, phi=0.0, temp=300.0, 
                 tol=1e-7, *args, **kwargs):
        """Parameters
        ----------
        element_lib : str
            Path to element library.
        data_library : str
            The data_library card (see ALARA user's guide).
        args : tuple, optional
            Other arguments ignored for compatibility with other Transmuters.
        kwargs : dict, optional
            Other keyword arguments ignored for compatibility with other 
            Transmuters.
        """
        super(AlaraTransmuter, self).__init__(t, phi, temp, tol)
        self.element_lib = element_lib
        self.data_lib = data_lib

    def transmute(self, x, t=None, phi=None, tol=None, log=None, *args, 
                  **kwargs):
        """Transmutes a material into its daughters.

        Parameters
        ----------
        x : Material or similar
            Input material for transmutation.
        t : float
            Transmutations time [sec].
        phi : float or array of floats
            Neutron flux vector [n/cm^2/sec].  Currently this must either be 
            a scalar or match the group structure of EAF.
        tol : float
            Tolerance level for chain truncation.
        element_lib : str
            Path to element library.
        data_library : str
            The data_library card (see ALARA user's guide).

        Returns
        -------
        y : Material
            The output material post-transmutation.

        """
        if t is not None:
            self.t = t
        if phi is nophi None:
            self.phi = phi
        if tol is not None:
            self.tol = tol
        if element_lib is not None:
            self.element_lib = element_lib
        if data_lib is not None:
            self.data_lib = data_lib

        mesh = Mesh(structured_coords=[[-1,0,1],[-1,0,1],[0,1]], structured=True, 
                    structured_ordering='zyx', mats={0: x})
        flux_tag = "flux"
        tag = mesh.mesh.createTag(flux_tag, 1, float)
        ve = list(mesh.structure_iterate_hex())[0]
        tag[ve] = phi
        transmute_mesh(mesh, flux_tag, t, tol, element_lib, data_lib, args, 
                       kwargs)
        y = mesh.mats[0]
        return y        

    def transmute_mesh(self, mesh, flux_tag, t=None, tol=None, element_lib=None, 
                       data_lib=None, *args, **kwargs):
        """Transmutes a mesh of material into their daughters.

        Parameters
        ----------
        t : float
            Transmutations time [sec].
        tol : float
            Tolerance level for chain truncation.
        mesh :
            The mesh to be transmuted
        flux_tag :
            The name of the tag related to flux on the mesh
        element_lib : str
            Path to element library.
        data_library : str
            The data_library card (see ALARA user's guide).
        """
        if t is not None:
            self.t = t
        if tol is not None:
            self.tol = tol
        if element_lib is not None:
            self.element_lib = element_lib
        if data_lib is not None:
            self.data_lib = data_lib

        fluxin = str(uuid.uuid4())
        geom = str(uuid.uuid4())
        matlib = str(uuid.uuid4())

        # build alara input
        alara.mesh_to_fluxin(mesh, flux_tag, fluxin, reverse=False)
        alara.mesh_to_geom(mesh, geom, matlib)
        irr_blocks = alara.irradiation_blocks(matlib, element_lib, data_lib,
                                              ["0 s"], fluxin, self.t, 
                                              output="number_density",
                                              truncation=self.tol)
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
