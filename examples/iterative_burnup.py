"""Top level wrapper script for iterative burnup calculation.
"""
from itaps import iMesh
import subprocess

from pyne.material import Material
from pyne.mesh import Mesh
from pyne.mcnp import Meshtal, Meshtally, mesh_to_geom
from pyne.alara import mesh_to_fluxin, mesh_to_geom, num_density_to_mesh

mcnp_data_cards = ("               ")

def burnup(mesh, irr_time, power, count)

    geom = mcnp.mesh_to_geom(mesh, 
           title_card="Input for iteration {0}".format(count)

    with open("MCNP_inp_{0}".count, 'w') as f:
        f.write(geom + mcnp_data_cards)

    run_mcnp()

    meshtal = Meshtal("meshtal_{0}".format(count))
    normalize_to_power(meshtal, power)
    mesh_to_fluxin(meshtal, "n_flux", "fluxin_{0}".format(count), reverse=False)
    alara.mesh_to_geom(mesh, "alara_geom_{0}".format(count), 
                       "alara_matlib_{0}".format(count))
    irr_blocks = irradiation_blocks("alara_matlib_{0}".format(count), "isolib", 
                   "FEINDlib CINDER CINDER90 THERMAL", ["0 s"], 
                   "fluxin_{0}".format(count), irr_time)

    with open("alara_geom_{0}".format(count), 'a') as f:
        f.write(irr_blocks)

    alara_out = subprocess.check_output("alara alara_geom_{0}".format(count))

    num_density_to_mesh("alara_out_{0}".format(count), shutdown, mesh)

    return mesh

def normalize_to_power(mesh, power)
    # calculate fission energy per source
    e_per_s = 0
    ves = mesh.structured_iterate_hex('xyz')
    tag = mesh.mesh.getTagHandle("n_flux")
    for ve in ves:
       fluxes = tag[ve]
       e_per_s += 
 
    print("The normalization is {0}".format(norm))   

def gen_reactor():



    

def main( arguments = None ) :

    #Instantiate option parser
    parser = OptionParser(usage='%prog <x_dim>  [options]')

    # 
    #parser.add_option('-o', dest='mesh_output', default='flux_mesh.h5m',\
    #                  help = 'Name of mesh output file, default=%default')

    (opts, args) = parser.parse_args(arguments)


    step = 0
    while step < num_steps
        mesh = burnup(mesh, irr_time, power, count)
        step += 1

if __name__ == __main__:
    main()
