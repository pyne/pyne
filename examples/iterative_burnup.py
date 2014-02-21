"""Top level wrapper script for iterative burnup calculation.
"""
from itaps import iMesh
import subprocess
from optparse import OptionParser

from pyne.material import Material, from_atom_frac
from pyne.mesh import Mesh
from pyne import mcnp
from pyne import alara
from pyne import nucname

J_PER_MEV = 1.602177E-22

mcnp_data_cards = (
"""MODE:n
IMP:n 1 26R 0
KCODE 1000 0.8 1 2
KSRC 1E-9 0 0
FMESH4:n origin=-5,-5,-5
      imesh = 5
      iints = 3
      jmesh = 5
      jints = 3
      kmesh = 5
      kints = 3
      $ lowest energy bin commented out, MCNP uses 0
      emesh=  $ 1.00000E-11
                 5.00000E-09 1.00000E-08 1.50000E-08 2.00000E-08
                 2.50000E-08 3.00000E-08 3.50000E-08 4.20000E-08
                 5.00000E-08 5.80000E-08 6.70000E-08 8.00000E-08
                 1.00000E-07 1.52000E-07 2.51000E-07 4.14000E-07
                 6.83000E-07 1.12500E-06 1.85500E-06 3.05900E-06
                 5.04300E-06 8.31500E-06 1.37100E-05 2.26000E-05
                 3.72700E-05 6.14400E-05 1.01300E-04 1.67000E-04
                 2.75400E-04 4.54000E-04 7.48500E-04 1.23400E-03
                 2.03500E-03 2.40400E-03 2.84000E-03 3.35500E-03
                 5.53100E-03 9.11900E-03 1.50300E-02 1.98900E-02
                 2.55400E-02 4.08700E-02 6.73800E-02 1.11100E-01
                 1.83200E-01 3.02000E-01 3.88700E-01 4.97900E-01
                 6.39279E-01 8.20850E-01 1.10803E+00 1.35335E+00
                 1.73774E+00 2.23130E+00 2.86505E+00 3.67879E+00
                 4.96585E+00 6.06500E+00 1.00000E+01 1.49182E+01
                 1.69046E+01 2.00000E+01 2.50000E+01
FMESH14:n origin=-5,-5,-5
      imesh = 5
      iints = 1
      jmesh = 5
      jints = 1
      kmesh = 5
      kints = 1
FM14 -1 0 -6 -8
""")

def burnup(mesh, irr_time, power, xsdir_path, count):

    # create m
    geom = mcnp.mesh_to_geom(mesh, 
           title_card="Input for iteration {0}".format(count))

    with open("MCNP_inp_{0}".format(count), 'w') as f:
        f.write(geom + mcnp_data_cards)

    # run mcnp
    #subprocess.call("mcnp5 i=MCNP_inp_{0} meshtal=meshtal_{0} xsdir={1}".format(count, xsdir_path), shell=True)
    subprocess.call("mcnp5 i=MCNP_inp_{0} meshtal=meshtal_{0}".format(count, xsdir_path), shell=True)

    # use mcnp out to create alara input files
    meshtal = mcnp.Meshtal("meshtal_{0}".format(count))
    normalize_to_power(meshtal, power)
    alara.mesh_to_fluxin(meshtal.tally[4], "n_result",
                         "fluxin_{0}".format(count), reverse=False)
    alara.mesh_to_geom(mesh, "alara_geom_{0}".format(count), 
                       "alara_matlib_{0}".format(count))
    irr_blocks = alara.irradiation_blocks("alara_matlib_{0}".format(count), "../isolib", 
                   "FEINDlib CINDER ../CINDER90 THERMAL", ["0 s"], 
                   "fluxin_{0}".format(count), irr_time, output="number_density",truncation=1E-6)

    with open("alara_geom_{0}".format(count), 'a') as f:
        f.write(irr_blocks)

    # run alara
    #p = subprocess.Popen(["alara", "alara_geom_{0}".format(count)])
    p = subprocess.Popen(["alara", "alara_geom_{0}".format(count)], stdout=subprocess.PIPE)
    #p = subprocess.Popen(["alara", "alara_geom_{0}".format(count)], shell=True, stderr=None)
    alara_out, err = p.communicate()
    # tag transmuted materials back to mesh
    alara.num_density_to_mesh(alara_out.split("\n"), 'shutdown', mesh)

    xsdir = mcnp.Xsdir(xsdir_path)
    #nucs = [nucname.id(x.name.split('.')[0]) for x in xsdir.tables if x.name.split('.')[0].isdigit()]
    nucs = [nucname.id(x.name.split('.')[0]) for x in xsdir.tables if nucname.isnuclide(x.name.split('.')[0])]
    bad_nucs = [60150000, 60140000, 80180000, 410960000, 410980000, 411000000, 421010000, 410970000]
    for bad_nuc in bad_nucs:
        if bad_nuc in nucs:
            nucs.remove(bad_nuc)

    for i, mat, ve in mesh:
       comp = mesh.mats[i].comp
       density = mesh.mats[i].density
       new_mat = Material(nucvec=dict(comp))
       mesh.mats[i]  = new_mat[nucs]
       mesh.mats[i].density = density

def normalize_to_power(meshtal, power):
    # calculate fission energy per source
    flux_mesh = meshtal.tally[4]
    e_mesh = meshtal.tally[14]

    e_ve = list(e_mesh.structured_iterate_hex())[0]
    tot_vol = e_mesh.structured_hex_volume(0, 0, 0)
    e = e_mesh.mesh.getTagHandle("n_result")[e_ve]*tot_vol
    print e
    norm = power/J_PER_MEV/e
    ves = flux_mesh.structured_iterate_hex('xyz')
    tag = flux_mesh.mesh.getTagHandle("n_result")
    for ve in ves:
       tag[ve] = [x*norm for x in tag[ve]]
 
    print("The normalization is {0} [n/s]".format(norm))   


def gen_reactor(mesh, fuel_idx):

    ves = mesh.structured_iterate_hex(mesh.structured_ordering)
    for i, ve in enumerate(ves):
       if i in fuel_idx:
           mesh.mats[i] = from_atom_frac({'U235': 0.045, 'U238': 0.955, 'O16': 2.0}, density=10.7)
       else:
           mesh.mats[i] = from_atom_frac({'H1': 2.0, 'O16': 1.0}, density=1.0)

def main(arguments=None):

    #Instantiate option parser
    parser = OptionParser(usage =
             '%prog <power [W]> <time [s]> <num steps> <xsdir path> [options]')

    # 
    #parser.add_option('-o', dest='mesh_output', default='flux_mesh.h5m',\
    #                  help = 'Name of mesh output file, default=%default')

    (opts, args) = parser.parse_args(arguments)
    power = float(args[0])
    time = float(args[1])
    num_steps = int(args[2])
    irr_time = str(time/num_steps) + ' s'
    xsdir_path = args[3]
    xsdir = mcnp.Xsdir(xsdir_path)

    # create reactor that spans [[-5, 5], [-5, 5],[-5, 5]] with fuel volume
    # elements near the center
    # big reactor
    #fuel_idx = range(399, 500)
    #coords = range(-5, 6)

    # small reactor
    fuel_idx = range(0,15)
    coords = [-5000, -300, 300, 500]
    mesh = Mesh(structured_coords = [coords, coords, coords], structured=True)
    gen_reactor(mesh, fuel_idx)

    nucs = xsdir.nucs()
    step = 0

    while step < num_steps:
        #for i, mat, ve in mesh:
        #    mesh.mats[i]  = mesh.mats[i][nucs]

        burnup(mesh, irr_time, power, xsdir_path, step)
        step += 1

main()
