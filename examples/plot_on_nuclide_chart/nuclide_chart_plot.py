import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyne.data as data
import pyne.nucname as nucname
import pyne.material as material

def get_nuclide_list():
    """
    Get list of all nuclides that PyNE knows about in its atomic mass table.
    """

    nuclide_list = []

    # perform a query to populate the atomic_mass_map
    if (data.atomic_mass(1001)):
        nuclide_list = [key for key, _ in data.atomic_mass_map.items()]

    return nuclide_list

def nuclide_chart_plot(fresh_material, full_nuclide_list, decay_time=0, response="inventory"):
    """Generate a plot of data, provided per nuclide, in a visualiation
    that looks like a typical chart of the nuclides.

    Based on implementation suggested by Nick Touran
    https://terrapower.github.io/armi/gallery/framework/run_chartOfNuclides.html
    
    input
    -----
    nuclide_datamap: dict(nucname: doube) 
                     provides data for each of a set of nuclides given as nucnames
    """
    material = fresh_material

    one_day = 24*60*60
    if decay_time > one_day:
        material = fresh_material.decay(decay_time)

    if response == "decay_heat":
        material = material.decay_heat()
    elif response == "dose":
        material = material.dose_per_g("ext_air")

    plot_data = []
    plot_no_data = []
    for nuclide in full_nuclide_list:
        Z = nucname.znum(nuclide)
        N = nucname.anum(nuclide) - Z
        if N > 0:
            if nuclide in material and material[nuclide] > 0:
                plot_data.append((Z, N, material[nuclide]))
            else:
                plot_no_data.append((Z, N, 0))
    z, n, c = zip(*plot_data)
    zn, nn, _ = zip(*plot_no_data)
    plt.figure(figsize=(15, 8))
    plt.scatter(nn, zn, c="gray", marker=".", s=1)
    plt.scatter(n, z, c=c, marker="s", s=6, 
                norm='log', vmin=np.max(c)/1e6,
                cmap=mpl.colormaps['Oranges'])
    plt.title("Chart of the nuclides")
    plt.xlabel("Number of neutrons (N)")
    plt.ylabel("Number of protons (Z)")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.colorbar()
    plt.show()

def test_code():

    all_nuclide_list = get_nuclide_list()

    U_nat = material.Material({92235: 0.05, 92238: 0.95})

    U_decay = U_nat

    nuclide_chart_plot(U_nat, all_nuclide_list, decay_time=1e19, response="decay_heat")

if __name__ == "__main__":
    test_code()