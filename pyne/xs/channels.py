from pyne import nucname
from pyne.xs.cache import xs_cache
from pyne.xs.models import group_collapse

def _prep_cache(E_g=None, E_n=None, phi_n=None):
    """Ensures that certain values are in the cache safely."""
    if E_n is not None:
        xs_cache['E_n'] = E_n

    # needs to follow E_n to calc collapse matrix properly
    if E_g is not None:
        xs_cache['E_g'] = E_g

    if phi_n is not None:
        xs_cache['phi_n'] = phi_n


def sigma_f(nuc, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron fission cross-section for a nuclide for a new, 
    lower resolution group structure using a higher fidelity flux.  Note that 
    g indexes G, n indexes N, and G < N.  If any of these are None-valued, 
    values from the cache are used.  The energy groups and fluxes are normally 
    ordered from highest-to-lowest energy.

    Parameters
    ----------
    nuc :
        A nuclide name for which to calculate the fission cross-section.
    E_g : array-like of floats, optional
        New, lower fidelity energy group structure [MeV] that is of length G+1. 
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : array-like of floats, optional
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over.  
        Length N.  

    Returns
    -------
    sigma_f_g : ndarray 
        An array of the collapsed fission cross-section.

    Notes
    -----
    This always pulls the fission cross-section out of nuc_data library.    

    """
    _prep_cache(E_g, E_n, phi_n)

    # Get the fission XS
    nuc_zz = nucname.zzaaam(nuc)
    sigma_f_n_nuc_zz = ('sigma_f_n', nuc_zz)
    sigma_f_g_nuc_zz = ('sigma_f_g', nuc_zz)

    # Don't recalculate anything if you don't have to
    if sigma_f_g_nuc_zz in xs_cache:
        return xs_cache[sigma_f_g_nuc_zz]
    else:
        sigma_f_n = xs_cache[sigma_f_n_nuc_zz]

    # Perform the group collapse, knowing that the right data is in the cache
    sigma_f_g = group_collapse(sigma_f_n, xs_cache['phi_n'], phi_g=xs_cache['phi_g'], 
                               partial_energies=xs_cache['partial_energy_matrix'])

    # Put this value back into the cache, with the appropriate label
    xs_cache[sigma_f_g_nuc_zz] = sigma_f_g

    return sigma_f_g

