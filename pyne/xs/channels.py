from pyne import nucname
from pyne.xs.cache import xs_cache


def sigma_f(nuc, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron fission cross-section for an isotope for a new, lower resolution
    group structure using a higher fidelity flux.  Note that g indexes G, n indexes N, and G < N.

    Parameters
    ----------
    nuc :
        A nuclide name for which to calculate the fission cross-section.

    Keyword Args:
        If any of these are None-valued, values from the cache are used.

        * E_g (sequence of floats): New, lower fidelity energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.
        * E_n (sequence of floats): higher resolution energy group structure [MeV]
          that is of length N+1. Ordered from lowest-to-highest energy.
        * phi_n (sequence of floats): The high-fidelity flux [n/cm^2/s] to collapse the fission 
          cross-section over.  Length N.  Ordered from lowest-to-highest energy.

    Returns:
        * sigma_f_g (numpy array): A numpy array of the collapsed fission cross-section.

    Notes
    -----
    This always pulls the fission cross-section out of the nuc_data library.    

    """
    # Ensure that the low-fidelity group structure is in the cache
    if E_n is not None:
        xs_cache['E_n'] = E_n

    if E_g is not None:
        xs_cache['E_g'] = E_g

    if phi_n is not None:
        xs_cache['phi_n'] = phi_n

    # Get the fission XS
    iso_zz = isoname.mixed_2_zzaaam(iso)
    sigma_f_n_iso_zz = 'sigma_f_n_{0}'.format(iso_zz)
    sigma_f_g_iso_zz = 'sigma_f_g_{0}'.format(iso_zz)

    # Don't recalculate anything if you don't have to
    if sigma_f_g_iso_zz in xs_cache:
        return xs_cache[sigma_f_g_iso_zz]
    else:
        sigma_f_n = xs_cache[sigma_f_n_iso_zz]

    # Perform the group collapse, knowing that the right data is in the cache
    sigma_f_g = partial_group_collapse(sigma_f_n)

    # Put this value back into the cache, with the appropriate label
    xs_cache[sigma_f_g_iso_zz] = sigma_f_g

    return sigma_f_g

