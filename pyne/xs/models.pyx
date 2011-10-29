"""This module provides physical cross-section models and helper functions."""

cimport numpy as np
import numpy as np

from pyne cimport nucname
from pyne import nucname

from scipy import constants
from scipy.special import erf

# Integration imports
#from scipy import integrate
#import metasci.mathematics.integrate as msmintegrate

# Bolzmann's constant in MeV/K
k = constants.physical_constants['Boltzmann constant in eV/K'][0] * (1.0E-6)

# Neutron mass in amu
m_n = constants.physical_constants['neutron mass in u'][0]


##############################
### Partial group collapse ###
##############################

def partial_energy_matrix_mono(np.ndarray[np.float64_t, ndim=1] E_g, np.ndarray[np.float64_t, ndim=1] E_n, int slope=-1):
    """Gerenates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.
    Here, both of the energy arrays must be monotonic. This is useful for performing group collapses.

    Parameters
    ----------
    E_g : 1d numpy float array 
        Lower resolution energy group structure [MeV] that is of length G+1. 
        Ordered based on slope.
    E_n : 1d numpy float array 
        Higher resolution energy group structure [MeV] that is of length N+1. 
        Ordered based on slope.
    slope : int, optional
        Gives the monotonicity of E_g and E_n.  If positive, then they are 
        monotonicly increasing (lowest-to-highest).  If negative, they are
        monotonicly decreasing (highest-to-lowest).

    Returns
    -------
    pem : 2d numpy float array of fractions 
        This is a GxN sized matrix that when dotted with a high-resolution 
        flux (or cross section) produces a low-resolution flux (or cross section).
    """
    # Some convienence parameters
    cdef Py_ssize_t G, N, g, lig, uig
    cdef np.ndarray[np.float64_t, ndim=2] pem

    G = E_g.shape[0] - 1
    N = E_n.shape[0] - 1

    index_E_n = np.arange(N+1)

    # Get the interior points for each gth group in n-space
    if slope < 0:
        inner_mask = np.array([(E_n <= E_g[g]) & (E_g[g+1] <= E_n) for g in range(G)])
    elif 0 < slope:
        inner_mask = np.array([(E_g[g] <= E_n) & (E_n <= E_g[g+1]) for g in range(G)])
    else:
        raise ValueError("slope must be positive or negative.")

    # Get the upper and lower nth index for every gth group
    lower_index = np.array([index_E_n[inner_mask[g]][0] for g in range(G)])
    upper_index = np.array([index_E_n[inner_mask[g]][-1] for g in range(G)])

    # Convert the mask to initialize the partial enery matrix
    # Hack off the last index of the mask to make the right size
    pem = np.array(inner_mask[:, :-1], dtype=float)

    # Check for partial contibutions at the edges
    for g in range(G):
        # Lower bound
        lig = lower_index[g]
        if lig != 0:
            pem[g,lig-1] = (E_n[lig] - E_g[g]) / (E_n[lig] - E_n[lig-1])

        # Upper bound
        uig = upper_index[g]
        if uig < N:
            pem[g,uig] = (E_g[g+1] - E_n[uig]) / (E_n[uig+1] - E_n[uig])

    return pem



def partial_energy_matrix(E_g, E_n):
    """Gerenates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.
    The group structures must have the same monotonicity. This is useful for performing group collapses.

    Parameters
    ----------
    E_g : sequence of floats
        Lower resolution energy group structure [MeV] that is of length G+1. 
    E_n : sequence of floats
        Higher resolution energy group structure [MeV] that is of length N+1. 

    Returns
    -------
    pem : 2d numpy float array of fractions 
        This is a GxN sized matrix that when dotted with a high-resolution 
        flux (or cross section) produces a low-resolution flux (or cross section).
    """
    cdef np.ndarray[np.float64_t, ndim=2] pem

    E_g = np.asarray(E_g, dtype=float)
    E_n = np.asarray(E_n, dtype=float)

    if (E_g[:-1] > E_g[1:]).all() and (E_n[:-1] > E_n[1:]).all():
        # Both energy arrays are monotonically decreasing
        assert E_g[0] <= E_n[0]
        assert E_n[-1] <= E_g[-1]
        pem = partial_energy_matrix_mono(E_g, E_n, -1)
    elif (E_g[:-1] < E_g[1:]).all() and (E_n[:-1] < E_n[1:]).all():
        # Both energy arrays are monotonically increasing
        assert E_n[0] <= E_g[0]
        assert E_g[-1] <= E_n[-1]
        pem = partial_energy_matrix_mono(E_g, E_n, 1)
    else:
        raise ValueError("E_g and E_n are not both monotonic in the same direction.")

    return pem



######################
### Group Collapse ###
######################
def phi_g(E_g, E_n, phi_n):
    """Calculates a lower resolution flux, phi_g, from a lower resolution group stucture E_g, 
    a higher resolution groups E_n, and a higher resolution flux phi_n.

    Parameters
    ----------
    E_g : sequence of floats 
        Lower resolution energy group structure [MeV] that is of length G+1. 
    E_n : sequence of floats 
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : sequence of floats
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over (length N).  

    Returns
    -------
    phi_g : numpy array of floats 
        The flux collapsed to G energy groups.
    """
    pem = partial_energy_matrix(E_g, E_n)
    phi_g = np.dot(pem, phi_n)
    return phi_g



#######################
### Physical models ###
#######################

def chi(E):
    """Calculates the fission neutron spectrum (frequency) at energy E.
    E may be either a float or an array of floats. This is based off of 
    the values for U-235, which are representative for other isotopes.  
    See Lamarsh or 'Comparison of prompt-fission neutron multiplicities 
    and energy spectra for intermediate energy proton-and neutron-induced 
    fission' --Oleg Batenkov, Georgy Boikov, Vilen Eismont, Mikhail Majorov, 
    Sergey Soloviev, Jan Blomgren, and Walter Loveland.
    """
    chi_E = 0.453 * np.exp(-1.036*E) * np.sinh(np.sqrt(2.29*E))
    return chi_E


#
# Scattering cross section
#

def alpha(E_prime, E, theta, M_A=1.0, T=300.0):
    """Scattering kernel alpha value.

    .. math::

        \\alpha = \\frac{E^\\prime + E - 2\\sqrt{E^\\prime E}\\cos\\theta}{\\frac{M_A}{m_n}kT}

    Parameters
    ----------
    E_prime : float (or array)
        The exiting energy of the neutron after scattering event [MeV].
    E : float (or array)
        The incident energy of the neutron prior to scattering event [MeV].
    theta : float (or array) 
        Scattering angle in [radians].
    M_A : float (or array), optional
        Atomic mass of the target nucleus [amu].
    T : float (or array), optional 
        Tempurature of the target material [kelvin].

    Returns
    -------
    a : float (or array)
        alpha value
    """
    a = (E_prime + E - 2 * np.sqrt(E_prime*E) * np.cos(theta)) / (k * T * M_A / m_n)
    return a


def beta(E_prime, E, T=300.0):
    """Scattering kernel beta value.

    .. math::

        \\beta = \\frac{E^\\prime - E}{kT}

    Parameters
    ----------
    E_prime : float (or array) 
        The exiting energy of the neutron after scattering event [MeV].
    E : float (or array) 
        The incident energy of the neutron prior to scattering event [MeV].
    T : float (or array), optional 
        Tempurature of the target material [kelvin].

    Returns
    -------
    b : float 
        beta value.
    """
    b = (E_prime - E) / (k*T)
    return b


def alpha_at_theta_0(E_prime, E, M_A=1.0, T=300.0):
    """Scattering kernel alpha value at the lower bound of the scattering angle.

    .. math::

        \\alpha_{\\theta=0} = \\frac{E^\\prime + E - 2\\sqrt{E^\\prime E}}{\\frac{M_A}{m_n}kT}

    Parameters
    ----------
    E_prime : float (or array)
        The exiting energy of the neutron after scattering event [MeV].
    E : float (or array)
        The incident energy of the neutron prior to scattering event [MeV].
    M_A : float (or array), optional
        Atomic mass of the target nucleus [amu].
    T : float (or array), optional 
        Tempurature of the target material [kelvin].

    Returns
    -------
    a : float (or array)
        alpha value with theta = 0.
    """
    a = (E_prime + E - 2 * np.sqrt(E_prime*E)) / (k * T * M_A / m_n)
    return a


def alpha_at_theta_pi(E_prime, E, M_A=1.0, T=300.0):
    """Scattering kernel alpha value at the upper bound of the scattering angle.

    .. math::

        \\alpha_{\\theta=\\pi} = \\frac{E^\\prime + E + 2\\sqrt{E^\\prime E}}{\\frac{M_A}{m_n}kT}

    Parameters
    ----------
    E_prime : float (or array)
        The exiting energy of the neutron after scattering event [MeV].
    E : float (or array)
        The incident energy of the neutron prior to scattering event [MeV].
    M_A : float (or array), optional
        Atomic mass of the target nucleus [amu].
    T : float (or array), optional 
        Tempurature of the target material [kelvin].

    Returns
    -------
    a : float (or array)
        alpha value with theta = pi.
    """
    a = (E_prime + E + 2 * np.sqrt(E_prime*E)) / (k * T * M_A / m_n)
    return a


def one_over_gamma_squared(E):
    """The inverse of the Lorentz factor sqared. Sometimes used as 
    a realitivistic correction factor for the bound scattering length.

    .. math::

        \\frac{1}{\\gamma^2} = \\left( 1 - \\frac{2E}{931.46 \\cdot m_n} \\right)

    Parameters
    ----------
    E : float (or array)
        The incident energy of the neutron prior to scattering event [MeV].

    Returns
    -------
    inv_g2 : float (or array)
        Inverse of gamma squared.
    """
    inv_g2 = 1.0 - (E / (465.73 * m_n))
    return inv_g2


#
# This needs more thought, much like all of the scattering models
#

#def d2sigma_s_dE_prime_dOmega(E_prime, E, theta, b=1.0, M_A=1.0, T=300.0):
#    """Computes the double differential total scattering cross section from the equation.  
#    This is likely only valid in the thermal region.
#
#    .. math::
#
#        \\frac{d\\sigma_s(E)}{dE^\\prime} = \\mbox{[1]}
#
#    FIXME: I am untested
#
#    Parameters
#    ----------
#    E_prime : float (or array)
#        The exiting energy of the neutron after scattering event [MeV].
#    E : float (or array)
#        The incident energy of the neutron prior to scattering event [MeV].
#    theta : float (or array) 
#        Scattering angle in [radians].
#    b : float (or array), optional
#        The bound scattering length of the target nucleus.
#    M_A : float (or array), optional
#        Atomic mass of the target nucleus [amu].
#    T : float (or array), optional 
#        Tempurature of the target material [kelvin].
#
#    Refs:
#        1. Mattes M, Keinert J. Thermal neutron scattering data for the moderator 
#           materials H2O, D2O and ZrHx in ENDF-6 format and as ACE library for 
#           MCNP (X) codes. IAEA report INDC (NDS)-0470. 2005;(April). Available at: 
#           http://200.136.52.101/reports-new/indc-reports/indc-nds/indc-nds-0470.pdf.
#    """
#    kT = k * T
#    rcf = one_over_gamma_squared(E)
#
#    _alpha = alpha(E_prime, E, theta, M_A, T) 
#    _beta = beta(E_prime, E, T)
#
#    power_term = (_beta/2.0) + (_alpha/4.0) + (_beta**2 / (4.0 * _alpha))
#
#    return rcf * (b**2 / kT) * np.sqrt((np.pi * E_prime) / (_alpha * E)) * np.exp(-power_term)


    
#
# Failed analytic soltion to the integral of the double differential 
# scattering cross-section over all solid angles.  This is because the 
# S(alpha, beta) term was only valid in the thermal regime but was used 
# even up to fast energies.
#

#def dsigma_s_dE_prime(E_prime, E, b=1.0, M_A=1.0, T=300.0):
#    """Computes the differential total scattering cross section from an analytic
#    solution to the integral of the double-differentional scattering cross section, 
#    integrated over all solid angles.
#
#    .. math::
#
#        \\frac{d\sigma_s(E)}{dE^\prime} = b^2 \left( 1 - \\frac{2E}{931.46 \\cdot m_n} \\right) 
#            \\frac{e^{-\\frac{\\beta + |\\beta|}{2}}}{2E} \\frac{M_A}{m_n} Q
#
#        Q = \left( \mbox{Erf}\left(\\frac{|\\beta| - \\alpha_{\\theta=0}}{2 \sqrt{\\alpha_{\\theta=0}}}\\right) 
#            - \mbox{Erf}\left(\\frac{|\\beta| - \\alpha_{\\theta=\pi}}{2 \sqrt{\\alpha_{\\theta=\pi}}}\\right) \\right) 
#            - e^{-\\frac{|\\beta|}{2}} \left( \mbox{Erf}\left(\\frac{|\\beta| + \\alpha_{\\theta=0}}{2 \sqrt{\\alpha_{\\theta=0}}}\\right) 
#            - \mbox{Erf}\left(\\frac{|\\beta| + \\alpha_{\\theta=\pi}}{2 \sqrt{\\alpha_{\\theta=\pi}}}\\right) \\right)
#
#    Args:
#        * E_prime (float): The exiting energy of the neutron after the 
#          scattering event [MeV].
#        * E (float): The incident energy of the neutron prior to the 
#          scattering event [MeV].
#
#    Keyword Args:
#        * b (float): The bound scattering length of the target nucleus.
#        * M_A (float): Atomic mass of the target nucleus [amu].
#        * T (float): Tempurature of the target material [kelvin].
#    """
#    kT = k * T
#    rcf = one_over_gamma_squared(E)
#
#    alpha_lower = alpha_given_theta_0(E_prime, E, M_A, T)
#    alpha_upper = alpha_given_theta_pi(E_prime, E, M_A, T)
#    _beta = beta(E_prime, E, T)
#    abs_beta = np.abs(_beta)
#
#    Q = erf((abs_beta - alpha_lower) / (2.0 * np.sqrt(alpha_lower))) - \
#        erf((abs_beta - alpha_upper) / (2.0 * np.sqrt(alpha_upper))) + \
#        np.exp(-abs_beta/2.0) * (erf((abs_beta + alpha_lower) / (2.0 * np.sqrt(alpha_lower))) - \
#        erf((abs_beta + alpha_upper) / (2.0 * np.sqrt(alpha_upper))))
#
#    deriv = (rcf * b**2) * (np.exp(-(_beta + abs_beta)/2.0) / (2.0*E)) * (M_A / m_n) * Q * np.exp(-np.sqrt(E))
#
#    return deriv


def E_prime_min(E, M_A=1.0):
    """The minimum possible exiting enegy of a neuron after a scattering collision. 
    This is based on the incident energy and the mass of the target.
    For a proof, use the conservation of energy and momentum.

    .. math::

        \\mbox{min}(E^\\prime) = \\left(\\frac{M_A - m_n}{M_A + m_n}\\right)^2 E

    
    Parameters
    ----------
    E : float (or array)
        The incident energy of the neutron prior to scattering event [MeV].
    M_A : float (or array), optional
        Atomic mass of the target nucleus [amu].

    Returns
    -------
    min_E_prime : float (or array)
        Minimum exiting energy.
    """
    alph = ((M_A - m_n)/(M_A + m_n))**2 * E
    min_E_prime = alph * E
    return min_E_prime


#
# Below was the average maximal value
#

#def E_prime_max(E, T=300.0):
#    """The maximum possible exiting enegy of a neuron after a scattering collision. 
#    This is based on the incident energy and the tempurature of the target.
#    The neutron can gain no more energy than the kinetic energy of the target.
#    In a macroscopic system, this is on average equal to kT.
#
#    .. math::
#
#        \mbox{max}(E^\prime) = E + kT
#
#    Args:
#        * E (float): The incident energy of the neutron prior to the 
#          scattering event [MeV].
#
#    Keyword Args:
#        * T (float): Tempurature of the target material [kelvin].
#    """
#    max_E = E + (k * T)
#    return max_E


def sigma_s_const(b):
    """Computes the constant scattering cross-section based on the 
    scattering length.

    .. math::

        \\sigma_s = 4 \\pi b^2

    Parameters
    ----------
    b : float (or array)
        The bound scattering length [cm] of the target nucleus.

    Returns
    -------
    sig_s : float (or array) 
        The micorscopic scattering cross-section [barns].

    See Also
    --------
    pyne.data.b : scattering length data.
    """
    sig_s = (4E24) * np.pi * (b**2)
    return sig_s


#
# These definitely need more thought
#

#def sigma_s_E(E, b=1.0, M_A=1.0, T=300.0):
#    """Computes the total scattering cross section from an empirical model.
#
#    .. math::
#
#        \\sigma_s(E) = 4 \\pi b^2 \\cdot \\left( 1 - \\frac{2E}{931.46 \\cdot m_n} \\right) \\cdot
#                      \\left( 1 + \\frac{m_n}{M_A} \\frac{kT}{E} \\cdot e^{-\\frac{M_A}{m_n}\\frac{E}{kT}} \\right) 
#                      \\cdot \\left( 1 - \\mbox{Exp}\\left[-\\sqrt{\\frac{0.1}{E}}\\right] \\right)
#
#    Args:
#        * E (float): The incident energy of the neutron prior to the 
#          scattering event [MeV].
#
#    Keyword Args:
#        * b (float): The bound scattering length of the target nucleus.
#        * M_A (float): Atomic mass of the target nucleus [amu].
#        * T (float): Tempurature of the target material [kelvin].
#    """
#    kT_over_AE = k * T / ((M_A / m_n) * E)
#
#    sig_s = sigma_s_const(b)
#    rcf = one_over_gamma_squared(E)
#    
#    sig_s_E = (rcf * sig_s) * (1.0 + kT_over_AE * np.exp(-1.0/kT_over_AE)) * (1.0 - np.exp(-np.sqrt(0.1/E))) 
#
#    return sig_s_E


#def P(E, E_prime, M_A=1.0, T=300.0):
#    """Calculates the neutron scattering transfer probability from
#
#    .. math::
#        P(E \\to E^\\prime) = \\frac{1}{E + kT - \\left( \\frac{M_A - m_n}{M_A + m_n} \\right)^2 E}
#
#    as long as 
#
#    .. math::
#        \\left( \\frac{M_A - m_n}{M_A + m_n} \\right)^2 E \\le E^\\prime \\le E + kT
#
#    The probability is zero outside of this range.
#
#    Args:
#        * E (float): The incident energy of the neutron prior to the 
#          scattering event [MeV].
#        * E_prime (float): The exiting energy of the neutron after the 
#          scattering event [MeV].
#
#    Keyword Args:
#        * M_A (float): Atomic mass of the target nucleus [amu].
#        * T (float): Tempurature of the target material [kelvin].
#
#    Returns:
#        * p (float): A transfer probablity.
#    """
#    # Find bounds
#    E_prime_lower = E_prime_min(E, M_A)    
#    E_prime_upper = E_prime_max(E, T)    
#
#    if E_prime_lower <= E_prime <= E_prime_upper:
#        p = 1.0 / (E_prime_upper - E_prime_lower)
#    else:
#        p = 0.0
#
#    return p
    
