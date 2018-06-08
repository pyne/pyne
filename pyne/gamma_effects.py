"""Purpose:

  This module includes the gamma effects compton edge, backscatter, single escape and double escape 
  By using E as the intial gamma energy, these gamma effects can be calculated.

"""

MEC_2 = 511. #keV # Rest mass energy of an electron #.511 MeV 

def compton_edge(E):
    """Returns compton edge energy from a given gamma energy"""  
    E_compton = E- E/(1+(2*E/MEC_2))
    return E_compton


def backscatter(E):
    """Returns backscatter energy from a given gamma energy"""
    E_backscatter = E/(1+(2*E/MEC_2))
    return E_backscatter


def single_escape(E):
    """Returns single escape energy from a given gamma energy"""
    if (E < 1022):	#keV
        raise RuntimeError('E needs to be equal to or greater than 1022 keV to be valid')
    E_single_escape = E-MEC_2
    return E_single_escape


def double_escape(E):
    """Returns double escape energy from a given gamma energy"""
    if (E < 1022):	#keV
        raise RuntimeError('E needs to be equal to or greater than 1022 keV to be valid')
    E_double_escape = E-(2*MEC_2)
    return E_double_escape
