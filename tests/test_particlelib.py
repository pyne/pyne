"""PyNE MCNP tools tests"""
import os
import unittest
import nose
import struct
import warnings

import nose.tools
from nose.tools import assert_almost_equal, assert_equal, assert_true, \
    assert_not_equal, assert_false, assert_raises
from nose.plugins.skip import SkipTest

from pyne.utils import VnVWarning

warnings.simplefilter("ignore", VnVWarning)
try:
    from pyne import mcnp
    from pyne.mcnp import read_mcnp_inp
except ImportError:
    raise SkipTest

from pyne.particlelibrary import ParticleLibrary

def test_new_lib():
    # at the time of writing there were 74 particle 
    # psuedonyms
    lib = ParticleLibrary()
    print len(lib)
    assert_equal(len(lib),74)

def test_pseudo_particles():
    lib = ParticleLibrary()
    pseudo_particle = particle_dict["pseduo"]
    assert_equal(pseudo_particle,"pseudo_particle")

def test_reactor_particles():
    lib = ParticleLibrary()
    # test the aliases of electron/anti-electron
    electron_name = lib["electron"]
    assert_equal(electron_name,"e")
    electron_name = lib["e-"]
    assert_equal(electron_name,"e")
    electron_name = lib["anti-electron"]
    assert_equal(electron_name,"anti-e")
    electron_name = lib["positron"]
    assert_equal(electron_name,"anti-e")
    electron_name = lib["e+"]
    assert_equal(electron_name,"anti-e")

    # test photons
    photon_name = lib["gamma"]  
    assert_equal(photon_name,"photon")
    photon_name = lib["x-ray"]
    assert_equal(photon_name,"photon")
    photon_name = lib["photon"]
    assert_equal(photon_name,"photon")

    # test "regular" baryons
    neutron_name = lib["neutron"]
    assert_equal(neutron_name,"neutron")
    neutron_name = lib["anti-neutron"] 
    assert_equal(neutron_name,"anti-neutron")
    proton_name = lib["proton"] 
    assert_equal(proton_name,"proton")
    proton_name = lib["protium"]
    assert_equal(proton_name,"proton")
    proton_name = lib["p+"] 
    assert_equal(proton_name,"proton")
    proton_name = lib["h"]  
    assert_equal(proton_name,"proton")
    proton_name = lib["anti-proton"] 
    assert_equal(proton_name,"anti-proton")

    # test reactor relevant baryons
    helium_name = lib["helium-4"] 
    assert_equal(helium_name,"alpha")
    helium_name = lib["alpha"] 
    assert_equal(helium_name,"alpha")
    helium_name = lib["helium-3"] 
    assert_equal(helium_name,"helium-3")
    helium_name = lib["helion"]
    assert_equal(helium_name,"helium-3")
    deuterium_name = lib["deuteron"] 
    assert_equal(deuterium_name,"deuteron")
    triton_name = lib["triton"] 
    assert_equal(triton_name,"triton")

def test_high_energy_leptons():
    lib = ParticleLibrary()
    muon_name =  lib["muon"]
    assert_equal(muon_name,"mu")
    muon_name = lib["anti-muon"]
    assert_equal(muon_name,"anti-mu")
    tauon_name = lib["tauon"] 
    assert_equal(tauon_name,"tau")
    tauon_name = lib["anti-tauon"]
    assert_equal(tauon_name,"anti-tau")
    neutrino_name = lib["electron-neutrino"]
    assert_equal(neutrino_name,"nu_e")
    neutrino_name =  lib["anti-electron-neutrino"]
    assert_equal(neutrino_name,"anti-nu_e")
    neutrino_name =  lib["muon-neutrino"] 
    assert_equal(neutrino_name,"nu_mu")
    neutrino_name =  lib["anti-muon-neutrino"]
    assert_equal(neutrino_name,"anti-nu_mu")
    neutrino_name =  lib["tauon-neutrino"]
    assert_equal(neutrino_name,"tau_mu")
    neutrino_name =  lib["anti-tauon-neutrino"]
    assert_equal(neutrino_name,"anti-tau_mu")
