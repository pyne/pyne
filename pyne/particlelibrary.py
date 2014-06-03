#!/usr/bin/env python

from warnings import warn
from pyne.utils import VnVWarning

"""Module for the IO of particle type related data, the idea being that there
are several functions to query the particle library

"""

class ParticleLibrary(object):
    def __init__(self, values=dict()):
        _lib = { }
        self._lib = self._library_setup()
        print len(_lib)

    def __contains__(self, key):
        return key in self._lib

    def __len__(self):
        return len(self._lib)

    def __iter__(self):
        return iter(self._lib)

    def __getitem__(self, key):
        return self._lib[key]

    def __delitem__(self, key):
        del self._lib[key]

    def _library_setup(self):
        """Setup the particle library 
        """
        particle_dict = { }
        # pseudo particle
        particle_dict["pseduo"] = "pseudo"
        
        # leptons
        particle_dict["electron"] = "e"
        particle_dict["anti-electron"]="anti-e"     
        particle_dict["positron"]="anti-e"     
        particle_dict["e-"]="e"     
        particle_dict["e+"]="anti-e"     
        particle_dict["muon"] = "mu"
        particle_dict["anti-muon"]="anti-mu"        
        particle_dict["tauon"] = "tau"
        particle_dict["anti-tauon"]="anti-tau"
        particle_dict["electron-neutrino"] = "nu_e"
        particle_dict["anti-electron-neutrino"]="anti-nu_e"
        particle_dict["muon-neutrino"] = "nu_mu"
        particle_dict["anti-muon-neutrino"]="anti-nu_mu"        
        particle_dict["tauon-neutrino"] = "tau_mu"
        particle_dict["anti-tauon-neutrino"]="anti-tau_mu"

        # gauge bosons
        particle_dict["gamma"] = "photon"
        particle_dict["x-ray"] = "photon"
        particle_dict["photon"] = "photon"
        particle_dict["optical-photon"] = "optical-photon"

        # baryons
        particle_dict["neutron"] = "neutron"
        particle_dict["anti-neutron"] = "anti-neutron"
        particle_dict["proton"] = "proton"
        particle_dict["protium"] = "proton"
        particle_dict["p+"] = "proton"
        particle_dict["h"] = "proton"
        particle_dict["anti-proton"] = "anti-proton"
        # special reactor physics baryons
        particle_dict["helium-4"] = "alpha"
        particle_dict["alpha"] = "alpha"
        particle_dict["helium-3"] = "helium-3"
        particle_dict["helion"] = "helium-3"
        particle_dict["deuteron"] = "deuteron"
        particle_dict["triton"] = "triton"
        

        # the baryon list beyond the "standard set" are those defaulted
        # to by the fluka naming scheme
        particle_dict["lambda"] = "lambda"
        particle_dict["anti-lambda"] = "anti-lambda"
        particle_dict["charmed-lambda"] = "lambda_c^+"
        particle_dict["bottom-lambda"] = "lambda_b^0"
        particle_dict["charmed-anti-lambda"] = "anti-lambda_c^-"

        particle_dict["kaon-long"] = "K-long"
        particle_dict["kaon-short"] = "K-short"
        particle_dict["kaon-zero"] = "K_0"
        particle_dict["anti-kaon-zero"] = "anti-K_0"

        particle_dict["sigma-zero"] = "sigma^0"
        particle_dict["sigma-plus"] = "sigma^+"
        particle_dict["sigma-minus"] = "sigma^-"
        particle_dict["anti-sigma-zero"] = "anti-sigma^0"
        particle_dict["anti-sigma-plus"] = "anti-sigma^+"
        particle_dict["anti-sigma-minus"] = "anti-sigma^-"

        particle_dict["xi-zero"] = "xi^0"
        particle_dict["anti-xi-zero"] = "anti-xi^0"
        particle_dict["xi-minus"] = "xi^-"
        particle_dict["xi-plus"] = "xi^+"
        particle_dict["charmed-xi-plus"] = "xi_c^+"
        particle_dict["charmed-xi-zero"] = "xi_c^0"
        particle_dict["charmed-xi-primed-plus"] = "xi-primed_c^+"
        particle_dict["charmed-xi-primed-zero"] = "xi-primed_c^0"
        particle_dict["charmed-anti-xi-plus"] = "anti-xi_c^+"
        particle_dict["charmed-anti-xi-zero"] = "anti-xi_c^0"
        particle_dict["charmed-anti-xi-primed-plus"] = "anti-xi-primed_c^+"
        particle_dict["charmed-anti-xi-primed-zero"] = "anti-xi-primed_c^0"

        particle_dict["omega-minus"] = "Omega^-"        
        particle_dict["anti-omega"] = "anti-Omega"        
        particle_dict["charmed-omega-zero"] = "Omega_c^0"        
        particle_dict["charmed-omega-minus"] = "Omega_c^-"        
        particle_dict["anti-omega-charmed-zero"] = "anti-Omega_c^0"        

        particle_dict["d-plus"] = "D^+"        
        particle_dict["d-minus"] = "D^-" 
        particle_dict["d-zero"] = "D^0" 
        particle_dict["anti-d-zero"] = "anti-D^0" 
        particle_dict["d_s-plus"] = "D_s^+" 
        particle_dict["d_s-minus"] = "D_s^-" 

        # mesons
        particle_dict["pion"] = "pi^-"
        particle_dict["pion-zero"] = "pi_0"
        particle_dict["anti-pion"]="pi^+"   

        return particle_dict
    
