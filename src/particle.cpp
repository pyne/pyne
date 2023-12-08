#ifndef PYNE_IS_AMALGAMATED
#include "particle.h"
#endif

std::string pyne::particle::_names[NUM_PARTICLES] = {
  // leptons
  "Electron",
  "Positron",
  "ElectronNeutrino",
  "ElectronAntiNeutrino",
  "Muon",
  "AntiMuon",
  "MuonNeutrino",
  "MuonAntiNeutrino",
  "Tauon",
  "AntiTauon",
  "TauNeutrino",
  "TauAntiNeutrino",
  // gauge bosons
  "Photon",
  // light mesons
  "Pion",
  "AntiPion",
  // strange mesons
  "Kaon",
  "AntiKaon",
  "KaonZeroShort",
  "KaonZero",
  "AntiKaonZero",
  // light baryons
  "Neutron",
  "AntiNeutron",
  "Proton",
  "AntiProton",
  // strange baryons
  "Lambda",
  "AntiLambda",
  "Sigma-",
  "AntiSigma-",
  "Sigma+",
  "AntiSigma+",
  "Sigma",
  "AntiSigmaZero"
  // Charmed baryons
};

int pyne::particle::_pdcids[NUM_PARTICLES] = {
  11,
  -11,
  12,
  -12,
  13,
  -13,
  14,
  -14,
  15,
  -15,
  16,
  -16,
  // gauge bosons
  22,
  // light mesons
  211,
  -211,
  // strange mesons
  321,
  -321,
  310,
  311,
  -311,
  // light baryons
  2112,
  -2112,
  2212,
  -2212,
  // strange Baryons
  3122,
  -3122,
  3112,
  3112,
  3222,
  -3222,
  3212,
  -3212
  // charmed baryons
};

std::set<std::string> pyne::particle::names(pyne::particle::_names,
              pyne::particle::_names+NUM_PARTICLES);

std::set<int> pyne::particle::pdc_nums(pyne::particle::_pdcids,
                 pyne::particle::_pdcids+NUM_PARTICLES);

std::map<std::string,int> pyne::particle::altnames;
std::map<int,std::string> pyne::particle::id_name;
std::map<std::string,int> pyne::particle::name_id;
std::map<std::string,std::string> pyne::particle::docs;

std::map<std::string,std::string> pyne::particle::part_to_fluka;
std::map<std::string,std::string> pyne::particle::part_to_mcnp;
std::map<std::string,std::string> pyne::particle::part_to_mcnp6;
std::map<std::string,std::string> pyne::particle::part_to_geant4;


void * pyne::particle::_fill_maps() {
  using std::make_pair;

  std::string _docs[NUM_PARTICLES] = {
    // leptons
    "Electron",
    "Positron",
    "Electron Neutrino",
    "Electron Anti Neutrino",
    "Muon Neutrino",
    "Anti Muon",
    "Muon Neutrino",
    "Muon Anti Neutrino",
    "Tauon",
    "Anti Tauon",
    "Tau Neutrino",
    "Tau Anti Neutrino",
    // gauge bosons
    "Photon",
    // light mesons
    "Pion",
    "Anti Pion",
    // strange mesons
    "Kaon",
    "Anti Kaon",
    "Kaon Zero Short",
    "Kaon Zero",
    "Anti Kaon Zero",
    // light baryons
    "Neutron",
    "Anti Neutron",
    "Proton",
    "Anti Proton",
    // strange baryons
    "Lambda",
    "Anti Lambda",
    "Sigma-",
    "Anti Sigma-",
    "Sigma+",
    "Anti Sigma+",
    "Sigma",
    "Anti Sigma Zero"
    // Charmed baryons
  };

  int pid;  // particle id
  for ( int i = 0 ; i < NUM_PARTICLES ; i++ ) {
    pid = _pdcids[i];
    // make id to name map
    id_name[pid] = _names[i];
    // make name to id map
    name_id[_names[i]] = pid;
    // make doc correspondence
    docs[_names[i]] = _docs[i];
  }

  // make the alternates
  altnames["Hydrogen"] = name_id["Proton"];
  altnames["Protium"] = name_id["Proton"];
  altnames["Beta"] = name_id["Electron"];
  altnames["Beta-"] = name_id["Electron"];
  altnames["Beta+"] = name_id["Positron"];
  altnames["Gamma"] = name_id["Photon"];
  altnames["X-Ray"] = name_id["Photon"];

  part_to_mcnp["Neutron"]="n";
  part_to_mcnp["Photon"]="p";
  part_to_mcnp["Electron"]="e";

  part_to_mcnp6["Neutron"]="n";
  part_to_mcnp6["Photon"]="p";
  part_to_mcnp6["Electron"]="e";
  part_to_mcnp6["Proton"]="h";

  part_to_fluka["Electron"]="ELECTRON";
  part_to_fluka["Positron"]="POSITRON";
  part_to_fluka["ElectronNeutrino"] ="NEUTRIE";
  part_to_fluka["ElectronAntiNeutrino"] ="ANEUTRIE";
  part_to_fluka["Muon"]="MUON+";
  part_to_fluka["AntiMuon"]="MUON-";
  part_to_fluka["MuonNeutrino"]="NEUTRIM";
  part_to_fluka["MuonAntiNeutrino"]="ANEUTRIM";
  part_to_fluka["Tauon"]="TAU+";
  part_to_fluka["Anti Tauon"]="TAU-";
  part_to_fluka["TauNeutrino"]="NEUTRIT";
  part_to_fluka["TauAntiNeutrino"]="ANEUTRIT";
  // gauge bosons
  part_to_fluka["Photon"]="PHOTON";
  // light mesons
  part_to_fluka["Pion"]="PION-";
  part_to_fluka["Anti Pion"]="PION+";
  // strange mesons
  part_to_fluka["Kaon"]="KAON+";
  part_to_fluka["AntiKaon"]="KAON-";
  part_to_fluka["KaonZero Short"]="KAONSHRT";
  part_to_fluka["KaonZero"]="KAONZERO";
  part_to_fluka["AntiKaonZero"]="AKAONZER";
  // light baryons
  part_to_fluka["Neutron"]="NEUTRON";
  part_to_fluka["AntiNeutron"]="ANEUTRON";
  part_to_fluka["Proton"]="PROTON";
  part_to_fluka["AntiProton"]="APROTON";
  // strange baryons
  part_to_fluka["Lambda"]="LAMBDA";
  part_to_fluka["AntiLambda"]="ALAMBDA";
  part_to_fluka["Sigma-"]="SIGMA-";
  part_to_fluka["Anti Sigma-"]="ASIGMA-";
  part_to_fluka["Sigma+"]="SIGMA+";
  part_to_fluka["Anti Sigma+"]="ASIGMA+";
  part_to_fluka["Sigma"]="SIGMAZER";
  part_to_fluka["AntiSigmaZero"]="ASIGMAZE";

  part_to_geant4["Electron"]="e-";
  part_to_geant4["Positron"]="e+";
  part_to_geant4["ElectronNeutrino"] ="nu_e";
  part_to_geant4["ElectronAntiNeutrino"] ="anti_nu_e";
  part_to_geant4["Muon"]="mu+";
  part_to_geant4["AntiMuon"]="mu-";
  part_to_geant4["MuonNeutrino"]="nu_mu";
  part_to_geant4["MuonAntiNeutrino"]="anti_nu_mu";
  part_to_geant4["Tauon"]="tau+";
  part_to_geant4["Anti Tauon"]="tau-";
  part_to_geant4["TauNeutrino"]="nu_tau";
  part_to_geant4["TauAntiNeutrino"]="anti_nu_tau";
  // gauge bosons
  part_to_geant4["Photon"]="gamma";
  // light mesons
  part_to_geant4["Pion"]="pi-";
  part_to_geant4["Anti Pion"]="pi+";
  // strange mesons
  part_to_geant4["Kaon"]="kaon+";
  part_to_geant4["AntiKaon"]="kaon-";
  part_to_geant4["KaonZero Short"]="kaon0S";
  part_to_geant4["KaonZero"]="kaon0";
  // light baryons
  part_to_geant4["Neutron"]="neutron";
  part_to_geant4["AntiNeutron"]="anti_neutron";
  part_to_geant4["Proton"]="proton";
  part_to_geant4["AntiProton"]="anti_proton";
  // strange baryons
  part_to_geant4["Lambda"]="lambda";
  part_to_geant4["AntiLambda"]="anti_lambda";
  part_to_geant4["Sigma-"]="sigma-";
  part_to_geant4["Anti Sigma-"]="anti_sigma-";
  part_to_geant4["Sigma+"]="sigma+";
  part_to_geant4["Anti Sigma+"]="anti_sigma+";
  part_to_geant4["Sigma"]="sigma0";
  part_to_geant4["AntiSigmaZero"]="anti_sigma0";
  return NULL;
}

void * pyne::particle::filler = pyne::particle::_fill_maps();

// is hydrogen
bool pyne::particle::is_hydrogen(int s) {
  if(s == name_id["Proton"])
    return true;
  if(pyne::particle::is_hydrogen(pyne::nucname::name(s)))
      return true;
  return false;
}

bool pyne::particle::is_hydrogen(char *s) {
  return pyne::particle::is_hydrogen(std::string(s));
}

bool pyne::particle::is_hydrogen(std::string s) {
  // check std name
  if(name_id[s] == name_id["Proton"])
    return true;
  if(altnames[s] == name_id["Proton"])
    return true;
  if(pyne::nucname::name(s).find("H1") != std::string::npos)
    return true;
  return false;
}
// heavy ion
bool pyne::particle::is_heavy_ion(int s) {
  return pyne::particle::is_heavy_ion(std::string(id_name[s]));
}

bool pyne::particle::is_heavy_ion(char *s) {
  return pyne::particle::is_heavy_ion(std::string(s));
}

bool pyne::particle::is_heavy_ion(std::string s) {
  if(pyne::nucname::isnuclide(s)) {
    if(pyne::particle::is_hydrogen(s))
      return false;
    else
      return true;
  }
  return false;
}

// is valid functions
bool pyne::particle::is_valid(int s) {
  if(pyne::nucname::isnuclide(s))
    return true;
  else
    return pyne::particle::is_valid(std::string(id_name[s]));
}

bool pyne::particle::is_valid(char *s) {
  return pyne::particle::is_valid(std::string(s));
}

bool pyne::particle::is_valid(std::string s) {
  // check std name
  if(0 < names.count(s))
    return true;
  // check alternative name
  if(0 < altnames.count(s))
    return true;
  // check if is a heavy ion
  if(pyne::nucname::isnuclide(s))
    return true;
  else
    return false;
}

// pdc functions
int pyne::particle::id(int s) {
  if (0 < pdc_nums.count(s))
    return s;
  else
    return 0;
}

int pyne::particle::id(char *s) {
  return pyne::particle::id(std::string(s));
}

int pyne::particle::id(std::string s) {
  if(pyne::nucname::isnuclide(s))
    {
      if(pyne::particle::is_hydrogen(s))
  return name_id["Proton"];
      if( pyne::particle::is_heavy_ion(s) )
  return 0;
    }

  if (0 < pdc_nums.count(name_id[s]))
    return name_id[s];
  if (0 < pdc_nums.count(altnames[s]))
    return altnames[s];
  return 0;
}

// name functions
std::string pyne::particle::name(int s) {
  if(s < 9999999)
    return pyne::particle::name(id_name[s]);
  if(pyne::nucname::isnuclide(s))
    return pyne::particle::name(pyne::nucname::name(s));
  return pyne::particle::name(id_name[s]);
}

std::string pyne::particle::name(char *s) {
  return pyne::particle::name(std::string(s));
}

std::string pyne::particle::name(std::string s) {
  // check if is a hydrogen
  if(pyne::nucname::isnuclide(s))
    {
      if(pyne::particle::is_hydrogen(s))
        return "Proton";
      if( pyne::particle::is_heavy_ion(s) )
        return s;
    }
  // check std name
  if(0 < names.count(s))
    return s;
  // check alternative name
  if(0 < altnames.count(s))
    return id_name[altnames[s]];
  // check for heavy ion
  else
    throw NotAParticle(s);
}

// convert name to mcnp id
std::string pyne::particle::mcnp(int s) {
  return pyne::particle::mcnp(pyne::particle::name(s));
}

std::string pyne::particle::mcnp(char *s) {
  return pyne::particle::mcnp(pyne::particle::name(s));
}

std::string pyne::particle::mcnp(std::string s) {
  if(0 < part_to_mcnp.count(pyne::particle::name(s)))
    return part_to_mcnp[pyne::particle::name(s)];
  else
    {
      std::cout << "Not a valid MCNP5 particle" << std::endl;
      return "?";
    }
}

// convert name to mcnp6 id
std::string pyne::particle::mcnp6(int s) {
  return pyne::particle::mcnp6(pyne::particle::name(s));
}

std::string pyne::particle::mcnp6(char *s) {
  return pyne::particle::mcnp6(pyne::particle::name(s));
}

std::string pyne::particle::mcnp6(std::string s) {
  if(0 < part_to_mcnp6.count(pyne::particle::name(s)))
    return part_to_mcnp6[pyne::particle::name(s)];
  else
      return "?";
}

// convert name to fluka id
std::string pyne::particle::fluka(int s) {
  return pyne::particle::fluka(pyne::particle::name(s));
}

std::string pyne::particle::fluka(char *s) {
  return pyne::particle::fluka(pyne::particle::name(s));
}

std::string pyne::particle::fluka(std::string s) {
  if (pyne::particle::is_heavy_ion(s))
    return "HEAVYION";
  else if(0 < part_to_fluka.count(pyne::particle::name(s)))
    return part_to_fluka[pyne::particle::name(s)];
  else
    {
      std::cout << "Not a valid Fluka particle" << std::endl;
      return "???????";
    }
}

// convert name to geant4 id
std::string pyne::particle::geant4(int s) {
  return pyne::particle::geant4(pyne::particle::name(s));
}

std::string pyne::particle::geant4(char *s) {
  return pyne::particle::geant4(pyne::particle::name(s));
}

std::string pyne::particle::geant4(std::string s) {
  if (pyne::particle::is_heavy_ion(s))
    return "GenericIon";
  else if(0 < part_to_geant4.count(pyne::particle::name(s)))
    return part_to_geant4[pyne::particle::name(s)];
  else
    {
      std::cout << "Not a valid Geant4 particle" << std::endl;
      return "???????";
    }
}


// describe functions
std::string pyne::particle::describe(int s) {
  if(pyne::nucname::isnuclide(s))
    return pyne::particle::describe(pyne::nucname::name(s));
  return pyne::particle::describe(id_name[s]);
}

std::string pyne::particle::describe(char *s) {
  return pyne::particle::describe(std::string(s));
}

std::string pyne::particle::describe(std::string s) {
  // check if is a hydrogen
  if (pyne::nucname::isnuclide(s))
    {
      if (pyne::particle::is_hydrogen(s))
  return docs[pyne::particle::name(s)];
      if (pyne::particle::is_heavy_ion(s))
  return "Is a heavy ion";
    }
  // check std name
  if(0 < names.count(s))
    return docs[s];
  // check alternative name
  if(0 < altnames.count(s))
    return docs[id_name[altnames[s]]];
  // check if is a heavy ion
  else
    throw NotAParticle(s);
}
