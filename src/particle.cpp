#include "particle.h"

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
  -2122,
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

std::map<std::string,int> pyne::particle::altnames;
std::map<int,std::string> pyne::particle::id_name;
std::map<std::string,int> pyne::particle::name_id;
std::map<std::string,std::string> pyne::particle::docs;

void * pyne::particle::_fill_maps() {
  using std::make_pair;

  std::string pyne::particle::_docs[NUM_PARTICLES] = {
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
};

// name functions
std::string pyne::particle::name(char *s) {
  return pyne::particle::name(std::string(s));
};

std::string pyne::particle::name(std::string s) {
  // check std name
  if(0 < names.count(s))
    return s;
  // check alternative name
  if(0 < altnames.count(s))
    return id_name[altnames[s]];
  // check if is a heavy ion
  if( pyne::nucname::isnuclide(s) )
    return s
  else
    throw NotAParticle(s,"???");
};

std::string pyne::particle::describe(char *s) {
  return pyne::particle::describe(std::string(s));
};

std::string pyne::particle::describe(std::string s) {
  // check std name
  if(0 < names.count(s))
    return docs[s];
  // check alternative name
  if(0 < altnames.count(s))
    return docs[altnames[s]];
  // check if is a heavy ion
  if( pyne::nucname::isnuclide(s) )
    return s
  else
    throw NotAParticle(s,"???");
};


// functions
bool pyne::particle::is_anti(int pdc_num) {
  if ( pdc_num < 0 )
    return true;
  else if ( pdc_num > 0 )
    return false;
  else
    throw NotAPdcNum(int pdcnum);
}


bool pyne::particle::is_boson(int pdc_num) {
  if ( pdc_num >= 21 && pdc_num <= 37 )
    return true;
  else
    return false;
}

bool pyne::particle::is_lepton(int pdc_num) {
  if ( pdc_num >= 11 && pdc_num <= 18 )
    return true;
  else
    return false;
}

