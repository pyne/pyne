// Converts between naming conventions for nuclides.
// zzaaam is for numerals only (923350).
// name is for letters  as well (U-235).
// MCNP is for numerals without the meta-stable flag (92235), as used in MCNP.

#include "nucname.h"

/*** Constructs the LL to zz Dictionary ***/
pyne::nucname::name_zz_t pyne::nucname::get_name_zz()
{
  pyne::nucname::name_zz_t lzd;

  lzd["BE"] = 04;
  lzd["BA"] = 56;
  lzd["BH"] = 107;
  lzd["BI"] = 83;
  lzd["BK"] = 97;
  lzd["BR"] = 35;
  lzd["RU"] = 44;
  lzd["RE"] = 75;
  lzd["RF"] = 104;
  lzd["RG"] = 111;
  lzd["RA"] = 88;
  lzd["RB"] = 37;
  lzd["RN"] = 86;
  lzd["RH"] = 45;
  lzd["TM"] = 69;
  lzd["H"] = 01;
  lzd["P"] = 15;	
  lzd["GE"] = 32;
  lzd["GD"] = 64;
  lzd["GA"] = 31;
  lzd["OS"] = 76;
  lzd["HS"] = 108;
  lzd["ZN"] = 30;
  lzd["HO"] = 67;
  lzd["HF"] = 72;
  lzd["HG"] = 80;
  lzd["HE"] = 02;
  lzd["PR"] = 59;
  lzd["PT"] = 78;
  lzd["PU"] = 94;
  lzd["PB"] = 82;
  lzd["PA"] = 91;
  lzd["PD"] = 46;
  lzd["PO"] = 84;
  lzd["PM"] = 61;
  lzd["C"] = 6;
  lzd["K"] = 19;
  lzd["O"] = 8;
  lzd["S"] = 16;
  lzd["W"] = 74;
  lzd["EU"] = 63;
  lzd["ES"] = 99;
  lzd["ER"] = 68;
  lzd["MD"] = 101;
  lzd["MG"] = 12;
  lzd["MO"] = 42;
  lzd["MN"] = 25;
  lzd["MT"] = 109;
  lzd["U"] = 92;
  lzd["FR"] = 87;
  lzd["FE"] = 26;
  lzd["FM"] = 100;
  lzd["NI"] = 28;
  lzd["NO"] = 102;
  lzd["NA"] = 11;
  lzd["NB"] = 41;
  lzd["ND"] = 60;
  lzd["NE"] = 10;
  lzd["ZR"] = 40;
  lzd["NP"] = 93;
  lzd["B"] = 05;
  lzd["CO"] = 27;
  lzd["CM"] = 96;
  lzd["F"] = 9;
  lzd["CA"] = 20;
  lzd["CF"] = 98;
  lzd["CE"] = 58;
  lzd["CD"] = 48;
  lzd["V"] = 23;
  lzd["CS"] = 55;
  lzd["CR"] = 24;
  lzd["CU"] = 29;
  lzd["SR"] = 38;
  lzd["KR"] = 36;
  lzd["SI"] = 14;
  lzd["SN"] = 50;
  lzd["SM"] = 62;
  lzd["SC"] = 21;
  lzd["SB"] = 51;
  lzd["SG"] = 106;
  lzd["SE"] = 34;
  lzd["YB"] = 70;
  lzd["DB"] = 105;
  lzd["DY"] = 66;
  lzd["DS"] = 110;
  lzd["LA"] = 57;
  lzd["CL"] = 17;
  lzd["LI"] = 03;
  lzd["TL"] = 81;
  lzd["LU"] = 71;
  lzd["LR"] = 103;
  lzd["TH"] = 90;
  lzd["TI"] = 22;
  lzd["TE"] = 52;
  lzd["TB"] = 65;
  lzd["TC"] = 43;
  lzd["TA"] = 73;
  lzd["AC"] = 89;
  lzd["AG"] = 47;
  lzd["I"] = 53;
  lzd["IR"] = 77;
  lzd["AM"] = 95;
  lzd["AL"] = 13;
  lzd["AS"] = 33;
  lzd["AR"] = 18;
  lzd["AU"] = 79;
  lzd["AT"] = 85;
  lzd["IN"] = 49;
  lzd["Y"] = 39;
  lzd["N"] = 07;
  lzd["XE"] = 54;
  lzd["CN"] = 112;
  lzd["FL"] = 114;
  lzd["LV"] = 116;

  return lzd;
};
pyne::nucname::name_zz_t pyne::nucname::name_zz = pyne::nucname::get_name_zz();


/*** Constructs zz to LL dictionary **/
pyne::nucname::zzname_t pyne::nucname::get_zz_name()
{
  zzname_t zld;
  for (name_zz_iter i = name_zz.begin(); i != name_zz.end(); i++)
  {
    zld[i->second] = i->first;
  }
  return zld;
};
pyne::nucname::zzname_t pyne::nucname::zz_name = pyne::nucname::get_zz_name();



/******************************************/
/*** Define useful elemental group sets ***/
/******************************************/

pyne::nucname::zz_group pyne::nucname::name_to_zz_group(pyne::nucname::name_group eg)
{
  zz_group zg;
  for (name_group_iter i = eg.begin(); i != eg.end(); i++)
  {
    zg.insert(name_zz[*i]);
  }
  return zg;
};

// Lanthanides
pyne::nucname::name_t pyne::nucname::LAN_array[15] = {"LA", "CE", "PR", "ND", "PM", "SM", "EU", \
                                        "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU"};
pyne::nucname::name_group pyne::nucname::LAN (pyne::nucname::LAN_array, pyne::nucname::LAN_array+15);
pyne::nucname::zz_group pyne::nucname::lan = pyne::nucname::name_to_zz_group(pyne::nucname::LAN);

// Actinides
pyne::nucname::name_t pyne::nucname::ACT_array[15] = {"AC", "TH", "PA", "U",  "NP", "PU", "AM", "CM", \
                                        "BK", "CF", "ES", "FM", "MD", "NO", "LR"};
pyne::nucname::name_group pyne::nucname::ACT (pyne::nucname::ACT_array, pyne::nucname::ACT_array+15);
pyne::nucname::zz_group pyne::nucname::act = pyne::nucname::name_to_zz_group(pyne::nucname::ACT);

// Transuarnics
pyne::nucname::name_t pyne::nucname::TRU_array[22] = {"NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", \
                                        "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", \
                                        "MT", "DS", "RG", "CN", "FL", "LV"};
pyne::nucname::name_group pyne::nucname::TRU (pyne::nucname::TRU_array, pyne::nucname::TRU_array+22);
pyne::nucname::zz_group pyne::nucname::tru = pyne::nucname::name_to_zz_group(pyne::nucname::TRU);

//Minor Actinides
pyne::nucname::name_t pyne::nucname::MA_array[10] = {"NP", "AM", "CM", "BK", "CF", "ES", "FM", "MD", \
                                       "NO", "LR"};
pyne::nucname::name_group pyne::nucname::MA (pyne::nucname::MA_array, pyne::nucname::MA_array+10);
pyne::nucname::zz_group pyne::nucname::ma = pyne::nucname::name_to_zz_group(pyne::nucname::MA);

//Fission Products
pyne::nucname::name_t pyne::nucname::FP_array[88] = {"AG", "AL", "AR", "AS", "AT", "AU", "B",  "BA", \
                                       "BE", "BI", "BR", "C",  "CA", "CD", "CE", "CL", \
                                       "CO", "CR", "CS", "CU", "DY", "ER", "EU", "F",  \
                                       "FE", "FR", "GA", "GD", "GE", "H",  "HE", "HF", \
                                       "HG", "HO", "I",  "IN", "IR", "K",  "KR", "LA", \
                                       "LI", "LU", "MG", "MN", "MO", "N",  "NA", "NB", \
                                       "ND", "NE", "NI", "O",  "OS", "P",  "PB", "PD", \
                                       "PM", "PO", "PR", "PT", "RA", "RB", "RE", "RH", \
                                       "RN", "RU", "S",  "SB", "SC", "SE", "SI", "SM", \
                                       "SN", "SR", "TA", "TB", "TC", "TE", "TI", "TL", \
                                       "TM", "V",  "W",  "XE", "Y",  "YB", "ZN", "ZR"};
pyne::nucname::name_group pyne::nucname::FP (pyne::nucname::FP_array, pyne::nucname::FP_array+88);
pyne::nucname::zz_group pyne::nucname::fp = pyne::nucname::name_to_zz_group(pyne::nucname::FP);


/***************************/
/*** isnuclide functions ***/
/***************************/

bool pyne::nucname::isnuclide(std::string nuc)
{
  int n;
  try {
    n = id(nuc);
  }
  catch(NotANuclide) {
    return false;
  };
  return isnuclide(n);
};

bool pyne::nucname::isnuclide(char * nuc)
{
  return isnuclide(std::string(nuc));
};

bool pyne::nucname::isnuclide(int nuc)
{
  int n;
  try {
    n = id(nuc);
  }
  catch(NotANuclide) {
    return false;
  };
  int zzz = n / 10000000;
  int aaa = (n % 10000000) / 10000;
  if (aaa < zzz)
    return false;
  return true;
};



/********************/
/*** id functions ***/
/********************/
int pyne::nucname::id(int nuc)
{
  int newnuc;
  int zzz = nuc / 10000000;     // ZZZ ?
  int aaassss = nuc % 10000000; // AAA-SSSS ?
  int aaa = aaassss / 10000;    // AAA ?
  //int ssss = aaassss % 10000;   // SSSS ?
  // Nuclide must already be in id form
  if (zzz <= aaa && aaa <= zzz * 7)
  {
    // Normal nuclide
    return nuc;
  }
  else if (aaassss == 0 && 0 < zz_name.count(zzz))
  {
    // Natural elemental nuclide:  ie for Urnaium = 920000000
    return nuc;
  }

  // Not in id form, try  ZZZAAAM form.
  zzz = nuc / 10000;     // ZZZ ?
  aaassss = nuc % 10000; // AAA-SSSS ?
  aaa = aaassss / 10;    // AAA ?
  if (zzz <= aaa && aaa <= zzz * 7)
  {
    // ZZZAAAM nuclide
    return (zzz*10000000) + (aaa*10000) + (nuc%10);
  }
  else if (aaa < zzz && 0 < zz_name.count(aaa))
  {
    // Cinder-form (aaazzzm), ie 2350920
    return (aaa*10000000) + (zzz*10000000) + (nuc%10);
  }
  //else if (aaassss == 0 && 0 == zz_name.count(nuc/1000) && 0 < zz_name.count(zzz))
  else if (aaassss == 0 && 0 < zz_name.count(zzz))
  {
    // zzaaam form natural nuclide
    return zzz * 10000000;
  }

  if (nuc >= 1000000){
    // From now we assume no metastable info has been given.
    throw IndeterminateNuclideForm(nuc, "");
  };

  // Nuclide is not in zzaaam form, 
  // Try MCNP form, ie zzaaa
  zzz = nuc / 1000;
  aaa = nuc % 1000; 
  if (zzz <= aaa)
  {
    if (aaa - 400 < 0)
    {
      if (nuc == 95242)
        return nuc * 10000 + 1;  // special case MCNP Am-242m
      else
        return nuc * 10000;  // Nuclide in normal MCNP form
    }
    else
    {
      // Nuclide in MCNP metastable form
      if (nuc == 95642)
        return (95642 - 400)*10000;  // special case MCNP Am-242
      nuc = ((nuc - 400) * 10000) + 1;
      while (3.0 < (float ((nuc/10000)%1000) / float (nuc/10000000)))
        nuc -= 999999;
      return nuc;
    }
  }
  else if (aaa == 0 && 0 < zz_name.count(zzz))
    // MCNP form natural nuclide
    return zzz * 10000000;

  // Not a normal nuclide, might be a 
  // Natural elemental nuclide.  
  // ie 92 for Urnaium = 920000
  if (0 < zz_name.count(nuc))
    return nuc * 10000000;
  throw IndeterminateNuclideForm(nuc, "");
};

int pyne::nucname::id(char * nuc)
{
  std::string newnuc (nuc);
  return zzaaam(newnuc);
};

int pyne::nucname::id(std::string nuc)
{
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int newnuc;

  // Get the string into a regular form
  std::string nucstr = pyne::to_upper(nuc);
  nucstr = pyne::remove_substring(nucstr, "-");
  int nuclen = nucstr.length();

  if (pyne::contains_substring(pyne::digits, nucstr.substr(0, 1)))
  {
    if (pyne::contains_substring(pyne::digits, nucstr.substr(nuclen-1, nuclen)))
    {
      // Nuclide must actually be an integer that 
      // just happens to be living in string form.
      newnuc = pyne::to_int(nucstr);
      newnuc = id(newnuc);
    }
    else
    {
      // probably in NIST-like form (242Am)
      // Here we know we have both digits and letters
      std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);
      newnuc = pyne::to_int(anum_str) * 10000;

      // Add the Z-number
      std::string elem_name = pyne::remove_characters(nucstr, pyne::digits);
      if (0 < name_zz.count(elem_name))
        newnuc = (10000000 * name_zz[elem_name]) + newnuc;
      else
        throw NotANuclide(nucstr, newnuc);
    };
  }
  else if (pyne::contains_substring(pyne::alphabet, nucstr.substr(0, 1)))
  {
    // Nuclide is probably in name form, or some variation therein
    std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);

    // natural element form, a la 'U' -> 920000000
    if (anum_str.empty() && (0 < name_zz.count(nucstr)))
      return 10000000 * name_zz[nucstr]; 

    int anum = pyne::to_int(anum_str);

    // bad form
    if (anum < 0)
      throw NotANuclide(nucstr, anum); 

    // Figure out if we are meta-stable or not
    std::string end_char = pyne::last_char(nucstr);
    if (end_char == "M")
      newnuc = (10000 * anum) + 1;
    else if (pyne::contains_substring(pyne::digits, end_char))
      newnuc = (10000 * anum);
    else
      throw NotANuclide(nucstr, newnuc);

    // Add the Z-number
    std::string elem_name = pyne::remove_characters(nucstr.substr(0, nuclen-1), 
                                                    pyne::digits);
    if (0 < name_zz.count(elem_name))
      newnuc = (10000000 * name_zz[elem_name]) + newnuc;
    else
      throw NotANuclide(nucstr, newnuc);
  }
  else
  {
    // Clearly not a nuclide
    throw NotANuclide(nuc, nucstr);
  }
  return newnuc;  
};



/**********************/
/*** name functions ***/
/**********************/
std::string pyne::nucname::name(int nuc)
{
  int nucid = id(nuc);
  std::string newnuc = "";

  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add LL
  newnuc += zz_name[zzz];

  // Add A-number
  if (0 < aaa)
    newnuc += pyne::to_str(aaa);

  // Add meta-stable flag
  if (0 < ssss)
    newnuc += "M";

  return newnuc;
};



std::string pyne::nucname::name(char * nuc)
{
  std::string newnuc (nuc);
  return name(newnuc);
}


std::string pyne::nucname::name(std::string nuc)
{
  return name(id(nuc));
}


/************************/
/*** zzaaam functions ***/
/************************/
int pyne::nucname::zzaaam(int nuc)
{
  int nucid = id(nuc);
  int zzzaaa = nucid / 10000;
  int ssss = nucid % 10000;
  if (10 <= ssss)
    ssss = 9;
  return zzzaaa*10 + ssss;
};


int pyne::nucname::zzaaam(char * nuc)
{
  std::string newnuc (nuc);
  return zzaaam(newnuc);
};


int pyne::nucname::zzaaam(std::string nuc)
{
  return zzaaam(id(nuc));
};


/**********************/
/*** mcnp functions ***/
/**********************/
int pyne::nucname::mcnp(int nuc)
{
  nuc = id(nuc);
  int ssss = nuc % 10000;
  int newnuc = nuc / 10000;

  // special case Am242(m)
  if (newnuc == 95242 && ssss < 2)
    ssss = (ssss + 1) % 2;

  // Handle the crazy MCNP meta-stable format
  if (0 != ssss && ssss < 10) 
    newnuc += 300 + (ssss * 100);

  return newnuc;
};



int pyne::nucname::mcnp(char * nuc)
{
  std::string newnuc (nuc);
  return mcnp(newnuc);
};



int pyne::nucname::mcnp(std::string nuc)
{
  return mcnp(id(nuc));
};





/*************************/
/*** serpent functions ***/
/*************************/
std::string pyne::nucname::serpent(int nuc)
{
  int nucid = id(nuc);
  std::string newnuc = "";

  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int zzz = nucid / 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add LL
  std::string llupper = pyne::to_upper(zz_name[zzz]);
  std::string lllower = pyne::to_lower(zz_name[zzz]);
  newnuc += llupper[0];
  for (int l = 1; l < lllower.size(); l++)
    newnuc += lllower[l];  

  // Add required dash
  newnuc += "-";

  // Add A-number
  if (0 < aaassss)
    newnuc += pyne::to_str(aaa);
  else if (0 == aaassss)
    newnuc += "nat";

  // Add meta-stable flag
  if (0 < ssss)
    newnuc += "m";

  return newnuc;
};


std::string pyne::nucname::serpent(char * nuc)
{
  std::string newnuc (nuc);
  return serpent(newnuc);
};


std::string pyne::nucname::serpent(std::string nuc)
{
  return serpent(id(newnuc));
};




/**********************/
/*** nist functions ***/
/**********************/
std::string pyne::nucname::nist(int nuc)
{
  int nucid = id(nuc);
  std::string newnuc = "";

  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add A-number
  if (0 < aaassss)
    newnuc += pyne::to_str(aaa);

  // Add name
  std::string name_upper = pyne::to_upper(zz_name[zzz]);
  std::string name_lower = pyne::to_lower(zz_name[zzz]);
  newnuc += name_upper[0];
  for (int l = 1; l < name_lower.size(); l++)
    newnuc += name_lower[l];  

  // Add meta-stable flag
  // No metastable flag for NIST, 
  // but could add star, by uncommenting below
  //if (0 < mod_10)
  //  newnuc += "*";

  return newnuc;
};


std::string pyne::nucname::nist(char * nuc)
{
  std::string newnuc (nuc);
  return nist(newnuc);
};


std::string pyne::nucname::nist(std::string nuc)
{
  return nist(id(nuc));
};




/************************/
/*** cinder functions ***/
/************************/
int pyne::nucname::cinder(int nuc)
{
  // cinder nuclides of form aaazzzm
  int nucid = id(nuc);
  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;
  if (10 <= ssss)
    ssss = 9;
  return (aaa*10000) + (zzz*10) + ssss;
};



int pyne::nucname::cinder(char * nuc)
{
  std::string newnuc (nuc);
  return cinder(newnuc);
};



int pyne::nucname::cinder(std::string nuc)
{
  return cinder(id(nuc));
};


/**********************/
/*** ALARA functions ***/
/**********************/
std::string pyne::nucname::alara(int nuc)
{
  int nucid = id(nuc);
  std::string newnuc = "";
  std::string ll = "";

  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucint);

  // Add LL, in lower case
  ll += zz_name[zzz];

  for(int i = 0; ll[i] != '\0'; i++)
    ll[i] = tolower(ll[i]);
  newnuc += ll;

  // Add A-number
  if (0 < mod_10000){
    newnuc += ":";
    newnuc += pyne::to_str(aaa);
  }

  // Note, ALARA input format does not use metastable flag
  return newnuc;
};


std::string pyne::nucname::alara(char * nuc)
{
  std::string newnuc (nuc);
  return alara(newnuc);
}


std::string pyne::nucname::alara(std::string nuc)
{
  return alara(id(nuc));
}
