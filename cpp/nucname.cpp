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
pyne::nucname::name_t pyne::nucname::TRU_array[19] = {"NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", \
                                        "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", \
                                        "MT", "DS", "RG"};
pyne::nucname::name_group pyne::nucname::TRU (pyne::nucname::TRU_array, pyne::nucname::TRU_array+19);
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



/********************/
/*** Current Form ***/
/********************/
std::string pyne::nucname::current_form(std::string nuc)
{
  // returns current form of a nuclide.
  using namespace pyne;

  nuc = to_upper(nuc);
  nuc = remove_substring(nuc, "-");

  if ( contains_substring(alphabet, nuc.substr(0,1)) )
    return "name";
  else
  {
    if (nuc.length() == 7)
      return "zzaaam";
    else if ( contains_substring("23456789", last_char(nuc)) )
      return "MCNP";

    if ( ternary_ge( to_int(nuc.substr(0,nuc.length()-4)), to_int(slice_from_end(nuc,-4,3)), to_int(nuc.substr(0,nuc.length()-4)) * 5) )
      return "zzaaam";
    else if ( ternary_ge( to_int(nuc.substr(0,nuc.length()-3)), to_int(slice_from_end(nuc,-3,3)), to_int(nuc.substr(0,nuc.length()-3)) * 5) )
      return "MCNP";
    else
      throw IndeterminateNuclideForm(nuc, "");
  };
};

std::string pyne::nucname::current_form(int nuc)
{
   return pyne::nucname::current_form(pyne::to_str(nuc));
};




/************************/
/*** zzaaam functions ***/
/************************/
int pyne::nucname::zzaaam(int nuc)
{
  int newnuc;

  int mod_10000 = nuc % 10000;
  int div_10000 = nuc / 10000;
  int mod_10000_div_10 = mod_10000 / 10;

  // Nuclide must already be in zzaaam form
  if (div_10000 <= mod_10000_div_10 && mod_10000_div_10 <= div_10000 * 7)
  {
    // Normal nuclide
    newnuc = nuc;
    return newnuc;
  }
  else if (mod_10000 == 0 && 0 < zz_name.count(div_10000))
  {
    // Natural elemental nuclide:  ie for Urnaium = 920000
    newnuc = nuc;
    return newnuc;
  }
  else if (mod_10000_div_10 < div_10000 && 0 < zz_name.count(mod_10000_div_10))
  {
    // Cinder-form (aaazzzm), ie 2350920
    newnuc = (mod_10000_div_10*10000) + (div_10000*10) + (nuc%10);
    return newnuc;
  }

  // Nuclide is not in zzaaam form, 
  // Try MCNP form, ie zzaaa
  int mod_1000 = nuc % 1000; 
  int div_1000 = nuc / 1000;

  if (div_1000 <= mod_1000)
  {
    if (mod_1000 - 400 < 0)
    {
      // Nuclide in normal MCNP form
      newnuc = nuc * 10;

      // special case MCNP Am-242m
      if (newnuc == 952420)
        newnuc += 1;
    }
    else
    {
      // Nuclide in MCNP metastable form
      newnuc = ((nuc - 400) * 10) + 1;
      while (3.0 < (float ((newnuc/10)%1000) / float (newnuc/10000)))
      {
        newnuc -= 999;
      };

      // special case MCNP Am-242
      if (newnuc == 952421)
        newnuc -= 1;
    }
    return newnuc;
  }


  // Not a normal nuclide, might be a 
  // Natural elemental nuclide.  
  // ie for Urnaium = 920000
  if (mod_10000 == 0 && 0 == zz_name.count(div_1000) && 0 < zz_name.count(div_10000))
  {
    // zzaaam form natural nuclide
    newnuc = nuc;
  }
  else if (mod_1000 == 0 && mod_10000 != 0 && 0 < zz_name.count(div_1000))
  {
    // MCNP form natural nuclide
    newnuc = nuc * 10;
  }
  else
  {
    newnuc = -1;
    throw IndeterminateNuclideForm(nuc, "");
  };

  return newnuc;
};



int pyne::nucname::zzaaam(char * nuc)
{
  std::string newnuc (nuc);
  return zzaaam(newnuc);
};



int pyne::nucname::zzaaam(std::string nuc)
{
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");

  int newnuc, nuclen;
  std::string nucstr;

  // Get the string into a regular form
  nucstr = pyne::to_upper(nuc);
  nucstr = pyne::remove_substring(nucstr, "-");
  nuclen = nucstr.length();

  if (pyne::contains_substring(pyne::digits, nucstr.substr(0, 1)))
  {
    if (pyne::contains_substring(pyne::digits, nucstr.substr(nuclen-1, nuclen)))
    {
      // Nuclide must actually be an integer that 
      // just happens to be living in string form.
      newnuc = pyne::to_int(nucstr);
      newnuc = zzaaam(newnuc);
    }
    else
    {
      // probably in NIST-like form (242Am)
      // Here we know we have both digits and letters
      std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);
      newnuc = pyne::to_int(anum_str) * 10;

      // Add the Z-number
      std::string elem_name = pyne::remove_characters(nucstr, pyne::digits);
      if (0 < name_zz.count(elem_name))
        newnuc = (10000 * name_zz[elem_name]) + newnuc;
      else
        throw NotANuclide(nucstr, newnuc);
    };
  }
  else if (pyne::contains_substring(pyne::alphabet, nucstr.substr(0, 1)))
  {
    // Nuclide is probably in name form, or some variation therein
    std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);

    // natural element form, a la 'U' -> 920000
    if (anum_str.empty() && (0 < name_zz.count(nucstr)))
      return 10000 * name_zz[nucstr]; 

    int anum = pyne::to_int(anum_str);

    // bad form
    if (anum < 0)
      throw NotANuclide(nucstr, anum); 

    // Figure out if we are meta-stable or not
    std::string end_char = pyne::last_char(nucstr);
    if (end_char == "M")
      newnuc = (10 * anum) + 1;
    else if (pyne::contains_substring(pyne::digits, end_char))
      newnuc = (10 * anum);
    else
      throw NotANuclide(nucstr, newnuc);

    // Add the Z-number
    std::string elem_name = pyne::remove_characters(nucstr.substr(0, nuclen-1), pyne::digits);
    if (0 < name_zz.count(elem_name))
      newnuc = (10000 * name_zz[elem_name]) + newnuc;
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
  int nucint = zzaaam(nuc);
  std::string newnuc = "";

  int mod_10 = nucint%10;
  int mod_10000 = nucint % 10000;
  int div_10000 = nucint / 10000;
  int mod_10000_div_10 = mod_10000 / 10;

  // Make sure the LL value is correct
  if (0 == zz_name.count(div_10000))
    throw NotANuclide(nuc, nucint);

  // Add LL
  newnuc += zz_name[div_10000];

  // Add A-number
  if (0 < mod_10000)
    newnuc += pyne::to_str(mod_10000_div_10);

  // Add meta-stable flag
  if (0 < mod_10)
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
  int newnuc = zzaaam(nuc);
  return name(newnuc);
}





/**********************/
/*** mcnp functions ***/
/**********************/
int pyne::nucname::mcnp(int nuc)
{
  int newnuc = zzaaam(nuc);
  int mod_10 = newnuc%10;

  newnuc = newnuc/10;

  // special case Am242(m)
  if (newnuc == 95242 && mod_10 < 2)
  {
    mod_10 = (mod_10 + 1) % 2;
  }

  // Handle the crazy MCNP meta-stable format
  if (0 != mod_10) 
  {
    newnuc += 300;
    newnuc += (mod_10 * 100);
  }

  return newnuc;
};



int pyne::nucname::mcnp(char * nuc)
{
  std::string newnuc (nuc);
  return mcnp(newnuc);
};



int pyne::nucname::mcnp(std::string nuc)
{
  int newnuc = zzaaam(nuc);
  return mcnp(newnuc);
};





/*************************/
/*** serpent functions ***/
/*************************/
std::string pyne::nucname::serpent(int nuc)
{
  int nucint = zzaaam(nuc);
  std::string newnuc = "";

  int mod_10 = nucint%10;
  int mod_10000 = nucint % 10000;
  int div_10000 = nucint / 10000;
  int mod_10000_div_10 = mod_10000 / 10;

  // Make sure the LL value is correct
  if (0 == zz_name.count(div_10000))
    throw NotANuclide(nuc, nucint);

  // Add LL
  std::string LLupper = pyne::to_upper(zz_name[div_10000]);
  std::string LLlower = pyne::to_lower(zz_name[div_10000]);
  newnuc += LLupper[0];
  for (int l = 1; l < LLlower.size(); l++)
    newnuc += LLlower[l];  

  // Add required dash
  newnuc += "-";

  // Add A-number
  if (0 < mod_10000)
    newnuc += pyne::to_str(mod_10000_div_10);
  else if (0 == mod_10000)
    newnuc += "nat";

  // Add meta-stable flag
  if (0 < mod_10)
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
  int newnuc = zzaaam(nuc);
  return serpent(newnuc);
};




/**********************/
/*** nist functions ***/
/**********************/
std::string pyne::nucname::nist(int nuc)
{
  int nucint = zzaaam(nuc);
  std::string newnuc = "";

  int mod_10 = nucint%10;
  int mod_10000 = nucint % 10000;
  int div_10000 = nucint / 10000;
  int mod_10000_div_10 = mod_10000 / 10;

  // Make sure the LL value is correct
  if (0 == zz_name.count(div_10000))
    throw NotANuclide(nuc, nucint);

  // Add A-number
  if (0 < mod_10000)
    newnuc += pyne::to_str(mod_10000_div_10);

  // Add name
  std::string name_upper = pyne::to_upper(zz_name[div_10000]);
  std::string name_lower = pyne::to_lower(zz_name[div_10000]);
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
  int newnuc = zzaaam(nuc);
  return nist(newnuc);
};




/************************/
/*** cinder functions ***/
/************************/
int pyne::nucname::cinder(int nuc)
{
  // cinder nuclides of form aaazzzm

  int newnuc = zzaaam(nuc);
  int mod_10000 = newnuc % 10000;
  int div_10000 = newnuc / 10000;
  int mod_10000_div_10 = mod_10000 / 10;
  int mod_10 = newnuc%10;

  newnuc = (mod_10000_div_10*10000) + (div_10000*10) + mod_10;
  return newnuc;
};



int pyne::nucname::cinder(char * nuc)
{
  std::string newnuc (nuc);
  return cinder(newnuc);
};



int pyne::nucname::cinder(std::string nuc)
{
  int newnuc = zzaaam(nuc);
  return cinder(newnuc);
};




