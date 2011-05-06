// Converts between naming conventions for nuclides.
// zzaaam is for numerals only (923350).
// LLAAAM is for letters  as well (U-235).
// MCNP is for numerals without the meta-stable flag (92235), as used in MCNP.

#include "nucname.h"

/*** Constructs the LL to zz Dictionary ***/
nucname::LLzz_t nucname::get_LLzz()
{
  nucname::LLzz_t lzd;

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
nucname::LLzz_t nucname::LLzz = nucname::get_LLzz();


/*** Constructs zz to LL dictionary **/
nucname::zzLL_t nucname::get_zzLL()
{
  zzLL_t zld;
  for (LLzz_iter i = LLzz.begin(); i != LLzz.end(); i++)
  {
    zld[i->second] = i->first;
  }
  return zld;
};
nucname::zzLL_t nucname::zzLL = nucname::get_zzLL();



/******************************************/
/*** Define useful elemental group sets ***/
/******************************************/

nucname::zz_group nucname::LL_to_zz_group(nucname::LL_group eg)
{
  zz_group zg;
  for (LL_group_iter i = eg.begin(); i != eg.end(); i++)
  {
    zg.insert(LLzz[*i]);
  }
  return zg;
};

// Lanthanides
nucname::LL_t nucname::LAN_array[15] = {"LA", "CE", "PR", "ND", "PM", "SM", "EU", \
                               "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU"};
nucname::LL_group nucname::LAN (nucname::LAN_array, nucname::LAN_array+15);
nucname::zz_group nucname::lan = nucname::LL_to_zz_group(nucname::LAN);

// Actinides
nucname::LL_t nucname::ACT_array[23] = {"AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", \
                               "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", \
                               "DB", "SG", "BH", "HS", "MT", "DS", "RG"};
nucname::LL_group nucname::ACT (nucname::ACT_array, nucname::ACT_array+23);
nucname::zz_group nucname::act = nucname::LL_to_zz_group(nucname::ACT);

// Transuarnics
nucname::LL_t nucname::TRU_array[19] = {"NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", \
                               "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", \
                               "MT", "DS", "RG"};
nucname::LL_group nucname::TRU (nucname::TRU_array, nucname::TRU_array+19);
nucname::zz_group nucname::tru = nucname::LL_to_zz_group(nucname::TRU);

//Minor Actinides
nucname::LL_t nucname::MA_array[18] = {"NP", "AM", "CM", "BK", "CF", "ES", "FM", "MD", \
                              "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", \
                              "DS", "RG"};
nucname::LL_group nucname::MA (nucname::MA_array, nucname::MA_array+18);
nucname::zz_group nucname::ma = nucname::LL_to_zz_group(nucname::MA);

//Fission Products
nucname::LL_group nucname::get_FP()
{
  // Returns the Fission Product group
  LL_group FPs;
  for (LLzz_iter i = LLzz.begin(); i != LLzz.end(); i++)
  {
    if (ACT.count(i->first) == 0)
      FPs.insert(i->first);
  }
  return FPs;
}

nucname::LL_t * nucname::get_FP_array(nucname::LL_group FPs)
{
  // Returns the Fission Product group as an array
  LL_t * fp_array = new LL_t [FPs.size()];
  int n = 0;

  for (LL_group_iter i = FPs.begin(); i != FPs.end(); i++)
  {
    fp_array[n] = *i;
    n++;
  }
  return fp_array;
};

nucname::LL_group nucname::FP = nucname::get_FP();
nucname::LL_t * nucname::FP_array = nucname::get_FP_array(nucname::FP);
nucname::zz_group nucname::fp = nucname::LL_to_zz_group(nucname::FP);



/********************/
/*** Current Form ***/
/********************/
std::string nucname::current_form(std::string nuc)
{
  // returns current form of a nuclide.
  using namespace pyne;

  nuc = to_upper(nuc);
  nuc = remove_substring(nuc, "-");

  if ( contains_substring(alphabet, nuc.substr(0,1)) )
    return "LLAAAM";
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

std::string nucname::current_form(int nuc)
{
   return nucname::current_form(pyne::to_str(nuc));
};




/************************/
/*** zzaaam functions ***/
/************************/
int nucname::zzaaam(int nuc)
{
  int newnuc;

  int mod_10000 = nuc % 10000;
  int div_10000 = nuc / 10000;
  int mod_10000_div_10 = mod_10000 / 10;

  // Nuclide must already be in zzaaam form
  if (div_10000 <= mod_10000_div_10 && mod_10000_div_10 <= div_10000 * 6)
  {
    // Normal nuclide
    newnuc = nuc;
    return newnuc;
  }
  else if (mod_10000 == 0 && 0 < zzLL.count(div_10000))
  {
    // Natural elemental nuclide:  ie for Urnaium = 920000
    newnuc = nuc;
    return newnuc;
  }

  // Nuclide is not in zzaaam form, 
  // Try MCNP form, ie zzaaa
  int mod_1000 = nuc % 1000; 
  int div_1000 = nuc / 1000;
  int mod_1000_div_10 = mod_10000 / 10;

  if (div_1000 <= mod_1000_div_10 && mod_1000_div_10 <= div_1000 * 6)
  {
    if (mod_1000 - 400 < 0)
    {
      // Nuclide in normal MCNP form
      newnuc = nuc * 10;
    }
    else
    {
      // Nuclide in MCNP metastable form
      newnuc = ((nuc - 400) * 10) + 1;
      while (2.5 < (float ((newnuc/10)%1000) / float (newnuc/10000)))
      {
        newnuc -= 999;
      };
    }
    return newnuc;
  }


  // Not a normal nuclide, might be a 
  // Natural elemental nuclide.  
  // ie for Urnaium = 920000
  if (mod_10000 == 0 && 0 == zzLL.count(div_1000) && 0 < zzLL.count(div_10000))
  {
    // zzaaam form natural nuclide
    newnuc = nuc;
  }
  else if (mod_1000 == 0 && mod_10000 != 0 && 0 < zzLL.count(div_1000))
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



int nucname::zzaaam(char * nuc)
{
  std::string newnuc (nuc);
  return zzaaam(newnuc);
};



int nucname::zzaaam(std::string nuc)
{
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");

  int newnuc;
  std::string nucstr;

  // Get the string into a regular form
  nucstr = pyne::to_upper(nuc);
  nucstr = pyne::remove_substring(nucstr, "-");

  if (pyne::contains_substring(pyne::digits, &nucstr[0]))
  {
    // Nuclide must actually be an integer that 
    // just happens to be living in string form.
    newnuc = pyne::to_int(nucstr);
    newnuc = zzaaam(newnuc);
  }
  else if (pyne::contains_substring(pyne::alphabet, &nucstr[0]))
  {
    // Nuclide is probably in LLAAAM form, or some variation therein
    int anum = pyne::to_int(pyne::remove_characters(nucstr, pyne::alphabet));
    if (anum < 1)
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
    std::string lnum = pyne::remove_characters(nucstr.substr(0, nucstr.length()-1), pyne::digits);
    if (0 < LLzz.count(lnum))
      newnuc = (10000 * LLzz[lnum]) + newnuc;
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




/************************/
/*** LLAAAM functions ***/
/************************/
std::string nucname::LLAAAM(int nuc)
{
  int nucint = zzaaam(nuc);
  std::string newnuc = "";

  int mod_10 = nucint%10;
  int mod_10000 = nuc % 10000;
  int div_10000 = nuc / 10000;
  int mod_10000_div_10 = mod_10000 / 10;

  // Make sure the LL value is correct
  if (0 == zzLL.count(div_10000))
    throw NotANuclide(nuc, nucint);

  // Add LL
  newnuc += zzLL[div_10000];

  // Add A-number
  if (0 < mod_10000)
    newnuc += pyne::to_str(mod_10000_div_10);

  // Add meta-stable flag
  if (0 < mod_10)
    newnuc += "M";

  return newnuc;
};



std::string nucname::LLAAAM(char * nuc)
{
  std::string newnuc (nuc);
  return LLAAAM(newnuc);
}


std::string nucname::LLAAAM(std::string nuc)
{
  int newnuc = zzaaam(nuc);
  return LLAAAM(newnuc);
}





/**********************/
/*** mcnp functions ***/
/**********************/
int nucname::mcnp(int nuc)
{
  int newnuc = zzaaam(nuc);
  int mod_10 = newnuc%10;

  newnuc = newnuc/10;

  // Handle the crazy MCNP meta-stable format
  if (0 != mod_10)
  {
    newnuc += 300;
    newnuc += (mod_10 * 100);
  }

  return newnuc;
};



int nucname::mcnp(char * nuc)
{
  std::string newnuc (nuc);
  return mcnp(newnuc);
};



int nucname::mcnp(std::string nuc)
{
  int newnuc = zzaaam(nuc);
  return mcnp(newnuc);
};





/*************************/
/*** serpent functions ***/
/*************************/
std::string nucname::serpent(int nuc)
{
  int nucint = zzaaam(nuc);
  std::string newnuc = "";

  int mod_10 = nucint%10;
  int mod_10000 = nuc % 10000;
  int div_10000 = nuc / 10000;
  int mod_10000_div_10 = mod_10000 / 10;

  // Make sure the LL value is correct
  if (0 == zzLL.count(div_10000))
    throw NotANuclide(nuc, nucint);

  // Add LL
  std::string LLupper = pyne::to_upper(zzLL[div_10000]);
  std::string LLlower = pyne::to_lower(zzLL[div_10000]);
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


std::string nucname::serpent(char * nuc)
{
  std::string newnuc (nuc);
  return serpent(newnuc);
};


std::string nucname::serpent(std::string nuc)
{
  int newnuc = zzaaam(nuc);
  return serpent(newnuc);
};




/************************/
/*** Helper Functions ***/
/************************/
double nucname::nuc_weight(int nuc)
{
    int nuc_zz = zzaaam(nuc);
    return (double) ((nuc_zz/10)%1000);
};

double nucname::nuc_weight(char * nuc)
{
    int nuc_zz = zzaaam(nuc);
    return nuc_weight(nuc_zz);
};

double nucname::nuc_weight(std::string nuc)
{
    int nuc_zz = zzaaam(nuc);
    return nuc_weight(nuc_zz);
};
