// Converts between naming conventions for nuclides.
// zzaaam is for numerals only (923350).
// name is for letters  as well (U-235).
// MCNP is for numerals without the meta-stable flag (92235), as used in MCNP.

#include "nucname.h"

/*** Constructs the LL to zz Dictionary ***/
nucname::name_zz_t nucname::get_name_zz()
{
  nucname::name_zz_t lzd;

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
nucname::name_zz_t nucname::name_zz = nucname::get_name_zz();


/*** Constructs zz to LL dictionary **/
nucname::zzname_t nucname::get_zz_name()
{
  zzname_t zld;
  for (name_zz_iter i = name_zz.begin(); i != name_zz.end(); i++)
  {
    zld[i->second] = i->first;
  }
  return zld;
};
nucname::zzname_t nucname::zz_name = nucname::get_zz_name();



/******************************************/
/*** Define useful elemental group sets ***/
/******************************************/

nucname::zz_group nucname::name_to_zz_group(nucname::name_group eg)
{
  zz_group zg;
  for (name_group_iter i = eg.begin(); i != eg.end(); i++)
  {
    zg.insert(name_zz[*i]);
  }
  return zg;
};

// Lanthanides
nucname::name_t nucname::LAN_array[15] = {"LA", "CE", "PR", "ND", "PM", "SM", "EU", \
                                        "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU"};
nucname::name_group nucname::LAN (nucname::LAN_array, nucname::LAN_array+15);
nucname::zz_group nucname::lan = nucname::name_to_zz_group(nucname::LAN);

// Actinides
nucname::name_t nucname::ACT_array[15] = {"AC", "TH", "PA", "U",  "NP", "PU", "AM", "CM", \
                                        "BK", "CF", "ES", "FM", "MD", "NO", "LR"};
nucname::name_group nucname::ACT (nucname::ACT_array, nucname::ACT_array+15);
nucname::zz_group nucname::act = nucname::name_to_zz_group(nucname::ACT);

// Transuarnics
nucname::name_t nucname::TRU_array[19] = {"NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", \
                                        "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", \
                                        "MT", "DS", "RG"};
nucname::name_group nucname::TRU (nucname::TRU_array, nucname::TRU_array+19);
nucname::zz_group nucname::tru = nucname::name_to_zz_group(nucname::TRU);

//Minor Actinides
nucname::name_t nucname::MA_array[10] = {"NP", "AM", "CM", "BK", "CF", "ES", "FM", "MD", \
                                       "NO", "LR"};
nucname::name_group nucname::MA (nucname::MA_array, nucname::MA_array+10);
nucname::zz_group nucname::ma = nucname::name_to_zz_group(nucname::MA);

//Fission Products
nucname::name_t nucname::FP_array[88] = {"AG", "AL", "AR", "AS", "AT", "AU", "B",  "BA", \
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
nucname::name_group nucname::FP (nucname::FP_array, nucname::FP_array+88);
nucname::zz_group nucname::fp = nucname::name_to_zz_group(nucname::FP);



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
    }
    else
    {
      // Nuclide in MCNP metastable form
      newnuc = ((nuc - 400) * 10) + 1;
      while (3.0 < (float ((newnuc/10)%1000) / float (newnuc/10000)))
      {
        newnuc -= 999;
      };
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



int nucname::zzaaam(char * nuc)
{
  std::string newnuc (nuc);
  return zzaaam(newnuc);
};



int nucname::zzaaam(std::string nuc)
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
std::string nucname::name(int nuc)
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



std::string nucname::name(char * nuc)
{
  std::string newnuc (nuc);
  return name(newnuc);
}


std::string nucname::name(std::string nuc)
{
  int newnuc = zzaaam(nuc);
  return name(newnuc);
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




/**********************/
/*** nist functions ***/
/**********************/
std::string nucname::nist(int nuc)
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


std::string nucname::nist(char * nuc)
{
  std::string newnuc (nuc);
  return nist(newnuc);
};


std::string nucname::nist(std::string nuc)
{
  int newnuc = zzaaam(nuc);
  return nist(newnuc);
};






/****************************/
/*** Nuc_weight Functions ***/
/****************************/
std::map<int, double> nucname::nuc_weight_map = std::map<int, double>();


void nucname::_load_nuc_weight_map()
{
  // Loads the importnat parts of atomic_wight table into nuc_weight_map

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool isH5 = H5::H5File::isHdf5(pyne::NUC_DATA_PATH);
  if (!isH5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  H5::CompType atomic_weight_desc(sizeof(atomic_weight_struct));
  atomic_weight_desc.insertMember("nuc_name", HOFFSET(atomic_weight_struct, nuc_name), H5::StrType(0, 6));
  atomic_weight_desc.insertMember("nuc_zz",   HOFFSET(atomic_weight_struct, nuc_zz),   H5::PredType::NATIVE_INT);
  atomic_weight_desc.insertMember("mass",     HOFFSET(atomic_weight_struct, mass),     H5::PredType::NATIVE_DOUBLE);
  atomic_weight_desc.insertMember("error",    HOFFSET(atomic_weight_struct, error),    H5::PredType::NATIVE_DOUBLE);
  atomic_weight_desc.insertMember("abund",    HOFFSET(atomic_weight_struct, abund),    H5::PredType::NATIVE_DOUBLE);

  // Open the HDF5 file
  H5::H5File nuc_data_h5 (pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY);

  // Open the data set
  H5::DataSet atomic_weight_set = nuc_data_h5.openDataSet("/atomic_weight");
  H5::DataSpace atomic_weight_space = atomic_weight_set.getSpace();
  int atomic_weight_length = atomic_weight_space.getSimpleExtentNpoints();

  // Read in the data
  atomic_weight_struct * atomic_weight_array = new atomic_weight_struct[atomic_weight_length];
  atomic_weight_set.read(atomic_weight_array, atomic_weight_desc);

  // close the nuc_data library, before doing anythng stupid
  nuc_data_h5.close();

  // Ok now that we have the array of stucts, put it in a the map
  for(int n = 0; n < atomic_weight_length; n++)
    nuc_weight_map[atomic_weight_array[n].nuc_zz] = atomic_weight_array[n].mass;
};


double nucname::nuc_weight(int nuc)
{
  // Find the nuclide;s weight in AMU
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = nuc_weight_map.find(nuc);
  nuc_end = nuc_weight_map.end();

  // First check if we already have the nuc weight in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (nuc_weight_map.empty())
  {
    // Don't fail if we can't load the library
    try
    {
      _load_nuc_weight_map();
      return nuc_weight(nuc);
    }
    catch(...){};
  };

  double aw;
  int nuc_zz = zzaaam(nuc);

  // If in an excited state, return the ground
  // state weight...not strictly true, but good guess.
  if (0 < nuc_zz%10)
  {
    aw = nuc_weight((nuc_zz/10)*10);
    nuc_weight_map[nuc] = aw;
    return aw;
  };

  // Finally, if none of these work, 
  // take a best guess based on the 
  // aaa number.
  aw = (double) ((nuc_zz/10)%1000);
  nuc_weight_map[nuc] = aw;
  return aw;
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
