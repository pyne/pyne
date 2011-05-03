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
LL_t nucname::LAN_array[15] = {"LA", "CE", "PR", "ND", "PM", "SM", "EU", \
                               "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU"};
nucname::LL_group nucname::LAN (nucname::LAN_array, nucname::LAN_array+15);
nucname::zz_group nucname::lan = nucname::LL_to_zz_group(nucname::LAN);

// Actinides
LL_t nucname::ACT_array[23] = {"AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", \
                               "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", \
                               "DB", "SG", "BH", "HS", "MT", "DS", "RG"};
nucname::LL_group nucname::ACT (nucname::ACT_array, nucname::ACT_array+23);
nucname::zz_group nucname::act = nucname::LL_to_zz_group(nucname::ACT);

// Transuarnics
LL_t nucname::TRU_array[19] = {"NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", \
                               "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", \
                               "MT", "DS", "RG"};
nucname::LL_group nucname::TRU (nucname::TRU_array, nucname::TRU_array+19);
nucname::zz_group nucname::tru = nucname::LL_to_zz_group(nucname::TRU);

//Minor Actinides
LL_t nucname::MA_array[18] = {"NP", "AM", "CM", "BK", "CF", "ES", "FM", "MD", \
                              "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", \
                              "DS", "RG"};
nucname::LL_group nucname::MA (nucname::MA_array, nucname::MA_array+18);
nucname::zz_group nucname::ma = nucname::LL_2_zz_group(nucname::MA);

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

LL_t * nucname::get_FP_array(nucname::LL_group FPs)
{
  // Returns the Fission Product group as an array
  LL_t * fp_array = new LL_t [FPs.size()];
  int n = 0;

  for (LL_GroupIter i = FPs.begin(); i != FPs.end(); i++)
  {
    fp_array[n] = *i;
    n++;
  }
  return fp_array;
};

nucname::LL_group nucname::FP = nucname::get_FP();
LL_t * nucname::FP_array = nucname::get_FP_array(nucname::FP);
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
      throw IndeterminateNuclideForm;
  };
};

std::string nucname::current_form(int nuc)
{
   return nucname::current_form(pyne::to_str(nuc));
};




/************************/
/*** zzaaam functions ***/
/************************/
int zzaaam(int nuc)
{
    int newnuc;
    //Converts nuclide from MCNP form to zzaaam form.
    if ( (nuc%1000)-400 < 0 )
        return nuc * 10;
    else
    {
        //Please make more general so that more that the first metastable state is returned...
        return (nuc - 400)*10 + 1;
    }
};




/****************************/
/*** LLAAAM_2_* Functions ***/
/****************************/
int nucname::LLAAAM_2_zzaaam(std::string nucstr)
{
    //Converts nuclide from LLAAAM form to zzaaam form.
    using namespace pyne;

    int newnuc;
    nucstr = ToUpper(nucstr);
    nucstr = Strip(nucstr, "-");

    int anum =  to_int( MultiStrip(nucstr, alphabet) );
    if (anum < 1)
        throw NotANuclide (nucstr, anum);

    if ( LastChar(nucstr) == "M" )
        newnuc = (10 * anum) + 1;
    else if ( SubInString(digits, LastChar(nucstr)) )
        newnuc = (10 * anum);
    else
        throw NotANuclide (nucstr, newnuc);

    std::string lnum = MultiStrip(nucstr.substr(0, nucstr.length()-1), digits);

    if ( 0 < LLzz.count(lnum) )
        newnuc = (10000 * LLzz[lnum]) + newnuc;
    else
        throw NotANuclide (nucstr, newnuc);
    return newnuc;
};


int nucname::LLAAAM_2_MCNP(std::string nucstr)
{
    //Converts nuclide from LLAAAM form to MCNP form.
    return nucname::zzaaam_2_MCNP( nucname::LLAAAM_2_zzaaam(nucstr) );
};

/***************************/
/*** zzaaam_2_* Functions **/
/***************************/
std::string nucname::zzaaam_2_LLAAAM(int nuc)
{
    //Converts nuclide from azzzm form to LLAAAM form.
//	using namespace nucname;
    using namespace pyne;

    std::string newnuc = "";
    std::string nucstr = to_str( nuc );
    nucstr = Strip(nucstr, "-");
    
    try
    {
        if ( LastChar(nucstr) == "1")
            newnuc = to_str(to_int(SubFromEnd(nucstr,-4,3))) + "M";
        else if ( LastChar(nucstr) == "0")
            newnuc = to_str(to_int(SubFromEnd(nucstr,-4,3)));
        else
            throw NotANuclide (nucstr, newnuc);
    }
    catch (NotANuclide& e)
    {
        throw NotANuclide (nucstr, newnuc);
    }

    int znum = to_int( nucstr.substr(0, nucstr.length()-4) );
    if ( 0 < zzLL.count(znum) )
        newnuc = zzLL[znum] + newnuc;
    else
    {
        newnuc = "LL" + newnuc;
        throw NotANuclide (nucstr, newnuc);
    }
    return newnuc;
};

int nucname::zzaaam_2_MCNP(int nuc)
{
    //Converts nuclide from zzaaam form to MCNP form.
    if ( nuc%10 == 0 )
        return nuc / 10;
    else
    {
        int znum = nuc / 10000;
        int anum = (nuc/10) - (znum*1000) + 300;
        anum = anum + ((nuc%10)*100);
        return (znum*1000) + anum;
    }
};

/**************************/
/*** MCNP_2_* Functions ***/
/**************************/
int nucname::MCNP_2_zzaaam(int nuc)
{
    //Converts nuclide from MCNP form to zzaaam form.
    if ( (nuc%1000)-400 < 0 )
        return nuc * 10;
    else
    {
        //Please make more general so that more that the first metastable state is returned...
        return (nuc - 400)*10 + 1;
    }
};

std::string nucname::MCNP_2_LLAAAM(int nuc)
{
    //Converts nuclide from MCNP form to LLAAAM form.
    return nucname::zzaaam_2_LLAAAM( nucname::MCNP_2_zzaaam(nuc) );
};

/****************************/
/*** mixed_2_*_ Functions ***/
/****************************/
int nucname::mixed_2_zzaaam(std::string nuc)
{
    //Converts nuclide from mixed form to zzaaam form.
//	using namespace nucname;
    using namespace pyne;

    nuc = ToUpper(nuc);
    nuc = Strip(nuc, "-");
    std::string currentform = CurrentForm(nuc);
    if (currentform == "zzaaam")
        return to_int(nuc);
    else if (currentform == "LLAAAM")
        return LLAAAM_2_zzaaam(nuc);
    else if (currentform == "MCNP")
        return MCNP_2_zzaaam( to_int(nuc) );
    else
        throw IndeterminateNuclideForm;
};

int nucname::mixed_2_zzaaam(int nuc)
{
    return nucname::mixed_2_zzaaam( pyne::to_str(nuc) );
};

std::string nucname::mixed_2_LLAAAM(std::string nuc)
{
    //Converts nuclide from mixed form to LLAAAM form.
//	using namespace nucname;
    using namespace pyne;

    nuc = ToUpper(nuc);
    nuc = Strip(nuc, "-");
    std::string currentform = CurrentForm(nuc);
    if (currentform == "zzaaam")
        return zzaaam_2_LLAAAM( to_int(nuc) );
    else if (currentform == "LLAAAM")
        return nuc;
    else if (currentform == "MCNP")
        return MCNP_2_LLAAAM( to_int(nuc) );
    else
        throw IndeterminateNuclideForm;
};

std::string nucname::mixed_2_LLAAAM(int nuc)
{
    return nucname::mixed_2_LLAAAM( pyne::to_str(nuc) );
};

int nucname::mixed_2_MCNP(std::string nuc)
{
    //Converts nuclide from mixed form to MCNP form.
//	using namespace nucname;
    using namespace pyne;

    nuc = ToUpper(nuc);
    nuc = Strip(nuc, "-");
    std::string currentform = CurrentForm(nuc);
    if (currentform == "zzaaam")
        return zzaaam_2_MCNP( to_int(nuc) );
    else if (currentform == "LLAAAM")
        return LLAAAM_2_MCNP(nuc);
    else if (currentform == "MCNP")
        return to_int(nuc);
    else
        throw IndeterminateNuclideForm;
};

int nucname::mixed_2_MCNP(int nuc)
{
    return nucname::mixed_2_MCNP( pyne::to_str(nuc) );
};

/************************/
/*** Helper Functions ***/
/************************/
double nucname::nuc_weight_zzaaam(int nuc)
{
    return (double) ((nuc/10)%1000);
};

double nucname::nuc_weight(int nuc)
{
    int nuc_zz = nucname::mixed_2_zzaaam(nuc);
    return nucname::nuc_weight_zzaaam(nuc_zz);
};

double nucname::nuc_weight(std::string nuc)
{
    int nuc_zz = nucname::mixed_2_zzaaam(nuc);
    return nucname::nuc_weight_zzaaam(nuc_zz);
};

