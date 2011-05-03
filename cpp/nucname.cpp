// Converts between naming conventions for nuclides.
// zzaaam is for numerals only (923350).
// LLAAAM is for letters  as well (U-235).
// MCNP is for numerals without the meta-stable flag (92235), as used in MCNP.

#include "isoname.h"

/*** Constructs the LL to zz Dictionary ***/
isoname::LLzzType isoname::get_LLzz()
{
    isoname::LLzzType lzd;
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
isoname::LLzzType isoname::LLzz = isoname::get_LLzz();


/*** Constructs zz to LL dictionary **/
isoname::zzLLType isoname::get_zzLL()
{
    using namespace isoname;
    zzLLType zld;
    for (LLzzIter i = LLzz.begin(); i != LLzz.end(); i++)
    {
        zld[i->second] = i->first;
    }
    return zld;
};
isoname::zzLLType isoname::zzLL = isoname::get_zzLL();


/******************************************/
/*** Define useful elemental group sets ***/
/******************************************/

isoname::zz_Group isoname::LL_2_zz_Group (isoname::LL_Group eg)
{
    using namespace isoname;
    zz_Group zg;
    for (LL_GroupIter i = eg.begin(); i != eg.end(); i++)
    {
        zg.insert( LLzz[*i] );
    }
    return zg;
}

//Lanthanides
std::string isoname::LAN_Array[15] = {"LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU"};
isoname::LL_Group isoname::LAN (isoname::LAN_Array, isoname::LAN_Array+15);
isoname::zz_Group isoname::lan = isoname::LL_2_zz_Group( isoname::LAN );

// Actinides
std::string isoname::ACT_Array[23] = {"AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG"};
isoname::LL_Group isoname::ACT (isoname::ACT_Array, isoname::ACT_Array+23);
isoname::zz_Group isoname::act = isoname::LL_2_zz_Group( isoname::ACT );

//Transuarnics
std::string isoname::TRU_Array[19] = {"NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG"};
isoname::LL_Group isoname::TRU (isoname::TRU_Array, isoname::TRU_Array+19);
isoname::zz_Group isoname::tru = isoname::LL_2_zz_Group( isoname::TRU );

//Minor Actinides
std::string isoname::MA_Array[18] = {"NP", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG"};
isoname::LL_Group isoname::MA (isoname::MA_Array, isoname::MA_Array+18);
isoname::zz_Group isoname::ma = isoname::LL_2_zz_Group( isoname::MA );

//Fission Products
isoname::LL_Group isoname::get_FP()
{
    using namespace isoname;
    LL_Group FPs;
    for (LLzzIter i = LLzz.begin(); i != LLzz.end(); i++)
    {
        if (ACT.count(i->first) == 0)
            FPs.insert(i->first);
    }
    return FPs;
}

std::string * isoname::get_FP_Array(isoname::LL_Group FPs)
{
    using namespace isoname;
    std::string *fp_array = new std::string[FPs.size()];
    int n = 0;
    for (LL_GroupIter i = FPs.begin(); i != FPs.end(); i++)
    {
        fp_array[n] = *i;
        n++;
    }
    return fp_array;
}
isoname::LL_Group isoname::FP = isoname::get_FP();
std::string * isoname::FP_Array = isoname::get_FP_Array( isoname::FP );
isoname::zz_Group isoname::fp = isoname::LL_2_zz_Group( isoname::FP );

/********************/
/*** Current Form ***/
/********************/
std::string isoname::CurrentForm(std::string nuc)
{
    //returns current form of a nuclide.
//	using namespace isoname;
    using namespace bright;

    nuc = ToUpper(nuc);
    nuc = Strip(nuc, "-");

    if ( SubInString(alphabet, nuc.substr(0,1)) )
        return "LLAAAM";
    else
    {
        if (nuc.length() == 7)
            return "zzaaam";
        else if ( SubInString("23456789", LastChar(nuc)) )
            return "MCNP";

        if ( ChainGreaterCompare( to_int(nuc.substr(0,nuc.length()-4)), to_int(SubFromEnd(nuc,-4,3)), to_int(nuc.substr(0,nuc.length()-4)) * 5) )
            return "zzaaam";
        else if ( ChainGreaterCompare( to_int(nuc.substr(0,nuc.length()-3)), to_int(SubFromEnd(nuc,-3,3)), to_int(nuc.substr(0,nuc.length()-3)) * 5) )
            return "MCNP";
        else
            throw IndeterminateNuclideForm;
    }
};

std::string isoname::CurrentForm(int nuc)
{
    return isoname::CurrentForm( bright::to_str(nuc) );
};

/****************************/
/*** LLAAAM_2_* Functions ***/
/****************************/
int isoname::LLAAAM_2_zzaaam(std::string nucstr)
{
    //Converts nuclide from LLAAAM form to zzaaam form.
//	using namespace isoname;
    using namespace bright;

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


int isoname::LLAAAM_2_MCNP(std::string nucstr)
{
    //Converts nuclide from LLAAAM form to MCNP form.
    return isoname::zzaaam_2_MCNP( isoname::LLAAAM_2_zzaaam(nucstr) );
};

/***************************/
/*** zzaaam_2_* Functions **/
/***************************/
std::string isoname::zzaaam_2_LLAAAM(int nuc)
{
    //Converts nuclide from azzzm form to LLAAAM form.
//	using namespace isoname;
    using namespace bright;

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

int isoname::zzaaam_2_MCNP(int nuc)
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
int isoname::MCNP_2_zzaaam(int nuc)
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

std::string isoname::MCNP_2_LLAAAM(int nuc)
{
    //Converts nuclide from MCNP form to LLAAAM form.
    return isoname::zzaaam_2_LLAAAM( isoname::MCNP_2_zzaaam(nuc) );
};

/****************************/
/*** mixed_2_*_ Functions ***/
/****************************/
int isoname::mixed_2_zzaaam(std::string nuc)
{
    //Converts nuclide from mixed form to zzaaam form.
//	using namespace isoname;
    using namespace bright;

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

int isoname::mixed_2_zzaaam(int nuc)
{
    return isoname::mixed_2_zzaaam( bright::to_str(nuc) );
};

std::string isoname::mixed_2_LLAAAM(std::string nuc)
{
    //Converts nuclide from mixed form to LLAAAM form.
//	using namespace isoname;
    using namespace bright;

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

std::string isoname::mixed_2_LLAAAM(int nuc)
{
    return isoname::mixed_2_LLAAAM( bright::to_str(nuc) );
};

int isoname::mixed_2_MCNP(std::string nuc)
{
    //Converts nuclide from mixed form to MCNP form.
//	using namespace isoname;
    using namespace bright;

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

int isoname::mixed_2_MCNP(int nuc)
{
    return isoname::mixed_2_MCNP( bright::to_str(nuc) );
};

/************************/
/*** Helper Functions ***/
/************************/
double isoname::nuc_weight_zzaaam(int nuc)
{
    return (double) ((nuc/10)%1000);
};

double isoname::nuc_weight(int nuc)
{
    int nuc_zz = isoname::mixed_2_zzaaam(nuc);
    return isoname::nuc_weight_zzaaam(nuc_zz);
};

double isoname::nuc_weight(std::string nuc)
{
    int nuc_zz = isoname::mixed_2_zzaaam(nuc);
    return isoname::nuc_weight_zzaaam(nuc_zz);
};

