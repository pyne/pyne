// Converts between naming conventions for nuclides.
// zzaaam is for numerals only (923350).
// LLAAAM is for letters  as well (U-235).
// MCNP is for numerals without the meta-stable flag (92235), as used in MCNP.

#ifndef _ISONAME_
#define _ISONAME_
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <exception>
#include <stdlib.h>
#include <stdio.h>

#include "bright.h"

namespace isoname 
{
    /*** Constructs the LL to zz Dictionary ***/
    typedef std::string LL;
    typedef int zz;
    typedef std::map<LL, zz> LLzzType; 
    typedef LLzzType::iterator LLzzIter; 
    LLzzType get_LLzz();
    extern LLzzType LLzz; 

    /*** Constructs zz to LL dictionary **/
    typedef std::map<zz, LL> zzLLType; 
    typedef zzLLType::iterator zzLLIter; 
    zzLLType get_zzLL();
    extern zzLLType zzLL;

    /******************************************/
    /*** Define useful elemental group sets ***/
    /******************************************/
    typedef std::string LL_Name;
    typedef std::set<LL_Name> LL_Group;
    typedef LL_Group::iterator LL_GroupIter;
    typedef int zz_Name;
    typedef std::set<zz_Name> zz_Group;
    typedef zz_Group::iterator zz_GroupIter;

    zz_Group LL_2_zz_Group (LL_Group);

    extern std::string LAN_Array[15];
    extern LL_Group LAN;
    extern zz_Group lan;

    extern std::string ACT_Array[23];
    extern LL_Group ACT;
    extern zz_Group act;

    extern std::string TRU_Array[19];
    extern LL_Group TRU;
    extern zz_Group tru;

    extern std::string MA_Array[18];
    extern LL_Group MA;
    extern zz_Group ma;

    LL_Group get_FP();
    extern LL_Group FP;
    std::string * get_FP_Array(LL_Group);
    extern std::string * FP_Array;
    extern zz_Group fp;


    /******************/
    /*** Exceptions ***/
    /******************/

    class NotANuclide : public std::exception
    {
    public:
        NotANuclide () {};
        ~NotANuclide () throw () {};
        NotANuclide(std::string  wasptr, std::string nowptr)
        {
            nucwas = wasptr;
            nucnow = nowptr;
        };
        NotANuclide(std::string wasptr, int nowptr)
        {
            nucwas = wasptr;
            nucnow = bright::to_str(nowptr);
        };
        NotANuclide(int wasptr, std::string nowptr)
        {
            nucwas = bright::to_str(wasptr);
            nucnow = nowptr;
        };
        NotANuclide(int wasptr, int nowptr)
        {
            nucwas = bright::to_str(wasptr);
            nucnow = bright::to_str(nowptr);
        };
        virtual const char* what() const throw()
        {
            std::string NaNEstr ("Not a Nuclide! ");
            if (!nucwas.empty())
                NaNEstr += nucwas;
            if (!nucnow.empty())
            {
                NaNEstr += " --> "; 
                NaNEstr += nucnow;
            }
            return (const char *) NaNEstr.c_str();
        };
    private:
        std::string nucwas;
        std::string nucnow;
    };

    static class IndeterminateNuclideFormException : public std::exception
    {
        virtual const char* what() const throw()
        {
            return "Nuclide Form Could Not Be Determined!\nPlease Ensure Nuclide is of zzaaam, LLAAAM, or MCNP form.\n";
        }
    } IndeterminateNuclideForm;

    /********************/
    /*** Current Form ***/
    /********************/
    std::string CurrentForm(std::string);
    std::string CurrentForm(int);

    /****************************/
    /*** LLAAAM_2_* Functions ***/
    /****************************/
    int LLAAAM_2_zzaaam(std::string);
    int LLAAAM_2_MCNP(std::string);

    /***************************/
    /*** zzaaam_2_* Functions **/
    /***************************/
    std::string zzaaam_2_LLAAAM(int);
    int zzaaam_2_MCNP(int);

    /**************************/
    /*** MCNP_2_* Functions ***/
    /**************************/
    int MCNP_2_zzaaam(int);
    std::string MCNP_2_LLAAAM(int);

    /****************************/
    /*** mixed_2_*_ Functions ***/
    /****************************/
    int mixed_2_zzaaam(std::string);
    int mixed_2_zzaaam(int);
    std::string mixed_2_LLAAAM(std::string);
    std::string mixed_2_LLAAAM(int);
    int mixed_2_MCNP(std::string);
    int mixed_2_MCNP(int);

    /************************/
    /*** Helper Functions ***/
    /************************/
    double nuc_weight_zzaaam(int);
    double nuc_weight(int);
    double nuc_weight(std::string);
}

#endif
