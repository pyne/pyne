#ifndef PYNE_E6JL4RB2RFE2RCJEM4HLIAHO6M
#define PYNE_E6JL4RB2RFE2RCJEM4HLIAHO6M

#include <sstream>
#include <fstream>
#include <list>

#ifndef PYNE_IS_AMALGAMATED
  #include "pyne.h"
#endif

namespace pyne {
  namespace endf {

    /*****************************************************************************/
    /************************* Structs for MT data *******************************/
    /*****************************************************************************/

    /// basic structure for an mt section
    typedef struct mt_base {
      int nuc_id; ///< ENDF style nuc id ZZZAAA
      double awr; ///< mass relative to neutron [double]
      int mat; ///< ENDF material
      int mf; ///< ENDF file number
      int mt; ///< ENDF reaction designation
      /// force structure to be polymorphic
      virtual ~mt_base() {};
    } mt_base;

    /// mt section for descriptive informaiton
    typedef struct mt451 : mt_base {
      int lrp; ///< resonance parameters given
      int lfi; ///< does material fission
      int nlib; ///< library identifier
      int nmod; ///< material modification
      double elis; ///< excitation energy
      int sta; ///< target stability
      int lis; ///< state number of the nucleas
      int liso; ///< isomeric state number
      int nfor; ///< library format
      double awi; ///< projectile mass relative to neutron
      double emax; ///< upper energy limit
      int lrel; ///< library release number
      int nsub; ///< sub library number
      int nver; ///< library version number
      double temp; ///< temperature in kelvin
      int ldrv; ///< derived material flag
      int nwd; ///< number of descriptove records
      int nxc; ///< number of mt record in this file
      std::vector<std::vector<int> > mt_list; ///< list of [mf,mt,lines,mod]
    } mt451;

    /// read an mt 451 information from a filestream
    /// \param infile An ENDF filestream
    mt451 read_mt451(std::ifstream &infile);

    /// A typedef struct for fission profuct yield data MT 454, 459
    typedef struct mtfpy_mf8 : mt_base {
      int le; ///< number of energy  dependent yields given
      std::vector<int> i; ///< interpolation to be used between E[i-1] and E[i]
      std::vector<double> e; ///< list of energies
      std::vector<std::vector<std::vector<double> > > yields; ///< yield data [zafp,fps,yi,dyi]
    } mtfpy_mf8;

    /// read an mt fission product yield information from a filestream
    /// \param infile An ENDF filestream
    mtfpy_mf8 read_mtfpy_mf8(std::ifstream &infile);

    /// A typedef struct for neutrons per fission
    typedef struct mt452_mf1 : mt_base {
      int lnu; ///< type of data in section
      /// if LNU = 1 this will be an nth order polynomial describing nu
      std::vector<double> poly; ///< polynomial describing neutrons per fission(E)
      /// if LNU = 2 this will contain the tabulated data parameters
      std::vector<int> nbt; ///< list of interpolation segments
      std::vector<int> intn; ///< list of interpolations to be used
      std::vector<double> eint; ///< Energy of the incident neutron
      std::vector<double> nu_e; ///< Neutrons per fission at the given energy
    } mt452_mf1;

    /// read neutrons per fission information from a filestream
    /// \param infile An ENDF filestream
    mt452_mf1 read_mt452_mf1(std::ifstream &infile);

    typedef struct mt455_mf1 : mt_base {
      int ldg; ///< energy dependence of decay constants
      int lnu; ///< data representation type

      std::vector<double> lambdas; ///< decay constant of ith precursor

      std::vector<int> nbt; ///< list of interpolation segments
      std::vector<int> intn; ///< list of interpolations to be used
      std::vector<int> ne; ///< number of energies
      std::vector<double> eint; ///< ith energy
      std::vector<int> einti; ///< energy interpolation scheme
      std::vector<double> nu_d; ///< average neutrons per delayed fission event

      /// decay constant of ith precursor with neutron energy energy E
      std::vector<std::vector<double> > lambda_arr;
      /// fraction of ith precursor with neutron energy E
      std::vector<std::vector<double> > alpha_arr;

    } mt455_mf1;

    /// read precursor information from a filestream
    /// \param infile An ENDF filestream
    mt455_mf1 read_mt455_mf1(std::ifstream &infile);

    /// structure containing data about prompt neutrons per fission as a function
    /// of energy
    typedef struct mt456_mf1 : mt_base {
      int lnu;

      /// If lnu = 1 nu contains an nth order polynomial describing neutrons as a
      /// function of energy.
      std::vector<double> nu; ///< polynomial describing neutrons per fission(E)
      /// If lnu = 2 a list of tabulated values are give in tab1 style format
      std::vector<int> nbt; ///< list of interpolation segments
      std::vector<int> intn; ///< list of interpolations to be used
      std::vector<double> eint; ///< ith energy
      std::vector<double> nu_e; ///< average neutrons per prompt fission event

    } mt456_mf1;

    /// read neutrons per fission data from a filestream
    /// \param infile An ENDF filestream
    mt456_mf1 read_mt456_mf1(std::ifstream &infile);

    /// Data structure for fission energy release
    typedef struct mt458_mf1 : mt_base {
      std::vector<std::pair<double,double> > efr; ///< kinetic energy from fission products
      std::vector<std::pair<double,double> > pen; ///< kinetic energy of prompt fission neutrons
      std::vector<std::pair<double,double> > den; ///< kinetic energy of delayed fission neutrons
      std::vector<std::pair<double,double> > egp; ///< total energy released by prompt gamma rays
      std::vector<std::pair<double,double> > egd; ///< total energy released by delayed gamma rays
      std::vector<std::pair<double,double> > eb; ///< total energy released by delayed betas
      std::vector<std::pair<double,double> > enu; ///< energy carried away by neutrinos
      std::vector<std::pair<double,double> > er; ///< total energy release not lost to neutrinos
      std::vector<std::pair<double,double> > et; ///< total energy release per fission

    } mt458_mf1;

    /// reads ENDF fission energy release data from an ifstream and returns
    /// structure containing associated data
    /// \param infile An ENDF filestream
    mt458_mf1 read_mt458_mf1(std::ifstream &infile);

    /// data structure for delayed photon data
    typedef struct mt460_mf1 : mt_base {

      int lo; ///< representation type: 1 if discrete, 2 if continuous
      /// lo = 1
      int ng; ///< number
      std::vector<double> elist; ///< energy of the ith photon in eV
      std::vector<std::vector<int> > nbt; ///< list^2 interpolation segments
      std::vector<std::vector<int> > intn; ///< list^2 of interpolations
      std::vector<std::vector<double> > tint; ///< time of ith multiplicity
      std::vector<std::vector<double> > t; ///< Time dependence of ith multiplicity
      /// lo = 2
      std::vector<double> lambdas; ///< decay constant for the ith precursor

    } mt460_mf1;

    /// read delayed photon data from a filestream
    /// \param infile An ENDF filestream
    mt460_mf1 read_mt460_mf1(std::ifstream &infile);


    /// data structure for ENDF decay data
    typedef struct mt457_mf8 : mt_base {
      int lis; ///< state of the original nuclide
      int liso; ///< isomeric state of the original nuclide
      int nst; ///< nucleus stability flag: 0 = radioactive
      int nsp; ///< total number of radiation types for which information is given
      std::vector<std::pair<double,double> > erel; ///< energy released
      std::vector<double> styp; ///< Decay radiation
      std::vector<int> lcon; ///< Continuum spectrum flag
      std::vector<int> ner; ///< number of tabulated discrete energies
      std::vector<int> nd; ///< number of decay modes given 
      std::vector<std::pair<double,double> > fd; ///< discrete spectrum normalization factor
      std::vector<std::pair<double,double> > eav; ///< average decay energy of radiation
      std::vector<std::pair<double,double> > fc; ///< continuum spectrum normalization
      std::vector<std::pair<double,double> > er; ///< discrete energy of radiation
      std::vector<double> rtyp; ///< decay mode of the nuclide
      std::vector<double> type; ///< type of beta/ec transition
      std::vector<std::pair<double,double> > ri; ///< intensity of discrete radiation
      std::vector<std::pair<double,double> > ris; ///< internal pair formation coefficient
      std::vector<std::pair<double,double> > ricc; ///< total internal conversion coefficient
      std::vector<std::pair<double,double> > rick; ///< k-shell internal conversion coefficient
      std::vector<std::pair<double,double> > ricl; ///< l-shell internal conversion coefficient
    } mt457_mf8;

    /// read decay data from a filestream
    /// \param infile An ENDF filestream
    mt457_mf8 read_mt457_mf8(std::ifstream &infile);

    /// data structure for cross section data
    typedef struct mf3 : mt_base {
      std::vector<int> nbt; ///< list of interpolation segments
      std::vector<int> intn; ///< list of interpolations to be used
      std::vector<double> energy; ///< ith energy
      std::vector<double> sigma; ///< ith cross section value
    } mf3;
    /// read cross section data from a filestream
    /// \param infile An ENDF filestream
    mf3 read_mf3(std::ifstream &infile);
    
    typedef struct mf2 : mt_base {
      std::vector<double> zai;
      std::vector<double> abn;
      std::vector<int> lfw;
      std::vector<int> ner;
      std::vector<std::vector<double> > el;
      std::vector<std::vector<double> > eh;
      std::vector<std::vector<int> > lru;
      std::vector<std::vector<int> > lrf;
      std::vector<std::vector<int> > nro;
      std::vector<std::vector<int> > naps;
      std::vector<std::vector<double> > spi;
      std::vector<std::vector<double> > ap;
      std::vector<std::vector<int> > lad;
      std::vector<std::vector<int> > nls;
      std::vector<std::vector<int> > nlsc;
      std::vector<std::vector<int> > l;
      std::vector<std::vector<int> > njs;
      std::vector<std::vector<int> > ifg;
      std::vector<std::vector<int> > krm;
      std::vector<std::vector<int> > njs;
      std::vector<std::vector<int> > krl;
      std::vector<std::vector<std::vector<int> > > lbk;
      std::vector<std::vector<std::vector<int> > > lps;


    }
  }
}

#endif
