#ifndef PYNE_E6JL4RB2RFE2RCJEM4HLIAHO6M
#define PYNE_E6JL4RB2RFE2RCJEM4HLIAHO6M

#include <sstream>
#include <fstream>
#include <list>

#ifndef PYNE_IS_AMALGAMATED
  #include "pyne.h"
#endif

namespace pyne
{


  /*****************************************************************************/
  /************************* Structs for MT data *******************************/
  /*****************************************************************************/

  /// basic structure for an mt section
  class mt_base {
    public:
      int nuc_id; ///< ENDF style nuc id ZZZAAA
      double awr; ///< mass relative to neutron [double]
      int mat; ///< ENDF material
      int mf; ///< ENDF file number
      int mt; ///< ENDF reaction designation
      /// force struct to be polymorphic
      virtual ~mt_base() {};
  };

  /// mt section for descriptive informaiton
  class mt451 : public mt_base {
  public:
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
  };

  /// read an mt 451 information from a filestream
  /// \param infile An ENDF filestream
  mt451 read_mt451(std::ifstream &infile);

  /// A struct for fission profuct yield data MT 454, 459
  class mtfpy_mf8 : public mt_base {
  public:
    int le; ///< number of energy  dependent yields given
    std::vector<int> i; ///< interpolation to be used between E[i-1] and E[i]
    std::vector<double> e; ///< list of energies
    std::vector<std::vector<std::vector<double> > > yields; ///< yield data [zafp,fps,yi,dyi]
  };

  /// read an mt fission product yield information from a filestream
  /// \param infile An ENDF filestream
  mtfpy_mf8 read_mtfpy_mf8(std::ifstream &infile);

  /// A struct for neutrons per fission
  class mt452_mf1 : public mt_base {
  public:
    int lnu; ///< type of data in section
    /// if LNU = 1 this will be an nth order polynomial describing nu
    std::vector<double> poly;///< polynomial describing neutrons per fission(E)
    /// if LNU = 2 this will contain the tabulated data parameters
    std::vector<int> nbt; ///< list of interpolation segments
    std::vector<int> intn; ///< list of interpolations to be used
    std::vector<double> eint; ///< Energy of the incident neutron
    std::vector<double> nu_e; ///< Neutrons per fission at the given energy
  };

  /// read neutrons per fission information from a filestream
  /// \param infile An ENDF filestream
  mt452_mf1 read_mt452_mf1(std::ifstream &infile);

  class mt455_mf1 : public mt_base {
  public:

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

  };

  /// read precursor information from a filestream
  /// \param infile An ENDF filestream
  mt455_mf1 read_mt455_mf1(std::ifstream &infile);

  /// structure containing data about prompt neutrons per fission as a function
  /// of energy
  class mt456_mf1 : public mt_base {
  public:
    int lnu;

    /// If lnu = 1 nu contains an nth order polynomial describing neutrons as a
    /// function of energy.
    std::vector<double> nu;///< polynomial describing neutrons per fission(E)
    /// If lnu = 2 a list of tabulated values are give in tab1 style format
    std::vector<int> nbt; ///< list of interpolation segments
    std::vector<int> intn; ///< list of interpolations to be used
    std::vector<double> eint; ///< ith energy
    std::vector<double> nu_e; ///< average neutrons per prompt fission event

  };

  /// read neutrons per fission data from a filestream
  /// \param infile An ENDF filestream
  mt456_mf1 read_mt456_mf1(std::ifstream &infile);

  /// Data structure for fission energy release
  class mt458_mf1 : public mt_base {
  public:

    std::vector<double> efr; ///< kinetic energy from fission products
    std::vector<double> defr; ///< error in efr
    std::vector<double> enp; ///< kinetic energy of prompt fission neutrons
    std::vector<double> denp; ///< error in denp
    std::vector<double> end; ///< kinetic energy of delayed fission neutrons
    std::vector<double> dend; ///< error in dend
    std::vector<double> egp; ///< total energy released by prompt gamma rays
    std::vector<double> degp; ///< error in degp
    std::vector<double> egd; ///< total energy released by delayed gamma rays
    std::vector<double> degd; ///< error in egd
    std::vector<double> eb; ///< total energy released by delayed betas
    std::vector<double> deb; ///< error in deb
    std::vector<double> enu; ///< energy carried away by neutrinos
    std::vector<double> denu; ///< error in enu
    std::vector<double> er; ///< total energy release not lost to neutrinos
    std::vector<double> der; ///< error in er
    std::vector<double> et; ///< total energy release per fission
    std::vector<double> det; ///< error in et

  };

  /// reads ENDF fission energy release data from an ifstream and returns
  /// structure containing associated data
  /// \param infile An ENDF filestream
  mt458_mf1 read_mt458_mf1(std::ifstream &infile);

  /// data structure for delayed photon data
  class mt460_mf1 : public mt_base {
  public:

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

  };

  /// read delayed photon data from a filestream
  /// \param infile An ENDF filestream
  mt460_mf1 read_mt460_mf1(std::ifstream &infile);


  /// data structure for ENDF decay data
  class mt457_mf8 : public mt_base {
  public:
    int lis; ///< state of the original nuclide
    int liso; ///< isomeric state of the original nuclide
    int nst; ///< nucleus stability flag: 0 = radioactive
    int nsp; ///< total number of radiation types for which information is given
    std::vector<std::pair<double,double> > erel; ///< energy released
    std::vector<double> styp; ///< Decay radiation
    std::vector<int> lcon;
    std::vector<int> ner;
    std::vector<int> nd;
    std::vector<std::pair<double,double> > fd;
    std::vector<std::pair<double,double> > eav;
    std::vector<std::pair<double,double> > fc;
    std::vector<std::pair<double,double> > er;
    std::vector<double> rtyp;
    std::vector<double> type;
    std::vector<std::pair<double,double> > ri;
    std::vector<std::pair<double,double> > ris;
    std::vector<std::pair<double,double> > ricc;
    std::vector<std::pair<double,double> > rick;
    std::vector<std::pair<double,double> > ricl;
  };

  /// read decay data from a filestream
  /// \param infile An ENDF filestream
  mt457_mf8 read_mt457_mf8(std::ifstream &infile);




}

#endif
