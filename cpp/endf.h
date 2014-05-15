
#include <sstream>
#include <fstream>

#ifndef PYNE_IS_AMALGAMATED
#include "pyne.h"
#endif

namespace pyne
{
namespace endf 
{

  /// a struct matching the ENDF control structure.
  typedef struct control_struct {
    double c1;  ///< First value in control struct ZZZAAA.M for HEAD [double]
    double c2;  ///< Second value in control struct AWR for HEAD [double]
    int l1; ///< first int in control struct
    int l2; ///< second int in control struct
    int n1; ///< third int in control struct
    int n2; ///< fourth int in control struct
    int mat; ///< material id
    int mf; ///< material file number
    int mt; ///< reaction number
  } control_struct; 
  
  /// a struct mathing the ENDF LIST structure
  typedef struct list_struct {
    double c1;  ///< First value in control struct ZZZAAA.M for HEAD [double]
    double c2;  ///< Second value in control struct AWR for HEAD [double]
    int l1; ///< first int in control struct
    int l2; ///< second int in control struct
    int npl; ///< number of data points
    int n2; ///< fourth int in control struct
    int mat; ///< material id
    int mf; ///< material file number
    int mt; ///< reaction number
    std::vector<double> data; ///<data as doubles
  } list_struct;   
  
  
    typedef struct mt_base_struct {
        int nuc_id; ///< 
        double awr; ///< mass relative to neutron [double]
        int mat; ///< ENDF material
        int mf; ///< ENDF file number
        int mt; ///< ENDF reaction designation
    } mt_base_struct;
    
    typedef struct mt_451_struct : mt_base_struct {
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
        int nwd;
        int nxc;
        std::vector<std::vector<int> > mt_list; ///< list of [mf,mt,lines,mod]
    } mt_451_struct; 
    
    typedef struct mt_fpy_struct : mt_base_struct {
        int le; ///< number of energy  dependent yields given minus 1
        int nfp; ///< number of fission products listed
        std::vector<int> i; ///< interpolation to be used
        std::vector<double> e; ///< list of energies
        std::vector<std::vector<double> > yields; ///< yield data [zafp,fps,yi,dyi]
    } mt_fpy_struct;
    
    typedef struct endf_library_struct {
        mt_451_struct mt_451;
        mt_fpy_struct mt_454;
        mt_fpy_struct mt_459;
        
        //std::list<mt_base_struct*> struct_list;
        ~endf_library_struct();
    } endf_library_struct;
    
    //std::vector<std::vector<int> > get_library_contents;
    
    //load_dataset_to_API(int mat, int mf, int mt)
    
    void read_endf(std::string filenm);
    void read_459(std::ifstream &infile);
    void read_454(std::ifstream &infile);
    void read_451(std::ifstream &infile);
    list_struct read_list(std::ifstream &infile);
    control_struct read_cont(std::ifstream &infile);
    
    extern endf_library_struct library;
}
}