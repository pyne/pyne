#include "material.h"
#define  PYNE_MAT "/filespace/people/z/zachman/.local/lib/python2.7/site-packages/pyne/nuc_data.h5"

int main (int argc, char* argv[])
{
   using namespace pyne;

   Material test_mat;
   std::string filename;
   filename = "/filespace/people/z/zachman/.local/lib/python2.7/site-packages/pyne/nuc_data.h5";

   // Week of May 12: Define a comp_map by hand in a valid way
   comp_map nucvec;
   nucvec.insert(std::pair<int,double>(10010000, 1.0)); 
   nucvec.insert(std::pair<int,double>(80160000, 1.0));   
   nucvec.insert(std::pair<int,double>(69169000, 1.0));   
   nucvec.insert(std::pair<int,double>(92235000, 1.0));   
   nucvec.insert(std::pair<int,double>(92238000, 1.0));   
   nucvec.insert(std::pair<int,double>(94239000, 1.0));   
   nucvec.insert(std::pair<int,double>(94241000, 1.0));   
   nucvec.insert(std::pair<int,double>(95242000, 1.0));   
   nucvec.insert(std::pair<int,double>(96244000, 1.0));   

   comp_map::iterator it;

   filename = PYNE_MAT;
   // filename="output.h5m";
   test_mat = Material(nucvec);

   // make an entir db from the stuff in the file
   // test_mat.from_hdf5(filename);
   // test_mat.from_hdf5(filename,"", 0,1);
   // This may populate the object with .attrs, 
   // test_mat.from_hdf5(filename,"/material", 0,1);

   // std::cout << "density = " << test_mat.density << std::endl;
   // std::cout << "mass = " << test_mat.mass << std::endl;
   // std::cout << "atoms_per_mol = " << test_mat.atoms_per_mol << std::endl;


   for (it=test_mat.comp.begin(); it != test_mat.comp.end(); ++it)
   {
       std::cout << it->first << ", " << it->second << std::endl;
   }
   // std::cout << test_mat.write_mcnp
   // std::cout << test_mat.metadata << std::endl;

}
