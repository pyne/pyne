#include "pyne/pyne.h"
#include <map>
#include <string>

class UWUW 
{
public:
  UWUW(); // empty constructor

  UWUW(char* filename); // c style filename

  UWUW(std::string filename); // normal constructor
 
  ~UWUW(); // destructor

  std::map<std::string, pyne::Material> material_library; // material library
  std::map<std::string, pyne::Tally> tally_library; // tally library
  std::string full_filepath;
  
private:
  // turns the filename string into the full file path
  std::string get_full_filepath(char *filename);
  // turns the filename string into the full file path
  std::string get_full_filepath(std::string filename);

  std::map<std::string, pyne::Material> load_pyne_materials(std::string filename);
  std::map<std::string, pyne::Tally> load_pyne_tallies(std::string filename);

};
