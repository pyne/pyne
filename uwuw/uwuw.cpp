#include <unistd.h>
#include "uwuw.hpp"

// Empty Constructor
UWUW::UWUW()
{
};

// Default constructor
UWUW::UWUW(char* file) 
{
  std::string filename(file);
  // turn the filename into a full filepath
  full_filepath = get_full_filepath(filename);

  std::cout << full_filepath << std::endl;
  // load materials
  material_library = load_pyne_materials(full_filepath);
  // load tallies
  tally_library = load_pyne_tallies(full_filepath);
};

// Default constructor
UWUW::UWUW(std::string filename) 
{
  // turn the filename into a full filepath
  full_filepath = get_full_filepath(filename);
  // load materials
  material_library = load_pyne_materials(full_filepath);
  // load tallies
  tally_library = load_pyne_tallies(full_filepath);
  
};

// Destructor
UWUW::~UWUW() 
{
};

// convert convert a filename into path+filename (for pyne)
std::string UWUW::get_full_filepath(char *filename)
{
  std::string file(filename);
  return UWUW::get_full_filepath(file);
}

// convert convert a filename into path+filename (for pyne)
std::string UWUW::get_full_filepath(std::string filename)
{
  char *current_dir = get_current_dir_name(); //current working dir needed for pyne load 
                                        
  std::string cwd(current_dir), file(filename);   // pyne needs absolute path

  std::string full_filepath = cwd + '/' + file;

  // get rid of all trailing white space
  full_filepath.erase(std::remove_if( full_filepath.begin(),
				      full_filepath.end(), ::isspace),
                                      full_filepath.end());    
  return full_filepath;
}

// loads all materials into map
std::map<std::string, pyne::Material> UWUW::load_pyne_materials(std::string filename) 
{
  std::map<std::string, pyne::Material> library; // material library
  
  bool end = false; // end of materials
  int i = -1;

  // neednt check for filename existance, since it is guarenteed to exist
  while( !end )
    {
      pyne::Material mat; // from file

      mat.from_hdf5(filename,"/materials",++i);
      // if already in the map we have looped through the materials
      // and need not continue
      if ( library.count(mat.metadata["name"].asString()) )
	{
	  end = true;  
	}
      else
	{
	  // renumber material number by position in the library
	  mat.metadata["mat_number"]=i+1;
	  library[mat.metadata["name"].asString()]=mat;
	}
    }

  // remove the first entry in the material library, its an artefact from the hdf5 import
  library.erase(library.begin());
  
  std::cout << "Materials present in the h5m file" << std::endl;
  for(std::map<std::string,pyne::Material>::const_iterator it = library.begin() ; it != library.end() ; ++it )
    {
      std::cout << it->first <<  std::endl;
    }
  
  return library;
}

// loads all tallies into map
std::map<std::string, pyne::Tally> UWUW::load_pyne_tallies(std::string filename) 
{
  std::map<std::string, pyne::Tally> library; // material library
  
  bool end = false; // end of materials
  int i = -1;

  // neednt check for filename existance, since it is guarenteed to exist
  while( !end )
    {
      pyne::Tally tally; // from file

      tally.from_hdf5(filename,"/tally",++i);
      // if already in the map we have looped through the materials
      // and need not continue
      if ( library.count(tally.tally_name) )
	{
	  end = true;  
	}

      library[tally.tally_name]=tally;
    }
  
  
  std::cout << "Tallies present in the h5m file" << std::endl;
  for(std::map<std::string,pyne::Tally>::const_iterator it = library.begin() ; it != library.end() ; ++it )
    {
      std::cout << it->first <<  std::endl;
    }

  return library;
}

