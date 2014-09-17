#define NUM_PARTICLES = 32

// defines the primary particle types that are allowed by most monte carlo codes
// some monte carlo codes allow us to score so called "heavy ions", in fact we
// define heavy ions to be particles with more than one neutron or proton

namespace pyne 
{
namespace particle
{
  extern std::string _names[NUM_PARTICLES];
  extern std::set<std::string> names;
  extern std::map<std::string,int> name_id;
  extern std::map<int,std::string> id_name;
  extern std::map<std::string,std::string> docs;
  extern std::map<std::string,int> altnames;


  class NotAParticle : public std::exception
  {
  public:
    /// Default constructor
    NotAParticle () {};

    /// Default destructor
    ~NotAParticle () throw () {};

    NotAParticle(std::string particle_name)
      {
	part_name = particle_name;
      }

    virtual const char* what() const throw()
    {
      std::string pname ("Not a valid particle name ");
      if(!part_name.empty())
	pname += part_name;
      return (const char *) pname.c_str();
    };

    private:
      std::string part_name;  


  };
};
};
  
