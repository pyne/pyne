// Header for general library file.

#if !defined(_PYNE_)
#define _PYNE_

//standard libraries
#include <string>
#include <string.h>
#include <sstream>
#include <cctype>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <exception>
#include <sys/stat.h> 
#include <vector>
#include <algorithm>

/*** Macros ***/
#define length_array(a) ( sizeof ( a ) / sizeof ( *a ) )

#if defined(__APPLE__)
#if (__GNUC__ >= 4)
  #include <cmath>
  #define isnan(x) std::isnan(x)
#else
  #include <math.h>
  #define isnan(x) __isnand((double)x)
#endif
#endif

#ifdef _WIN32_MSVC
    #define isnan(x) ((x) != (x))
#endif

namespace pyne {

  // PyNE Globals
  void pyne_start ();

  extern std::string PYNE_DATA;
  extern std::string NUC_DATA_PATH;

  // String Transformations
  static std::string digits = "0123456789";
  static std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  static std::string words = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_";

  std::string to_str (int);
  std::string to_str (double);
  std::string to_str (bool);

  int to_int (std::string);

  double to_dbl (std::string);

  std::string to_upper(std::string);

  std::string to_lower(std::string);

  std::string get_flag(char [], int);

  std::string remove_substring(std::string, std::string);

  std::string remove_characters(std::string, std::string);

  std::string replace_all_substrings(std::string, std::string, std::string);

  std::string last_char(std::string);

  std::string slice_from_end(std::string, int = -1, int = 1);

  bool ternary_ge(int, int, int);

  bool contains_substring(std::string, std::string);

  std::string natural_naming(std::string);


  // Math Helpers
  extern const double pi;
  extern const double N_A;	//Avagardo's Number
  extern const double barns_per_cm2; 	// barns per cm^2
  extern const double cm2_per_barn; 	// cm^2 per barn
  extern const double sec_per_day; 	// seconds per day

  double slope (double, double, double, double);
  double solve_line (double, double, double, double, double);

  double tanh(double);
  double coth(double);


  // File Helpers
  bool file_exists(std::string); 



  /******************/
  /*** Exceptions ***/
  /******************/

  class FileNotFound : public std::exception
  {
  public:
    FileNotFound () {};

    ~FileNotFound () throw () {};

    FileNotFound(std::string fname)
    {
      filename = fname;
    };

    virtual const char* what() const throw()
    {
      std::string FNFstr ("File not found: ");
      if (!filename.empty())
        FNFstr += filename;

      return (const char *) FNFstr.c_str();
    };

  private:
    std::string filename;
  };


// End PyNE namespace
};

#endif
