/// \brief This is the base PyNE library.
///
/// It contains a lot of utility functions and constants that are globaly useful
/// through out the rest of the PyNE infrastructure.
///

// Header for general library file.

#ifndef PYNE_KMMHYNANYFF5BFMEYIP7TUNLHA
#define PYNE_KMMHYNANYFF5BFMEYIP7TUNLHA

//standard libraries
#include <string>
#include <string.h>
#include <sstream>
#include <cctype>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <vector>
#include <algorithm>

#if (__GNUC__ >= 4)
  #include <cmath>
  #define isnan(x) std::isnan(x)
#else
  #include <math.h>
  #define isnan(x) __isnand((double)x)
#endif

#ifdef _WIN32
#define isnan(x) ((x) != (x))
#endif

#ifndef JSON_IS_AMALGAMATION
  #define JSON_IS_AMALGAMATION
#endif

/// The 'pyne' namespace all PyNE functionality is included in.
namespace pyne {

  void pyne_start (); ///< Initializes PyNE based on environment.

  /// Path to the directory containing the PyNE data.
  extern std::string PYNE_DATA;
  extern std::string NUC_DATA_PATH; ///< Path to the nuc_data.h5 file.
  extern std::string VERSION; ///< PyNE version number

  // String Transformations
  /// string of digit characters
  static std::string digits = "0123456789";
  /// uppercase alphabetical characters
  static std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  /// string of all valid word characters for variable names in programing languages.
  static std::string words = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_";

  /// \name String Conversion Functions
  /// \{
  /// Converts the variables of various types to their C++ string representation.
  std::string to_str(int t);
  std::string to_str(unsigned int t);
  std::string to_str(double t);
  std::string to_str(bool t);
  /// \}

  int to_int(std::string s);  ///< Converts a string of digits to an int using atoi().

  double to_dbl(std::string s);  ///< Converts a valid string to a float using atof().

  /// Converts a string from ENDF format to a float. Only handles E-less format
  /// but is roughly 5 times faster than endftod.
  double endftod_cpp(char * s);
  double endftod_f(char * s); ///< Converts a string from ENDF format to a float.
  extern  double (*endftod)(char * s); ///< endftod function pointer. defaults to fortran

  void use_fast_endftod();/// switches endftod to fast cpp version

  /// Returns an all upper case copy of the string.
  std::string to_upper(std::string s);

  /// Returns an all lower case copy of the string.
  std::string to_lower(std::string s);

  /// Returns a capitalized copy of the string.
  std::string capitalize(std::string s);

  /// Forms and returns the wrapped lines with a lenght up to line_lenght.
  std::string comment_line_wrapping(std::string line, std::string comment_prefix = "",
                                           int line_length = 79);

  /// Finds and returns the first white-space delimited token of a line.
  /// \param line a character array to take the first token from.
  /// \param max_l an upper bound to the length of the token.  Must be 11 or less.
  /// \returns a the flag as a string
  std::string get_flag(char line[], int max_l);

  /// Creates a copy of \a s with all instances of \a substr taken out.
  std::string remove_substring(std::string s, std::string substr);

  /// Removes all characters in the string \a chars from \a s.
  std::string remove_characters(std::string s, std::string chars);

  /// Replaces all instance of \a substr in \a s with \a repstr.
  std::string replace_all_substrings(std::string s, std::string substr,
                                                    std::string repstr);

  /// Returns the last character in a string.
  std::string last_char(std::string s);

  /// Returns the slice of a string \a s using the negative index \a n and the
  /// length of the slice \a l.
  std::string slice_from_end(std::string s, int n=-1, int l=1);

  /// Returns true if \a a <= \a b <= \a c and flase otherwise.
  bool ternary_ge(int a, int b, int c);

  /// Returns true if \a substr is in \a s.
  bool contains_substring(std::string s, std::string substr);

  /// Calculates a version of the string \a name that is also a valid variable name.
  /// That is to say that the return value uses only word characters.
  std::string natural_naming(std::string name);

  // split a string into a vector of string using a delimiter
  std::vector<std::string> split_string(std::string lists, std::string delimiter = " ");

  // join the vector element into a string, each values will be delimited ny the delimiter
  template<typename T>
  std::string join_to_string(std::vector<T> vect, std::string delimiter = " "){
    std::stringstream out;
    out << std::setiosflags(std::ios::fixed) << std::setprecision(6);

    // ensure there is at least 1 element in the vector
    if (vect.size() == 0)
      return out.str();
    // no delimiter needed before the first element
    out << vect[0];
    for( int i = 1; i < vect.size(); i++)
      out << delimiter << vect[i];
    return out.str();
  };
  template std::string join_to_string(std::vector<int> vect,
                                      std::string delimiter);
  template std::string join_to_string(std::vector<double> vect,
                                      std::string delimiter);
  template std::string join_to_string(std::vector<std::string> vect, std::string delimiter);

  /// Finds the slope of a line from the points (\a x1, \a y1) and (\a x2, \a y2).
  double slope (double x2, double y2, double x1, double y1);

  /// Solves the equation for the line y = mx + b, given \a x and the points that
  /// form the line: (\a x1, \a y1) and (\a x2, \a y2).
  double solve_line (double x, double x2, double y2, double x1, double y1);

  double tanh(double x);  ///< The hyperbolic tangent function.
  double coth(double x);  ///< The hyperbolic cotangent function.


  // File Helpers
  /// Returns true if the file can be found.
  bool file_exists(std::string strfilename);

  // turns the filename string into the full file path
  std::string get_full_filepath(char* filename);
  // turns the filename string into the full file path
  std::string get_full_filepath(std::string filename);

  // Message Helpers
  extern bool USE_WARNINGS;
  /// Toggles warnings on and off
  bool toggle_warnings();

  /// Prints a warning message.
  void warning(std::string s);

  /// Custom exception to be thrown in the event that a required file is not able to
  /// be found.
  class FileNotFound : public std::exception
  {
  public:

    /// default constructor
    FileNotFound () {};

    /// default destructor
    ~FileNotFound () throw () {};

    /// constructor with the filename \a fname.
    FileNotFound(std::string fname)
    {
      FNF_message = "File not found";
      if (!fname.empty())
        FNF_message += ": " + fname;
    };

    /// Creates a helpful error message.
    virtual const char* what() const throw()
    {
      return FNF_message.c_str();
    };

  private:
    std::string FNF_message; /// Message for exception
  };

  /// Exception representing value errors of all kinds
  class ValueError : public std::exception
  {
  public:

    /// default constructor
    ValueError () {};

    /// default destructor
    ~ValueError () throw () {};

    /// constructor with the filename \a fname.
    ValueError(std::string msg)
    {
      message = "ValueError: " + msg;
    };

    /// Creates a helpful error message.
    virtual const char* what() const throw()
    {
      const char* msgstr_rtn = message.c_str();
      return msgstr_rtn;
    };

  protected:
    std::string message; ///< extra message for the user.
  };


// End PyNE namespace
}

#endif  // PYNE_KMMHYNANYFF5BFMEYIP7TUNLHA
