// General Library
#ifndef PYNE_IS_AMALGAMATED
extern "C" double endftod_(char *str, int len);
#endif
#include <iomanip>

#ifdef _WIN32
  #include <filesystem>
#endif


#ifndef PYNE_IS_AMALGAMATED
#include "utils.h"
#include "pyne_version.h"
#endif



// PyNE Globals

std::string pyne::PYNE_DATA = "";
std::string pyne::NUC_DATA_PATH = "";
std::string pyne::VERSION = PYNE_VERSION_STRING;

void pyne::pyne_start() {
#if defined _WIN32
  char * tmpPYNE_DATA;
  size_t lenPYNE_DATA;
  errno_t errPYNE_DATA = _dupenv_s(&tmpPYNE_DATA, &lenPYNE_DATA, "PYNE_DATA");
  if (errPYNE_DATA)
    tmpPYNE_DATA = (char *) "<NOT_FOUND>";
  PYNE_DATA = (std::string) tmpPYNE_DATA;

  char * tmpNUC_DATA_PATH;
  size_t lenNUC_DATA_PATH;
  errno_t errNUC_DATA_PATH = _dupenv_s(&tmpNUC_DATA_PATH, &lenNUC_DATA_PATH, "NUC_DATA_PATH");
  if (errPYNE_DATA)
    tmpNUC_DATA_PATH = (char *) "<NOT_FOUND>";
  NUC_DATA_PATH = (std::string) tmpNUC_DATA_PATH;
#else
  char * tmppath;
  tmppath = getenv("PYNE_DATA");
  if (tmppath == NULL)
      tmppath = (char *) "<NOT_FOUND>";
  PYNE_DATA = std::string(tmppath);

  tmppath = getenv("NUC_DATA_PATH");
  if (tmppath == NULL)
      tmppath = (char *) "<NOT_FOUND>";
  NUC_DATA_PATH = std::string(tmppath);
#endif
  return;
}



// String Transformations
std::string pyne::to_str(int t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(unsigned int t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(double t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(bool t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}


int pyne::to_int(std::string s) {
  return atoi( s.c_str() );
}

double pyne::to_dbl(std::string s) {
  return strtod( s.c_str(), NULL );
}

double pyne::endftod_cpp(char * s) {
  // Converts string from ENDF only handles "E-less" format but is 5x faster
  int pos, mant, exp;
  double v, dbl_exp;

  mant = exp = 0;
  if (s[2] == '.') {
    // Convert an ENDF float
    if (s[9] == '+' || s[9] == '-') {
      // All these factors of ten are from place values.
      mant = s[8] + 10 * s[7] + 100 * s[6] + 1000 * s[5] + 10000 * s[4] + \
             100000 * s[3] + 1000000 * s[1] - 1111111 * '0';
      exp = s[10] - '0';
      // Make the right power of 10.
      dbl_exp = exp & 01? 10.: 1;
      dbl_exp *= (exp >>= 1) & 01? 100.: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e4: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e8: 1;
      // Adjust for powers of ten from treating mantissa as an integer.
      dbl_exp = (s[9] == '-'? 1/dbl_exp: dbl_exp) * 1.0e-6;
      // Get mantissa sign, apply exponent.
      v = mant * (s[0] == '-'? -1: 1) * dbl_exp;
    } else {
      mant = s[7] + 10 * s[6] + 100 * s[5] + 1000 * s[4] + 10000 * s[3] + \
             100000 * s[1] - 111111 * '0';
      exp = s[10] + 10 * s[9] - 11 * '0';
      dbl_exp = exp & 01? 10.: 1;
      dbl_exp *= (exp >>= 1) & 01? 100.: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e4: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e8: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e16: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e32: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e64: 1;
      dbl_exp = (s[8] == '-'? 1/dbl_exp: dbl_exp) * 1.0e-5;
      v = mant * (s[0] == '-'? -1: 1) * dbl_exp;
    }
  }

  // Convert an ENDF int to float; we start from the last char in the field and
  // move forward until we hit a non-digit.
  else {
    v = 0;
    mant = 1; // Here we use mant for the place value about to be read in.
    pos = 10;
    while (s[pos] != '-' && s[pos] != '+' && s[pos] != ' ' && pos > 0) {
      v += mant * (s[pos] - '0');
      mant *= 10;
      pos--;
    }
    v *= (s[pos] == '-'? -1: 1);
  }
  return v;
}

double pyne::endftod_f(char * s) {
#ifdef PYNE_IS_AMALGAMATED
  return endftod_cpp(s);
#else
  return endftod_(s, 12);
#endif
}

double (*pyne::endftod)(char * s) = &pyne::endftod_f;

void pyne::use_fast_endftod() {
  pyne::endftod = &pyne::endftod_cpp;
}

std::string pyne::to_upper(std::string s) {
  // change each element of the string to upper case.
  for(unsigned int i = 0; i < s.length(); i++)
    s[i] = toupper(s[i]);
  return s;
}

std::string pyne::to_lower(std::string s) {
  // change each element of the string to lower case
  for(unsigned int i = 0; i < s.length(); i++)
    s[i] = tolower(s[i]);
  return s;
}

std::string pyne::comment_line_wrapping(std::string line,
                                               std::string comment_prefix,
                                               int line_length) {
  std::ostringstream oss;

  line_length -= comment_prefix.length();
    
  // Include as is if short enough
  while (line.length() > line_length) {
    oss << comment_prefix << line.substr(0, line_length) << std::endl;
    line.erase(0, line_length);
  }

  if (line.length() > 0) {
    oss << comment_prefix << line << std::endl;
  }

  return oss.str();
}

std::string pyne::capitalize(std::string s) {
  unsigned int slen = s.length();
  if (slen == 0)
    return s;
  // uppercase the first character
  s[0] = toupper(s[0]);
  // change each subsequent element of the string to lower case
  for(unsigned int i = 1; i < slen; i++)
    s[i] = tolower(s[i]);
  return s;
}


std::string pyne::get_flag(char line[], int max_l) {
  char tempflag [10];
  for (int i = 0; i < max_l; i++)
  {
    if (line[i] == '\t' || line[i] == '\n' || line[i] == ' ' || line[i] == '\0')
    {
      tempflag[i] = '\0';
      break;
    }
    else
      tempflag[i] = line[i];
  }
  return std::string (tempflag);
}



std::string pyne::remove_substring(std::string s, std::string substr) {
  // Removes a substring from the string s
  int n_found = s.find(substr);
  while ( 0 <= n_found ) {
    s.erase( n_found , substr.length() );
    n_found = s.find(substr);
  }
  return s;
}


std::string pyne::remove_characters(std::string s, std::string chars) {
  // Removes all characters in the string chars from the string s
  for (int i = 0; i < chars.length(); i++ ) {
    s = remove_substring(s, chars.substr(i, 1) );
  }
  return s;
}


std::string pyne::replace_all_substrings(std::string s, std::string substr, std::string repstr) {
  // Replaces all instance of substr in s with the string repstr
  int n_found = s.find(substr);
  while ( 0 <= n_found ) {
    s.replace( n_found , substr.length(), repstr );
    n_found = s.find(substr);
  }
  return s;
}



std::string pyne::last_char(std::string s) {
    // Returns the last character in a string.
    return s.substr(s.length()-1, 1);
}


std::string pyne::slice_from_end(std::string s, int n, int l) {
  // Returns the slice of a string using negative indices.
  return s.substr(s.length()+n, l);
}


bool pyne::ternary_ge(int a, int b, int c) {
  // Returns true id a <= b <= c and flase otherwise.
  return (a <= b && b <= c);
}


bool pyne::contains_substring(std::string s, std::string substr) {
  // Returns a boolean based on if the sub is in s.
  int n = s.find(substr);
  return ( 0 <= n && n < s.length() );
}


std::string pyne::natural_naming(std::string name) {
  // Calculates a version on the string name that is a valid
  // variable name, ie it uses only word characters.
  std::string nat_name (name);

  // Replace Whitespace characters with underscores
  nat_name = pyne::replace_all_substrings(nat_name, " ",  "_");
  nat_name = pyne::replace_all_substrings(nat_name, "\t", "_");
  nat_name = pyne::replace_all_substrings(nat_name, "\n", "_");

  // Remove non-word characters
  int n = 0;
  while ( n < nat_name.length() ) {
    if ( pyne::words.find(nat_name[n]) == std::string::npos )
      nat_name.erase(n, 1);
    else
      n++;
  }

  // Make sure that the name in non-empty before continuing
  if (nat_name.length() == 0)
    return nat_name;

  // Make sure that the name doesn't begin with a number.
  if ( pyne::digits.find(nat_name[0]) != std::string::npos)
    nat_name.insert(0, "_");

  return nat_name;
}


std::vector<std::string> pyne::split_string(std::string particles_list, std::string delimiter) {
  std::vector<std::string> output_vector;
  size_t prev_pos = 0; //item start position
  size_t pos = 0; //item end position
 
  while( (pos = particles_list.find(delimiter, prev_pos)) != std::string::npos){
    output_vector.push_back(particles_list.substr(prev_pos, pos));
    prev_pos = pos + delimiter.length();
  }
  // catch list with a single particle
  if (pos == std::string::npos && prev_pos == 0 && particles_list.length() >0)
    output_vector.push_back(particles_list);

  return output_vector;
}


//
// Math Helpers
//

double pyne::slope(double x2, double y2, double x1, double y1) {
  // Finds the slope of a line.
  return (y2 - y1) / (x2 - x1);
}


double pyne::solve_line(double x, double x2, double y2, double x1, double y1) {
  return (slope(x2,y2,x1,y1) * (x - x2)) + y2;
}


double pyne::tanh(double x) {
  return std::tanh(x);
}

double pyne::coth(double x) {
  return 1.0 / std::tanh(x);
}



// File Helpers

bool pyne::file_exists(std::string strfilename) {
  // Thank you intarwebz for this function!
  // Sepcifically: http://www.techbytes.ca/techbyte103.html
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strfilename.c_str(), &stFileInfo);

  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  }
  else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }

  return(blnReturn);
}

// convert convert a filename into path+filename (for pyne)
std::string pyne::get_full_filepath(char* filename) {
  std::string file(filename);
  return pyne::get_full_filepath(file);
}

// convert convert a filename into path+filename (for pyne)
std::string pyne::get_full_filepath(std::string filename) {
  // remove all extra whitespace
  filename = pyne::remove_characters(" " , filename);
  // use stdlib call
#ifndef _WIN32
  const std::string full_filepath = realpath(filename.c_str(), NULL);
#else
  const std::string  full_filepath = std::filesystem::canonical(filename.c_str()).string();
#endif
  return full_filepath;
}

// Message Helpers

bool pyne::USE_WARNINGS = true;

bool pyne::toggle_warnings(){
  USE_WARNINGS = !USE_WARNINGS;
  return USE_WARNINGS;
}

void pyne::warning(std::string s){
  // Prints a warning message
  if (USE_WARNINGS){
    std::cout << "\033[1;33m WARNING: \033[0m" << s << "\n";
  }
}




