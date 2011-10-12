// General Library 

#include "pyne.h"


// PyNE Globals

std::string pyne::PYNE_DATA = "";
std::string pyne::NUC_DATA_PATH = "";

void pyne::pyne_start()
{
#ifdef _WIN32
  char * tmpPYNE_DATA;
  size_t lenPYNE_DATA;
  errno_t errPYNE_DATA = _dupenv_s(&tmpPYNE_DATA, &lenPYNE_DATA, "PYNE_DATA");
  if (errPYNE_DATA)
    std::cout << "PYNE_DATA enviromental variable could not be found\n";
  PYNE_DATA = (std::string) tmpPYNE_DATA;

  char * tmpNUC_DATA_PATH;
  size_t lenNUC_DATA_PATH;
  errno_t errNUC_DATA_PATH = _dupenv_s(&tmpNUC_DATA_PATH, &lenNUC_DATA_PATH, "NUC_DATA_PATH");
  if (errPYNE_DATA)
    std::cout << "NUC_DATA_PATH enviromental variable could not be found\n";
  NUC_DATA_PATH = (std::string) tmpNUC_DATA_PATH;
#else
  PYNE_DATA = getenv("PYNE_DATA");
  NUC_DATA_PATH = getenv("NUC_DATA_PATH");
#endif
  return;
};



// String Transformations
std::string pyne::to_str (int t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str (double t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str (bool t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}


int pyne::to_int (std::string s)
{
  return atoi( s.c_str() );
}

double pyne::to_dbl (std::string s)
{
  return strtod( s.c_str(), NULL );
}



std::string pyne::to_upper(std::string s)
{
  // change each element of the string to upper case.
  for(unsigned int i = 0; i < s.length(); i++)
  {
    s[i] = toupper(s[i]);
  }
  return s;
}

std::string pyne::to_lower(std::string s)
{
  // change each element of the string to lower case
  for(unsigned int i = 0; i < s.length(); i++)
  {
    s[i] = tolower(s[i]);
  }
  return s;
}


std::string pyne::get_flag(char line[], int max_l)
{
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



std::string pyne::remove_substring(std::string s, std::string substr)
{
  // Removes a substring from the string s
  int n_found = s.find(substr);
  while ( 0 <= n_found )
  {
    s.erase( n_found , substr.length() );
    n_found = s.find(substr);
  }
  return s;
}


std::string pyne::remove_characters(std::string s, std::string chars)
{
  // Removes all characters in the string chars from thr string s
  for (int i = 0; i < chars.length(); i++ )
  {
    s = remove_substring(s, chars.substr(i, 1) );
  }
  return s;
}


std::string pyne::replace_all_substrings(std::string s, std::string substr, std::string repstr)
{
  // Replaces all instance of substr in s with the string repstr
  int n_found = s.find(substr);
  while ( 0 <= n_found )
  {
    s.replace( n_found , substr.length(), repstr );
    n_found = s.find(substr);
  }
  return s;
};



std::string pyne::last_char(std::string s)
{
    // Returns the last character in a string.
    return s.substr(s.length()-1, 1);
}


std::string pyne::slice_from_end(std::string s, int n, int l)
{
  // Returns the slice of a string using negative indices.
  return s.substr(s.length()+n, l);
}


bool pyne::ternary_ge(int a, int b, int c)
{
  // Returns true id a <= b <= c and flase otherwise.
  return (a <= b && b <= c); 
}


bool pyne::contains_substring(std::string s, std::string substr)
{
  // Returns a boolean based on if the sub is in s.
  int n = s.find(substr);
  return ( 0 <= n && n < s.length() );
}


std::string pyne::natural_naming(std::string name)
{
  // Calculates a version on the string name that is a valid
  // variable name, ie it uses only word characters.
  std::string nat_name (name);

  // Replace Whitespace characters with underscores
  nat_name = pyne::replace_all_substrings(nat_name, " ",  "_");
  nat_name = pyne::replace_all_substrings(nat_name, "\t", "_");
  nat_name = pyne::replace_all_substrings(nat_name, "\n", "_");

  // Remove non-word characters
  int n = 0;
  while ( n < nat_name.length() )
  {
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
};


//
// Math Helpers
//

const double pyne::pi = 3.14159265359;
const double pyne::N_A = 6.0221415 * pow(10.0, 23);
const double pyne::barns_per_cm2 = pow(10.0, 24); 
const double pyne::cm2_per_barn = pow(10.0, -24); 
const double pyne::sec_per_day = 24.0 * 3600.0; 

double pyne::slope(double x2, double y2, double x1, double y1)
{
  // Finds the slope of a line.
  return (y2 - y1) / (x2 - x1);
};


double pyne::solve_line(double x, double x2, double y2, double x1, double y1)
{
  return (slope(x2,y2,x1,y1) * (x - x2)) + y2;
};


double pyne::tanh(double x)
{
  return std::tanh(x);
};

double pyne::coth(double x)
{
  return 1.0 / std::tanh(x);
};




// File Helpers

bool pyne::file_exists(std::string strfilename) 
{
  // Thank you intarwebz for this function!
  // Sepcifically: http://www.techbytes.ca/techbyte103.html
  struct stat stFileInfo; 
  bool blnReturn; 
  int intStat; 

  // Attempt to get the file attributes 
  intStat = stat(strfilename.c_str(), &stFileInfo); 

  if(intStat == 0) 
  { 
    // We were able to get the file attributes 
    // so the file obviously exists. 
    blnReturn = true; 
  } 
  else 
  { 
    // We were not able to get the file attributes. 
    // This may mean that we don't have permission to 
    // access the folder which contains this file. If you 
    // need to do that level of checking, lookup the 
    // return values of stat which will give you 
    // more details on why stat failed. 
    blnReturn = false; 
  } 
   
  return(blnReturn); 
};
