// General Library

#include "pyne.h"


// PyNE Globals

std::string pyne::PYNE_DATA = "";
std::string pyne::NUC_DATA_PATH = "";

void pyne::pyne_start()
{
#if defined __WIN_MSVC__
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
};



// String Transformations
std::string pyne::to_str (int t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(unsigned int t)
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

double pyne::endftod (char * s)
{
  // Converts string from ENDF to float64.
  int pos, mant, exp, mantsign, expsign;
  double v, dbl_exp;

  mant = exp = 0;
  mantsign = expsign = 1;
  if (s[2] == '.')
  // Convert an ENDF float
    {
      // if (s[0] == '-')
      //   mantsign = -1;
      if (s[9] == '+' or s[9] == '-')
        {
          // if (s[9] == '-')
          //   expsign = -1;
          // mantsize = 8;
          mant = s[8] + 10 * s[7] + 100 * s[6] + 1000 * s[5] + 10000 * s[4] + \
            100000 * s[3] + 1000000 * s[1] - 1111111 * '0';
          exp = s[10] - '0';
          dbl_exp = (exp & 01? 10.: 1) * \
                    ((exp >>= 1) & 01? 100.: 1) * \
                    ((exp >>= 1) & 01? 1.0e4: 1) * \
                    ((exp >>= 1) & 01? 1.0e8: 1) * \
                    ((exp >>= 1) & 01? 1.0e16: 1) * \
                    ((exp >>= 1) & 01? 1.0e32: 1) * \
                    ((exp >>= 1) & 01? 1.0e64: 1);
          dbl_exp = (s[9] == '-'? 1/dbl_exp: dbl_exp) * 1.0e-6;
          // 3.12345+6 = 3123450 = 312345 * 10^(6-5)
          // 3.12345-1 = 0.312345 = 312345 * 10^-(1+5)
          // 3.12345+4 = 31234.5 = 312345 * 10^(4-5)
          // exp = (s[9] == '-'? -1: 1) * (s[10] - '0') - 6;
          v = mant * (s[0] == '-'? -1: 1) * dbl_exp;
        }
      else
        {
          if (s[8] == '-')
            expsign = -1;
          // mantsize = 7;
          mant = s[7] + 10 * s[6] + 100 * s[5] + 1000 * s[4] + 10000 * s[3] + \
            100000 * s[1] - 111111 * '0';
          // exp = (s[8] == '-'? -1: 1) * (s[10] + 10 * s[9] - 11 * '0') - 5;
          exp = s[10] + 10 * s[9] - 11 * '0';
          dbl_exp = (exp & 01? 10.: 1) * \
                    ((exp >>= 1) & 01? 100.: 1) * \
                    ((exp >>= 1) & 01? 1.0e4: 1) * \
                    ((exp >>= 1) & 01? 1.0e8: 1) * \
                    ((exp >>= 1) & 01? 1.0e16: 1) * \
                    ((exp >>= 1) & 01? 1.0e32: 1) * \
                    ((exp >>= 1) & 01? 1.0e64: 1);
          dbl_exp = (s[8] == '-'? 1/dbl_exp: dbl_exp) * 1.0e-5;
          // dbl_exp *= s[9] == '-'? 1: 1.0e-5;
          // v = mant * (s[0] == '-'? -1: 1) * ((exp << 6) < 0? 1/dbl_exp: dbl_exp);
          v = mant * (s[0] == '-'? -1: 1) * dbl_exp;
        }
      // for (pos = 1; pos <= mantsize; pos++)
      //   {
      //     c = s[pos];
      //     if (c != '.')
      //       {
      //         mant = 10*mant + (c - '0');
      //       }
      //   }
      // for (pos = mantsize + 2; pos <= 10; pos++) // +2 because of sign chars
      //   {
      //     c = s[pos];
      //     exp = 10*exp + (c - '0');
      //   }
      // if (expsign == -1)
      //   exp = -exp - (mantsize - 2); // -2 because of mant sign and '.' char
      // else
      //   exp = exp - (mantsize - 2);
      // if (exp < 0)
      //   {
      //     expsign = -1;
      //     exp *= -1;
      //   }
      // v = mant * mantsign;
      // v = mant * (s[0] == '-'? -1: 1);
      // dbl_exp = (exp & 01? 10.: 1) * \
      //           ((exp >>= 1) & 01? 100.: 1) * \
      //           ((exp >>= 1) & 01? 1.0e4: 1) * \
      //           ((exp >>= 1) & 01? 1.0e8: 1) * \
      //           ((exp >>= 1) & 01? 1.0e16: 1) * \
      //           ((exp >>= 1) & 01? 1.0e32: 1) * \
      //           ((exp >>= 1) & 01? 1.0e64: 1);
      // for (d = powers_of_10; exp != 0; exp >>= 1, d += 1)
      //   {
      //     if (exp & 01)
      //       dbl_exp *= *d;
      //   }
      // if (expsign == -1)
      //   v /= dbl_exp;
      // else
      //   v *= dbl_exp;
    }

  // Convert an ENDF int to float
  else
    {
      mant = s[10] - '0';
      pos = 9;
      while (s[pos] != '-' and s[pos] != '+' and s[pos] != ' ' and pos > 0)
        {
          mant = mant + 10*(s[pos] - '0');
          pos--;
        }
      if (s[pos] == '-')
        mantsign = -1;
      // pos = v = 0;
      // for (; pos <=10; pos++)
      //   {
      //     c = s[pos];
      //     if (c == ' ')
      //       {
      //       }
      //     else if (c == '-')
      //       mantsign = -1;
      //     else
      //       mant = 10*mant + (c - '0');
      //   }
      v = mant * mantsign;
    }
  return v;
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
  // Removes all characters in the string chars from the string s
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
const double pyne::N_A = 6.0221415e+23;
const double pyne::barns_per_cm2 = 1e24; 
const double pyne::cm2_per_barn = 1e-24; 
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
