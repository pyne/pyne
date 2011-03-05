// General Library 

#include "bright.h"

//Bright Globals

std::string bright::BRIGHT_DATA = "";

void bright::bright_start()
{
#ifdef _WIN32
    char * tmpBRIGHT_DATA;
    size_t lenBRIGHT_DATA;
    errno_t errBRIGHT_DATA = _dupenv_s(&tmpBRIGHT_DATA, &lenBRIGHT_DATA, "BRIGHT_DATA");
    if (errBRIGHT_DATA) std::cout << "BRIGHT_DATA Enviromental Variable could not be found\n";
    BRIGHT_DATA = (std::string) tmpBRIGHT_DATA;
#else
    BRIGHT_DATA = getenv("BRIGHT_DATA");
#endif
    return;
};

//String Transformations
std::string bright::to_str (int t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

std::string bright::to_str (double t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

std::string bright::to_str (bool t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

int bright::to_int (std::string s)
{
    return atoi( s.c_str() );
}

double bright::to_dbl (std::string s)
{
    return strtod( s.c_str(), NULL );
}

std::string bright::ToUpper(std::string strToConvert)
{
    //change each element of the string to upper case.
    for(unsigned int i = 0; i < strToConvert.length(); i++)
    {
        strToConvert[i] = toupper(strToConvert[i]);
    }
    return strToConvert;
}

std::string bright::ToLower(std::string strToConvert)
{
    //change each element of the string to lower case
    for(unsigned int i=0;i<strToConvert.length();i++)
    {
        strToConvert[i] = tolower(strToConvert[i]);
    }
    return strToConvert;
}

std::string bright::getFlag(char line[], int max_l)
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

std::string bright::Strip(std::string strToConvert, std::string strStrip)
{
    //change each element of the string to lower case
    int n_found = strToConvert.find(strStrip);
    while ( 0 <= n_found )
    {
        strToConvert.erase( n_found , strStrip.length() );
        n_found = strToConvert.find(strStrip);
    }
    return strToConvert;
}

std::string bright::MultiStrip(std::string strToConvert, std::string strMultiStrip)
{
    //change each element of the string to lower case
    for (int i = 0; i < strMultiStrip.length(); i++ )
    {
        strToConvert = Strip(strToConvert, strMultiStrip.substr(i, 1) );
    }
    return strToConvert;
}

std::string bright::StrictReplaceAll(std::string strToAlter, std::string subToRemove, std::string subToPlace)
{
    int n_found = strToAlter.find(subToRemove);
    while ( 0 <= n_found )
    {
        strToAlter.replace( n_found , subToRemove.length(), subToPlace );
        n_found = strToAlter.find(subToRemove);
    }
    return strToAlter;
};

std::string bright::LastChar(std::string s)
{
    //Returns the last character in a string.
    return s.substr(s.length()-1, 1);
}

std::string bright::SubFromEnd(std::string s, int n, int l)
{
    //Returns the splice of a string using negative indices.
    return s.substr(s.length()+n, l);
}

bool bright::ChainGreaterCompare(int a, int b, int c)
{
    //returns true id a <= b <= c and flase otherwise.
    return (a <= b && b <= c); 
}

bool bright::SubInString(std::string s, std::string sub)
{
    //Returns a boolean based on if the sub is in s.
    int n = s.find(sub);
    return ( 0 <= n && n < s.length() );
}

std::string bright::natural_naming(std::string name)
{
    std::string nat_name (name);

    //Replace Whitespace characters with underscores
    nat_name = bright::StrictReplaceAll(nat_name, " ",  "_");
    nat_name = bright::StrictReplaceAll(nat_name, "\t", "_");
    nat_name = bright::StrictReplaceAll(nat_name, "\n", "_");

    //Remove non-word characters
    int n = 0;
    while ( n < nat_name.length() )
    {
        if ( bright::words.find(nat_name[n]) == std::string::npos )
            nat_name.erase(n, 1);
        else
            n++;
    }

    //Make sure that the name in non-empty before continuing
    if (nat_name.length() == 0)
        return nat_name;

    // Make sure that the name doesn't begin with a number.
    if ( bright::digits.find(nat_name[0]) != std::string::npos)
        nat_name.insert(0, "_"); 

    return nat_name;
};



/*
 *  Vectorized functions
 */

std::vector<double> bright::normalized_delta(double x, std::vector<double> vec)
{
    // This functions find the normalized 
    // value of (x - vec[i]) for all elements i 
    // in the vector.
    //
    // This is equivelent to the fraction:
    //       (x - vec[i])
    //    -------------------
    //    (vec_max - vec_min) 

    std::vector<double> nd (vec.size(), 0.0);

    // Get the min and max out of the vector
    double vec_min = *std::min_element(vec.begin(), vec.end());
    double vec_max = *std::max_element(vec.begin(), vec.end());

    double min_max_delta = vec_max - vec_min;

    // Calculate the normalized delta for 
    // all i elements.
    for(int i = 0; i < vec.size(); i++)
    {
        nd[i] = (x - vec[i]) / min_max_delta;
    };

    return nd;
};






/* 
 * Array Helpers
 */

int bright::find_index_char(char * val, char ** arr, int arr_len)
{
    //Finds an element 'val' in array 'arr'
    //returns the index of val's first location
    //returns -1 if not found.
    //For Arrays of char strings

    if (arr_len < 0)
        arr_len = length_array(arr);

    for (int n = 0; n < arr_len; n++)
    {
        if (strcmp(arr[n], val) == 0)
            return n;
    };

    return -1;
};


//Math Helpers
const double bright::pi = 3.14159265359;
const double bright::N_A = 6.0221415 * pow(10.0, 23);
const double bright::bpcm2 = pow(10.0, -24); 

double bright::slope (double x2, double y2, double x1, double y1)
{
    //Finds the slope of a line.
    return (y2 - y1) / (x2 - x1);
};

double bright::SolveLine (double x, double x2, double y2, double x1, double y1)
{
    return (slope(x2,y2,x1,y1) * (x - x2)) + y2;
};

double bright::TANH(double x)
{
    return tanh(x);
};

double bright::COTH(double x)
{
    return 1.0 / tanh(x);
};


// File Helpers

bool bright::FileExists(std::string strFilename) 
{
    // Thank you intarwebz for this function!
    // Sepcifically: http://www.techbytes.ca/techbyte103.html
    struct stat stFileInfo; 
    bool blnReturn; 
    int intStat; 

    // Attempt to get the file attributes 
    intStat = stat(strFilename.c_str(), &stFileInfo); 

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
