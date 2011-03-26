// Header for general library file.

#if !defined(_bright_)
#define _bright_

//standard libraries
#include <string>
#include <string.h>
#include <sstream>
#include <cctype>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <exception>
#include <sys/stat.h> 
#include <vector>
#include <algorithm>

/*** Macros ***/
#define length_array(a) ( sizeof ( a ) / sizeof ( *a ) )

#ifdef _WIN32
    #define isnan(x) ((x) != (x))
#endif

namespace bright {

//Bright Globals
void bright_start ();

extern std::string BRIGHT_DATA;

//String Transformations
static std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static std::string digits = "0123456789";
static std::string words = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_";


std::string to_str (int);
std::string to_str (double);
std::string to_str (bool);

int to_int (std::string);

double to_dbl (std::string);

std::string ToUpper(std::string);

std::string ToLower(std::string);

std::string getFlag(char [], int);

std::string Strip(std::string, std::string);

std::string MultiStrip(std::string, std::string);

std::string StrictReplaceAll(std::string, std::string, std::string);

std::string LastChar(std::string);

std::string SubFromEnd(std::string, int = -1, int = 1);

bool ChainGreaterCompare(int, int, int);

bool SubInString(std::string, std::string);

std::string natural_naming(std::string);


// Vectorized Functions
std::vector<double> delta_vector(double, std::vector<double>);
std::vector<double> normalized_delta(double, std::vector<double>);

bool sorted_index_comparator(std::pair<int, double>, std::pair<int, double>);
std::vector<int> sorted_index(std::vector<double>);


//Array Methods
template <class T>
int find_index(T val, T * arr, int arr_len = -1)
{
    //Finds an element 'val' in array 'arr'
    //returns the index of val's first location
    //returns -1 if not found.

    if (arr_len < 0)
        arr_len = sizeof(arr) / sizeof(T);

    for (int n = 0; n < arr_len; n++)
    {
        if (val == arr[n])
            return n;
    };

    return -1;
};

int find_index_char(char *, char **, int = -1);

//Math Helpers
extern const double pi;
extern const double N_A;	//Avagardo's Number
extern const double bpcm2; 	//barns per cm^2

double slope (double, double, double, double);
double SolveLine (double, double, double, double, double);

double TANH(double);
double COTH(double);


// File Helpers
bool FileExists(std::string); 



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




class BadFuelForm : public std::exception
{
//Exception for valid fuel form.
public:
    BadFuelForm () {};
    ~BadFuelForm () throw () {};

    static char * name ()
    {
        return (char *) "BadFuelForm";
    };

    virtual const char* what() const throw()
    {
        std::string BFFstr ("FUEL COMPOSITION NOT COMPUTABLE!");
        return (const char *) BFFstr.c_str();
    };
};




class BisectionMethodNotPerformed : public std::exception
{
//Exception for when the bisection method is not calculated.
public:
    BisectionMethodNotPerformed ()
    {
        errstr = "Bisection method was not performed.";
    };
    BisectionMethodNotPerformed (std::string calctype)
    {
        errstr = "Bisection method durring " + calctype + " calculation was not performed.";
    };
    ~BisectionMethodNotPerformed () throw () {};

    static char * name ()
    {
        return (char *) "BisectionMethodNotPerformed";
    };

    virtual const char* what() const throw()
    {
        return (const char *) errstr.c_str();
    };
private:
    std::string errstr;
};


// End bright namespace
};

#endif
