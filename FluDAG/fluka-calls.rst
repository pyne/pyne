
Fluka Calls
-----------
The collection of functions that can be substituted from external libraries is documented by FLUKA and 
summarized here.

==============    ========================  ========  ===========  ========   ==========
Function          Comment                   Language  Early Stage   Working    Complete
==============    ========================  ========  ===========  ========   ==========
flukam            Fluka main                Fortran                              yes
flgfwr            sets FLUKA geometry flag  C++                                  yes
g1wr              main tracking routine     C++                       yes
g1rtwr            dummy routine                                                  yes
lkdbwr, lkfxwr    dummy routines                                                 yes
conhwr            not implemented           Fortran                              yes
inihwr            History initialization                                         yes
jomiwr            Geometry initialization   C++                                  yes
lkwr              Particle localization     C++                                  yes
fldwr             Return B-field value      C++          yes                        
nrmlwr            Normal vector at surface  C++                                  yes
rg2nwr            Region name to number     C++                                  yes
rgrpwr            Scoring with history                   yes         
==============    ========================  ========  ===========  ========   ==========

FLUKA does not make its source code available; testing is underway to ensure the correct 
implementation of these calls.


