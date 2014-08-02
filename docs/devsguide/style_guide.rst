.. _devsguide_styleguide:

===========
Style Guide
===========
PyNE is a polyglot project about a technical subject with many independent developers
and users. To keep our heads on straight we have adopted the following styles and 
conventions.  We use these throughout the code base to ensure consistency. 

----------------------------------
Rules to Write By
----------------------------------
It is important to refer to things and concepts by their most specific name.
When writing PyNE code or documentation please use technical terms appropriately.
The following rules help provide needed clarity.

***********
Terminology
***********
* The terms *element* or *elemental* refer specifically to chemical elements,
  oxygen, uranium, etc.  Use *element* only when the physics relies on manipulating 
  elements and not their underlying isotopes (ie aqueous reprocessing).
* The term *isotope* refers only to a collection of nuclides containing the 
  same Z-number.  **Do not confuse this with *nuclide*!**
* The term *isotone* refers only to a collection of nuclides containing the 
  same neutron number.
* The term *isobar* refers only to a collection of nuclides containing the 
  same A-number.
* The term *isomer* refers only to a collection of nuclides containing the 
  same A-number and the same Z-number but possessing differing internal energy 
  states.
* The term *nuclide* may refer to any species with a nucleus. This is the most
  general term and encompasses isotopes, isotones, isobars, and isomers.

**********
Interfaces
**********
* Always provide units! 
* Use square-brace notation to mark units in text, ie [sec].
* User-facing APIs should be as generic and robust as possible.  
* Nuclear data belongs in nuc_data.h5.
* Views belong in the ``gui`` sub-package.
* Tests belong in the top-level ``tests`` directory.
* Documentation belongs in the top-level ``docs`` directory.
* The capitalized project name is "PyNE" while the lowercase project name is "pyne".
* Write code in whatever language you want but make sure that it is exposed to Python.
  Python is the glue that holds the universe together. It is the ``NULL`` and the 
  ``INT_MAX``.

*************************
Variable Name Conventions
*************************
Please use the following patterns when writing pyne code. These rules should 
only be ignored if they would cause name conflicts. These variable names are 
considered to have the same semantic meaning throughout the entire code base.

* ``xs`` stands for "cross section".
* ``rx`` stands for "reaction".
* ``ve`` stands for "volume element".
* ``nuc`` stands for a "nuclide" in id form.
* ``nuc_name`` stands for a "nuclide" in a string form.
* ``iso`` stands for an "isotope" in id form.
* ``iso_name`` stands for an "isotope" in string form.
* ``mat`` stands for "material".
* ``m`` and ``mesh`` are used for Mesh instances.
* ``idx`` stands for "indices". When used with meshes it is the canonical integer 
  identifier for an entity. This is not necessarily iteration order.
* ``i`` is the iteration variable of mesh indices, ie ``for i in idx``.
* Generator names should start with the prefix ``iter``.

***************
Canonical Forms
***************
* We have canonical forms; use them! For example, if a function accepts a nuclide 
  as an argument then you should use nucname to ensure that it accepts all possible 
  nuclide name spellings. If a function accepts a reaction name then use the rxname
  module. Turn materials specifications into Material objects.  And so on...
* Mesh is the gentle giant of canonical forms. Use its strong, kind arms when dealing
  with geometries.
* If a canonical form doesn't exist, follow these steps:

    1. invent one, and
    2. make it awesome.

* Making a canonical form great may take time and many iterations. Do not give up.
* Interim working solutions are better than the best solution never.
* HDF5 is the preferred persistence format for structured data.

************
Expectations
************
* Code must have associated tests and adequate documentation.  
* The *only* exceptions to not having tests and documentation are when merging in and
  slowly integrating legacy code or code not originally originally written for pyne.
* Without both tests and documentation, the code must be marked as experimental.
  It should issue a ``pyne.utils.VnVWarning``.
* Have *extreme* empathy for your users.
* Be selfish. Since you will be writing tests you will be your first user.
* Nothing says "I <3 PyNE" quite like an ASCII art dragon.

-------------------
Python Style Guide 
-------------------
PyNE uses `PEP8`_ for all Python code.  The following rules apply where `PEP8`_
is open to interpretation.

* Use absolute imports (``import pyne.material``) rather than explicit relative imports
  (``import .material``). Implicit relative imports (``import material``) are never
  allowed.
* Use 'single quotes' for string literals, and """triple double quotes""" for 
  docstrings. Double quotes are allowed to prevent single quote escaping, 
  e.g. "Y'all c'mon o'er here!"
* We use sphinx with the numpydoc extension to autogenerate API documentation. Follow 
  the numpydoc standard for docstrings `described here <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.
* Simple functions should have simple docstrings.
* Lines should be at most 80 characters long. The 72 and 79 character recommendations
  from PEP8 are not required here.
* All Python code should be compliant with Python 2.7 and Python 3.3+.  At some 
  unforeseen date in the future, Python 2.7 support *may* be dropped.
* Tests should be written with nose using a procedural style. Do not use unittest
  directly or write tests in an object-oriented style.
* Test generators make more dots and the dots must flow!

-------------------
Cython Style Guide 
-------------------
Cython as a super-set language of Python should follow `PEP8`_ for all syntax 
that the two languages share.  Cython-specific syntax should follow these additional
rules.

***************************
cdefs, cpdefs, & ctypedefs
***************************
Separate ``cdef``, ``cpdef``, and ``ctypedef`` statements from the following type by 
exactly one space. In turn, separate the type from the variable name by exactly 
one space. Only ``ctypedef`` one variable per line. You may ``cdef`` or ``cpdef``
multiple variables per line as long as these are simple declarations - multiple 
assignment, references, or pointers are not allowed. Grouping ``cdef`` statements 
is allowed.  For example,

.. code-block:: cython

    # Good
    cdef int n
    cdef char* s
    cpdef int i, j, k
    cdef Material mat = Material()
    cdef int true_enough(x):
        return 1

    # Bad
    cdef  char *s
    cdef char * s, * t, * u, * v
    cdef double x=42, y=x+1, z=x*y 
    cdef Material     mat   = Material()
    cdef   int   falsified(x):
        return 0

Inside of a function, place all ``cdef`` statements at the top of the function body.

.. code-block:: cython

    # Good
    cdef int true_enough(x):
        cdef int i = x
        cdef int rtn
        rtn = i + 42
        return rtn 

    # Bad
    cdef int falsified(x):
        cdef int i = x, j = -42
        j += i
        cdef int rtn = j / j - 1
        return rtn 

****************************
cimport & include statements
****************************
The ``cimports`` should follow the same rules defined in `PEP8`_ for 
``import`` statements.  If a module is both imported and cimported, the 
cimport should come before the import.

Do not use ``include`` statements.

*******************
Error return values
*******************
When declaring an error return value with the ``except`` keyword, use one 
space on both sides of the ``except``. If in a function definition, there should 
be no spaces between the error return value and the colon ``:``.  Avoid ``except *``
unless it is needed for functions returning ``void``. 

.. code-block:: cython

    # Good
    cdef void redwood() except *
    cdef int sequoia(x) except +:
        ...

    # Bad
    cdef char * spruce(x) except *:
    cdef int fir(x)    except   +  :
        ...


*********************
Pointers & References
*********************
Pointers and references may be either zero or one space away from the type name.
If followed by a variable name, they must be one space away from the variable name.
Do not put any spaces between the reference operator ``&`` and the variable name.

.. code-block:: cython

    # Good
    cdef int& i
    cdef char * s
    i = &j

    # Bad
    cdef int &i
    cdef char *s
    i = & j


*******
Casting
*******
When casting a variable there must be no whitespace between the opening ``<`` and
the type.  There must one space between the closing ``>`` and the variable.

.. code-block:: cython

    # Good
    <float> i
    <void *> s

    # Bad
    < float >i
    <void*>  s

*****
Loops
*****
Use Python loop syntax - ``for i in range(10):``.  Other for-loop constructs are 
deprecated and must be avoided.

****************
Property Keyword
****************
Properties are great! There should be exactly one space between the ``property``
keyword and the attribute name.  There may be no spaces between the attribute 
name and the colon ``:``.  All properties should have docstrings. There should 
be no blank lines between the property declaration line and the following line.

.. code-block:: cython

    # Good
    property has_cone:
        """This class has a cone.
        """
        def __get__(self):
            ...

    # Bad
    property    has_cone :

        def __get__(self):
            ...

**************************************************
Type Declarations, Extern, Public, API, & Readonly
**************************************************
Type declarations, the ``extern`` keyword, the ``public`` keyword, the ``api`` 
keyword, and the ``readonly`` keyword should always be followed by a single space.

.. code-block:: cython

    # Good
    cdef extern void * v
    cdef public api int i
    def sequoia(int x):
        ...

    # Bad
    cdef extern         void * v
    cdef public  api    int    i
    def spruce(int   x):
        ...

**************
Documentation
**************
In addition to following the numpydoc convention, also include the function or method 
signature as the first line in the docstring.  This helps sphinx print out the 
signature.  Include type information in this signature if available and relevant.

-------------------
C/C++ Style Guide 
-------------------
As software that is meant to be exposed to Python, C/C++ code written for pyne
has special needs.  Existing single-language style guides are non-idiomatic across 
the language barrier.  This style guide attempts to rectify this impedance 
mismatch by defining a hierarchy of style guides and special rules to follow that
make C/C++ more PyNEthonic. Legacy codes not originally written for pyne in these 
languages need not be migrated to this style.  While a custom style may not be 
ideal in terms of leveraging linters and style checker tools, the benefits 
in readability and portability outweigh this cost.  

The aim is to have all languages be as similar and have as idiomatic of APIs for that 
language as possible. 

Except as noted below, C/C++ code should adhere to the rules laid out in the 
following style guides in order of preference:

1. `PEP8`_
2. `The Linux Kernel Coding Style <http://www.maultech.com/chrislott/resources/cstyle/LinuxKernelCodingStyle.txt>`_
3. `The Google C++ Style Guide <http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml>`_

This section was forked from the `ROS C++ Style Guide <http://wiki.ros.org/CppStyleGuide>`_.
If you require clarification on a particular syntax or idiom, please ask!

*****
Files
*****
Files may have under_scores.

C source files have the extension ``.c``.

C++ source files have the extension ``.cpp``.

Header files have the extension ``.h``.

If the file primarily implements a class, name the file after the class.

****************************
Classes, Typedefs, & Structs
****************************
Class names are CapCased:

.. code-block:: c++

    class ExampleClass;

**Exception:** if the class name contains a short acronym, the acronym itself 
should be all capitals:

.. code-block:: c++

    class HokuyoURGLaser;

Name the class after what it is. If you can't think of what it is, perhaps you 
have not thought through the design well enough.

Class names should be nouns. 

Typedef names should be lowercase_with_underscores, like primitive C/C++ and 
Python types.

Struct names should be CapCased if they have non-trivial member functions
and are more class-like.  

However, if a struct is meant to be used primarily as compound data type 
it should have a lowercase_with_underscores name, like typedefs.

*********
Functions
*********
Functions and their arguments are lowercase_with_underscores:

.. code-block:: c++

    int example_func(int example_arg);

Functions usually performs an action, so the name should make clear what it does.
Function names thus should be verbs.

*********
Variables
*********
Variable names are lowercase_with_underscores.

Integral iterator variables can be very short, such as i, j, k. Be consistent in 
how you use iterators (e.g., i on the outer loop, j on the next inner loop).

STL iterator variables should indicate what they are iterating over:

.. code-block:: c++

    std::list<int> pid_list;
    std::list<int>::iterator pid_it;

*********
Constants
*********
Constants, wherever they are used, are ALL_CAPITALS.

****************
Member Variables
****************
Variables that are members of a class are lowercase_with_underscores.
Private and protected member variables start with a single leading underscore.
Public member variables do not have a leading underscore.

.. code-block:: c++

    int public_x;
    int _protected_y;
    int _private_z;

****************
Global Variables
****************
Global variables should never be used. 

**Exception:** a file may contain a main() function. 

**********
Namespaces
**********
Namespace names, like Python module names, are lowercase *without* underscores.

Everything should be in a namespace.  Anonymous namespaces are encouraged to help
meet this requirement.

The bodies of namespace declaration and definition are not indented. This is 
the same as the `GCSG`_.

Never use a ``using namespace`` directive. Using-declarations inside of class 
or function scope, which only grab the names you intend to use, are allowed.

.. code-block:: c++

    // Good
    using std::list;    // I want to refer to std::list as list
    using std::vector;  // I want to refer to std::vector as vector

    // Bad, because it imports all names from std::
    using namespace std;  

***************
Access Patterns
***************
We are all adults here. Everything should be public.  Use private and protected 
variables only when absolutely necessary.

*************************
Accessors/Mutator Pattern
*************************
Avoid getter and setter member functions. This pattern increases code volume, 
inlining is not guaranteed, and slows down run times.

Use this pattern only if implementing a Python/Cython-like property where
getting or setting a member variable is non-trivial. In these cases, the 
storage variable should be named with a leading underscore (even though it may be 
public) and the get/set names should have the same name as the variable but without
the leading underscore:

.. code-block:: c++

    class WithAnX {
     public:
      // storage variable
      int _x;

      // getter
      int x();

      // setter
      void x(int value);
    }


**********
Formatting
**********
Indent each block by 2 spaces. Never insert literal tab characters.

The contents of a namespace are not indented.

We are all friends here! Braces should be `cuddled <http://gskinner.com/blog/archives/2008/11/curly_braces_to.html>`_:

.. code-block:: c++

    if (a < b) {
      ...
    } else {
      ...
    }

Braces may be omitted if the enclosed block is a single-line statement:

.. code-block:: c++

    if (a < b)
      x = 2*a;

Only single line comments should be used.  Multi-line comments are inconsistent
and not allowed.

.. code-block:: c++

    // This is OK

    /* This is not OK */

    /* What is even going on here?!
     * All I can see are the stars...
     */

***********
Line Length
***********
Maximum line length is 80 characters.

**************
Include Guards
**************
All headers must be protected against multiple inclusion by #ifndef guards.
These guards ought to be UUIDs:

.. code-block:: c++

    #ifndef PYNE_W7WGLJVRGRDH7G47RDHRLLCP2A
    #define PYNE_W7WGLJVRGRDH7G47RDHRLLCP2A
    ...
    #endif

Use this command for generating UUIDs:

.. code-block:: bash

    $ python -c "import uuid; import base64; print('PYNE_' + base64.b32encode(uuid.uuid4().bytes).decode().strip('='))"

This guard should begin before any other code and should end at the end of the file.



*************
Documentation
*************
All code must be documented. We use doxygen to auto-document our code. 
All functions, methods, classes, variables, enumerations, and constants 
should be documented.

***************
Console Output
***************
Avoid printf if in C++.  Use ``std::cout`` instead.

******
Macros
******
Avoid preprocessor macros whenever possible. Unlike in-line functions and const 
variables, macros are neither typed nor scoped.

***********
Inheritance
***********
When overriding a virtual method in a subclass always declare it to be virtual
so that the reader knows what's going on.

**********
Exceptions
**********
Built-in exceptions are the preferred error-reporting mechanism, 
as opposed to returning integer error codes or custom exception mechanisms.

Do not throw exceptions from destructors.

Do not throw exceptions from callbacks that you don't invoke directly.

**************
Calling exit()
**************
Only call ``exit()`` at a well-defined exit point for the application.

Never call ``exit()`` in a library.

***********
Portability
***********
Portability counts. 

Do not use uint as a type. Instead use unsigned int.

Call ``isnan()`` from within the std namespace, i.e.: ``std::isnan()``.

.. _PEP8: http://www.python.org/dev/peps/pep-0008/
.. _GCSG: http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml
