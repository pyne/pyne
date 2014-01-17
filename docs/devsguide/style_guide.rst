.. _devsguide_styleguide:

===========
Style Guide
===========
PyNE is a polyglot project about a technical subject with many indpendent developers
and users. To keep our heads on straight we have adopted the following styles and 
conventions.  We use these throughout the code base to ensure consistency. 

----------------------------------
Rules to Write By
----------------------------------
It is important to refer to things and concpets by their most specific name.
When writing PyNE code or documentation please use technical terms approriately.
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
* The capilized project name is "PyNE" while the lowercase project name is "pyne".
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

***************
Canonical Forms
***************
* We have canonical forms; use them! For example, if a function accepts a nuclide 
  as an argument then you should use nucname to ensure that it accepts all possible 
  nuclide name spellings. If a function accepts a rection name then use the rxname
  module. Turn materials specifications into Material objects.  And so on...
* Mesh is the gentle giant of canonical forms. Use its strong, kind arms when dealing
  with geometries.
* If a canonical form doesn't exist, follow these steps:

    1. invent one, and
    2. make it awesome.

* Making a canonical form great may take time and many iterations. Do not give up.
* Interim working solutions are better than the best solution never.
* HDF5 is the preffered persistance format for structured data.

************
Expectations
************
* Code must have associated tests and adequate documentation.  
* The *only* exceptions to not having tests and documentation are when merging in and
  slowly integrating legacy code or code not originally originally written for pyne.
* Without both tests and documentation, the code must be marked as experimental.
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
one space. Only ``cdef``, ``cpdef``, or ``ctypedef`` one variable per line. 
Grouping ``cdef`` statements is allowed.  For example,

.. code-block:: cython

    # Good
    cdef int n
    cdef char* s
    cdef Material mat = Material()
    cdef int true_enough(x):
        return 1

    # Bad
    cdef  char *s
    cpdef int i, j, k
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
Pointers and refernces may be either zero or one space away from the type name.
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
be no blank lines between the propert declaration line and the following line.

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


.. _PEP8: http://www.python.org/dev/peps/pep-0008/
