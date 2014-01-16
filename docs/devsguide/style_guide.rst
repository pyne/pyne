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

* The terms 'element' or 'elemental' refer specifically to chemical elements,
  oxygen, uranium, etc.  Use 'element' only when the physics relies on manipulating 
  elements and not their underlying isotopes (ie aqueous reprocessing).
* The term 'isotope' refers only to a collection of nuclides containing the 
  same Z-number.  **Do not confuse this with 'nuclide'**!
* The term 'isotone' refers only to a collection of nuclides containing the 
  same neutron number.
* The term 'isobar' refers only to a collection of nuclides containing the 
  same A-number.
* The term 'isomer' refers only to a collection of nuclides containing the 
  same A-number and the same Z-number but possessing differing internal energy 
  states.
* The term 'nuclide' may refer to any species with a nucleus. This is the most
  general term and encompasses isotopes, isotones, isobars, and isomers.
* Always provide units! 
* Use square-brace notation to mark units in text, ie [sec].
* User-facing APIs should be as generic and robust as possible.  
* We have canonical forms; use them! For example, if a function accepts a nuclide 
  as an argument then you should use nucname to ensure that it accepts all possible 
  nuclide name spellings. If a function accepts a rection name then use the rxname
  module. Turn materials specifications into Material objects.  And so on...
* Mesh is the gentle giant of canonical forms. Use its strong, kind arms when dealing
  with geometries.
* If a canonical form doesn't exist, follow these steps:

    1. invent one, and
    2. make it awesome.

* Making a canonical form great may take time and many iterations.  Don't give up.
* Interim working solutions are better than the best solution never.
* Nuclear data belongs in nuc_data.h5.
* Views belong in the `gui` sub-package.
* Tests belong in the top-level `tests` directory.
* Documentation belongs in the top-level `docs` directory.
* Write code in whatever language you want but make sure that it is exposed to Python.
  Python is the glue that holds the universe together. It is the NULL and the INT_MAX.
* Code must have associated tests and adequate documentation.  
* The *only* exceptions to not having tests and documentation are when merging in and
  slowly integrating legacy code or code not originally originally written for pyne.
* Without both tests and documentation, the code must be marked as experimental.
* The capilized project name is "PyNE" while the lowercase project name is "pyne".
* Nothing says "I <3 PyNE" quite like an ASCII art dragon.


*************************
Variable Name Conventions
*************************
Please use the following patterns when writing pyne code. These rules should 
only be ignored if they would cause name conflicts. These variable names are 
considered to have the same semantic meaning throughout the entire code base.

* `xs` stands for "cross section".
* `rx` stands for "reaction".
* `ve` stands for "volume element".
* `nuc` stands for a "nuclide" in id form.
* `nuc_name` stands for a "nuclide" in a string form.
* `iso` stands for an "isotope" in id form.
* `iso_name` stands for an "isotope" in string form.
* `mat` stands for "material".

-------------------
Python Style Guide 
-------------------
PyNE uses PEP8 for all Python code.  The following rules apply where PEP8 
is open to interpretation.

* Use absolute imports (`import pyne.material`) rather than explicit relative imports
  (`import .material`). Implicit relative imports (`import material`) are never
  allowed.
* Use 'single quotes' for string literals, and """triple double quotes""" for 
  docstrings. Double quotes are allowed to prevent single quote escaping, 
  e.g. "Y'all c'mon o'er here!"
