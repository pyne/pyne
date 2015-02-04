#. **id (zas)**: This type places the charge of the nucleus out front, then has three
   digits for the atomic mass number, and ends with four state digits (0 = ground,
   1 = first excited state, 2 = second excited state, etc).  Uranium-235 here would
   be expressed as '922350000'.  This is th canonical form for nuclides.
#. **name**: This is the more common, human readable notation.  The chemical symbol
   (one or two characters long) is first, followed by the atomic weight.  Lastly if
   the nuclide is metastable, the letter *M* is concatenated to the end.  For
   example, 'H-1' and 'Am242M' are both valid.  Note that nucname will always
   return name form with the dash removed and all letters uppercase.
#. **zzaaam**: This type places the charge of the nucleus out front, then has three
   digits for the atomic mass number, and ends with a metastable flag (0 = ground,
   1 = first excited state, 2 = second excited state, etc).  Uranium-235 here would
   be expressed as '922350'.
#. **zzzaaa**: This type places the charge of the nucleus out front, then has three
   digits for the atomic mass.  It contains no information about the excited state.
#. **zzllaaam**: The ZZLLAAAM naming convention is similar to name form.  However, it 
   is preceded by the nuclides two AA numbers, followed by the two LL characters.  
   Of the two LL characters, only the first letter in the chemical symbol is uppercase, 
   the dash is always present, and the the meta-stable flag is lowercase.  For 
   instance, '95-Am-242m' is the valid serpent notation for this nuclide.
#. **SZA**: This type places three state digits out front, the charge of the nucleus in 
   the middle, and then has three digits for the atomic mass number. Uranium-235M here 
   would be expressed as '1092235'.  
#. **MCNP**: The MCNP format for entering nuclides is unfortunately
   non-standard.  In most ways it is similar to zzaaam form, except that it
   lacks the metastable flag.  For information on how metastable isotopes are
   named, please consult the MCNPX documentation for more information.
#. **Serpent**: The serpent naming convention is similar to name form.
   However, only the first letter in the chemical symbol is uppercase, the
   dash is always present, and the the meta-stable flag is lowercase.  For
   instance, 'Am-242m' is the valid serpent notation for this nuclide.
#. **NIST**: The NIST naming convention is also similar to the Serpent form.
   However, this convention contains no metastable information.  Moreover, the
   A-number comes before the element symbol.  For example, '242Am' is the
   valid NIST notation.
#. **CINDER**: The CINDER format is similar to zzaaam form except that the
   placement of the Z- and A-numbers are swapped. Therefore, this format is
   effectively aaazzzm.  For example, '2420951' is the valid cinder notation
   for 'AM242M'.
#. **ALARA**: In ALARA format, elements are denoted by the lower case atomic symbol. 
   Nuclides are specified by appending a semicolon and A-number. For example, "fe" 
   and "fe:56" represent elemental iron and iron-56 respectively. No metastable 
   flag exists.
#. **Groundstate**:  In Groundstate format, the nuclide is stored in a form similar 
   to the standard id form, but the last four digits are zero to eliminate the 
   information about the nuclide's state.  
#. **state_id**: The state id naming convention uses the form zzzaaassss. It is 
   different from the canonical zzzaaassss form in that ssss refers to a list 
   of states by ordered by energy. This is derived from the levels listed in the 
   ENSDF files for a given nuclide. Using this form is dangerous as it may change 
   with new releases of ENSDF data.
