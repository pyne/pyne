#!/usr/bin/env python

def get_isotopic_abundances():
    """Creates a dictionary of isotopic abundances based off of the IUPAC
    Technical Report entitled "Isotopic Compositions of the Elements 1997".

    Returns
    -------
    abundance : dict
        A dictionary where each key is an integer equal to the ZAID identifier
        for the isotope and the value is the atom fraction for the specified
        isotope relative to all atoms of that element.
    """

    abundance_file = open('abundances.txt', 'r')
    
    # Skip first line describing file
    abundance_file.readline()

    # Create dictionary
    abundance = {}
    abundance_by_Z = {i: [] for i in range(1,93)}

    # Read data
    for line in abundance_file:
        words = line.split()
        assert len(words) == 4

        # Read atomic number and mass number
        Z = int(words[0])
        A = int(words[2])
        zaid = Z*1000 + A
        
        # Read value and add to dictionary
        val = float(words[3])
        abundance[zaid] = val
        abundance_by_Z[Z].append((zaid,val))


    # Check data
    approx_equal = lambda a, b: abs(a - b) < 1e-8
    for Z in abundance_by_Z:
        total = 0.0

        # Skip elements with no stable isotopes
        if not abundance_by_Z[Z]:
            continue
    
        # Add abundance
        for zaid, val in abundance_by_Z[Z]:
            total += val

        # Check if sum is 100.0 (within floating point error)
        try:
            assert approx_equal(total, 100.0)
        except AssertionError:
            print("Atomic number: {0}".format(Z))
            for zaid in abundance_by_Z[Z]:
                print("  {0} = {1}".format(zaid, abundance[zaid]))
            print("  Total = {0}".format(total))
            raise
    
    # Return dictionary
    return abundance, abundance_by_Z
