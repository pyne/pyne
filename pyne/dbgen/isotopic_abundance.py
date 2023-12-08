#!/usr/bin/env python
from __future__ import print_function
import pkgutil
from pyne.utils import QA_warn

QA_warn(__name__)


def get_isotopic_abundances():
    """Creates a dictionary of isotopic abundances based off of [1].

    [1] M. Berglund, M. Wieser: Isotopic compositions of the elements 2009
        (IUPAC Technical Report).
        Pure Appl. Chem., 2011, Vol. 83, No. 2, pp. 397-410
        http://dx.doi.org/10.1351/PAC-REP-10-06-02

    Returns
    -------
    abundance : dict
        A dictionary where each key is an integer equal to the ZAID identifier
        for the isotope and the value is the atom fraction for the specified
        isotope relative to all atoms of that element.
    """

    abundance_file = (
        pkgutil.get_data("pyne.dbgen", "abundances.txt").decode().split("\n")
    )

    # Create dictionary
    abundance = {}
    abundance_by_Z = dict([(i, []) for i in range(1, 93)])

    # Read data
    for line in abundance_file:
        if line.startswith("#") or not line:
            continue
        words = line.split()
        assert len(words) == 4

        # Read atomic number and mass number
        Z = int(words[0])
        A = int(words[2])
        nuc = (Z * 1000 + A) * 10000
        name = "-".join(words[1:3])

        # Read value and add to dictionary
        val = 0.01 * float(words[3])
        abundance[nuc] = val
        abundance_by_Z[Z].append((name, val))

    # Check data
    approx_equal = lambda a, b: abs(a - b) < 1e-8
    for Z in abundance_by_Z:
        total = 0.0

        # Skip elements with no stable isotopes
        if not abundance_by_Z[Z]:
            continue

        # Add abundance
        for name, val in abundance_by_Z[Z]:
            total += val

        # Check if sum is 1.0 (within floating point error)
        try:
            assert approx_equal(total, 1.0)
        except AssertionError:
            print("Atomic number: {0}".format(Z))
            for name, val in abundance_by_Z[Z]:
                print("  {0} = {1}".format(name, val))
            print("  Total = {0}".format(total))
            raise

    # Return dictionary
    return abundance
