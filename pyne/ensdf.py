import re

from pyne import nucname



def half_life(ensdf):
    """Grabs the half-lives from an ENSDF file.

    Parameters
    ----------
    ensdf : str or file-like object
        ENSDF file to inspect for half-life data

    Returns
    -------
    hl : list of 4-tuples
        List of tuples where the indices match

        1. from_nuclide, int (zzaaam)
        2. to_nuclide,  int (zzaaam)
        3. half_life, float (seconds)
        4. branch_ratio, float (frac)

    """
    opened_here = False
    if isinstance(path, bsestring):
        ensdf = open(ensdf, 'r')
        opened_here = True

    hl = []

    # Run through the file
    lines = ensdf.readlines()
    for i, line in enumerate(lines):
        pass

    if opened_here:
        ensdf.close()
