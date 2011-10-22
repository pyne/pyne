import numpy as np
import re

_if_idx_str = ("""if (exist("idx", "var"));\n"""
              """  idx = idx + 1;\n"""
              """else;\n"""
              """  idx = 1;\n"""
              """end;"""
              )

_num_pattern = "([0-9]+[.]?[0-9]*[Ee]?[+-]?[0-9]*)"

_numpy_array_pattern = r"\[[0-9\sEe+-.,%]*\]"
_matlab_array_pattern = r"\[[0-9\sEe+-.%]*\]"

_comment_array_pattern = r"\[([\d\sEe+-.,]*\s*[%#]+\s*[%#\w\s.,+-]*)\]"
_comment_line_pattern = r"([\d\s\tEe+-.]*)\s*([%#]*)\s*([#\w\s]*\n)"

_lhs_variable_pattern = r"(\w+)\s*(\(idx.*?\))"
_rhs_variable_pattern = r"(\w+)\s*\(idx.*?\)\s*=\s*(.*)"

_zeros_pattern = r"(zeros)\((.*)\)"

_detector_pattern = r"(DET\w+)\s*=\s*np.array\("

def _replace_comments(s):
    """Replaces matlab comments with python arrays in string s."""
    s = s.replace('%', '#')
    return s

def _replace_semicolons(s):
    """Replaces matlab semicolons with nothing in string s."""
    s = s.replace(';', '')
    return s

def _replace_arrays(s):
    """Replaces matlab arrays with numpy arrays in string s."""

    # Replace matlab arrays with python lists
    arrays = re.findall(_matlab_array_pattern, s)
    for a in arrays:
        new_a = re.sub(_num_pattern, lambda mo: mo.group(0) + ',', a)
        s = s.replace(a, new_a)

    # Encapsulate python lists in numpy arrays
    s = re.sub(_numpy_array_pattern, lambda mo: 'np.array(' + mo.group(0) + ')', s)

    return s


def parse_res(resfile, write_py=False):
    """Converts a serpent results *_res.m output file to a dictionary (and 
    optionally to a *_res.py file).

    Parameters
    ----------
    resfile : str or file-like object
        Path to results file or a res file handle.
    write_py : bool, optional
        Flag for whether to write the res file to an analogous python file.

    Returns
    -------
    res : dict
        Dictionary of the parsed results.  Please see the Serpent manual for
        a complete description of contents.

    """
    if isinstance(resfile, basestring):
        with open(resfile, 'r') as mfile:
            f = mfile.read()
    else:
        f = resfile.read()

    # Keep comments around
    f = _replace_comments(f)

    # Grab the number of 'if' statements
    IDX = f.count(_if_idx_str)

    # Replace if statements with something more meaningful
    fpart = f.partition(_if_idx_str)
    f = fpart[0] + "idx = 0" + fpart[2]
    f = f.replace(_if_idx_str, 'idx += 1')

    # Replace matlab Arrays
    f = _replace_arrays(f)

    # Add imports to header
    header = "import numpy as np\n\n"

    # Find all variables and shape
    vars_shape = np.unique( re.findall(_lhs_variable_pattern, f) )
    vars_dtype = dict( re.findall(_rhs_variable_pattern, f) )
    # Initialize variables to zero
    header = header + "# Initialize variables\n"
    for vs in vars_shape:
        # Determine shape
        s = re.search(r'\[.*:(.*?)\]', vs[1])
        if s == None:
            vs_shape = ""
        else:
            vs_shape = s.group(1)
            vs_shape = vs_shape.split()
            vs_shape = ", ".join(vs_shape)

        # Determine Data type
        rhs = vars_dtype[vs[0]]
        if ("\'" in rhs) or ("\"" in rhs):
            dt = "'S{0}'".format(int( s.group(1) ))
            vs_shape = ""
        elif ('.' in rhs) or ('E' in rhs) or ('e' in rhs):
            dt = "float"
        else:
            dt = "int"

        zero_line = "{0} = np.zeros([{1}, {2}], dtype={3})\n".format(vs[0], IDX, vs_shape, dt)
        header = header + zero_line

    # Add IDx to file
    if 0 < IDX:
        header = header + "\n\n# Maximum Index\n\nIDX = {0}\n\n".format(IDX)

    # Add header to file
    f = header + f

    # Replace variable overrides
    vars = np.unique( re.findall("(" + _lhs_variable_pattern + ")", f) )
    for v in vars:
        f = f.replace(v[0], "{0}[idx] ".format(v[1]))

    # Remove semicolons
    f = _replace_semicolons(f)

    # Write the file out
    if write_py:
        if isinstance(resfile, basestring):
            new_filename = resfile.rpartition('.')[0] + '.py'
        else:
            new_filename = resfile.name.rpartition('.')[0] + '.py'
        with open(new_filename, 'w') as pyfile:
            pyfile.write(f)

    # Execute the adjusted file
    res = {}
    exec(f, {}, res)

    return res



def parse_dep(depfile, write_py=False, make_mats=True):
    """Converts a serpent depletion *_dep.m output file to a dictionary (and 
    optionally to a *_dep.py file).

    Parameters
    ----------
    depfile : str or file-like object
        Path to depletion file or a dep file handle.
    write_py : bool, optional
        Flag for whether to write the dep file to an analogous python file.
    make_mats : bool, optional
        Flag for weather or not to build Materials out of mass data.

    Returns
    -------
    dep : dict
        Dictionary of the parsed depletion information.  Please see the Serpent 
        manual for a complete description of contents.

    """
    if isinstance(depfile, basestring):
        with open(depfile, 'r') as mfile:
            f = mfile.read()
    else:
        f = depfile.read()

    # Keep comments around
    f = _replace_comments(f)

    # Replace matlab Arrays
    f = _replace_arrays(f)

    # Now to find and convert arrays that have comments in them
    comment_arrays = re.findall("(" + _comment_array_pattern + ")", f)
    for ca in comment_arrays:
        new_ca = ca[0]
        comment_lines = re.findall("(" + _comment_line_pattern + ")", ca[1])
        for cl in comment_lines:
            new_cl = re.sub(_num_pattern, lambda mo: mo.group(0) + ',', cl[1])
            if new_cl[0] == '\n':
                new_cl = "\n    [" + new_cl.strip() + "], "
            else:
                new_cl = "    [" + new_cl.strip() + "], "

            new_ca = new_ca.replace(cl[1], new_cl)

        new_ca = 'np.array( ' + new_ca + ' )'    
        f = f.replace(ca[0], new_ca)

    # Indent close of array
    f = f.replace("\n] )", "\n    ] )")

    # Replace MatLab zeros with numpy zeros
    f = re.sub(_zeros_pattern, lambda mo: "np.zeros((" + mo.group(2) + "))", f)

    # Replace some math operators
    f = f.replace('.*', "*")
    f = f.replace('./', "/")

    # Remove semicolons
    f = _replace_semicolons(f)

    # Add imports to header
    header = "import numpy as np\n"
    if make_mats:
        header += "from pyne.material import Material\n"
    header += "\n"

    # Add materials
    footer = ""
    if make_mats:
        mat_gen_line = "{name}MATERIAL = [{name}VOLUME * Material(dict(zip(ZAI[:-2], {name}MDENS[:-2, col]))) for col in range(len(DAYS))]\n"
        footer += "\n\n# Construct materials\n"
        base_names = re.findall('(MAT_\w*_)MDENS = ', f)
        for base_name in base_names:
            footer += mat_gen_line.format(name=base_name)
        footer += "TOT_MATERIAL = [Material(dict(zip(ZAI[:-2], TOT_MASS[:-2, col]))) for col in range(len(DAYS))]\n"

    # Add header & footer to file
    f = header + f + footer

    # Write the file out
    if write_py:
        if isinstance(depfile, basestring):
            new_filename = depfile.rpartition('.')[0] + '.py'
        else:
            new_filename = depfile.name.rpartition('.')[0] + '.py'
        with open(new_filename, 'w') as pyfile:
            pyfile.write(f)

    # Execute the adjusted file
    dep = {}
    exec(f, {}, dep)

    return dep



def parse_det(detfile, write_py=False):
    """Converts a serpent detector *_det.m output file to a dictionary (and 
    optionally to a *_det.py file).

    Parameters
    ----------
    detfile : str or file-like object
        Path to detector file or a det file handle.
    write_py : bool, optional
        Flag for whether to write the det file to an analogous python file.

    Returns
    -------
    det : dict
        Dictionary of the parsed detector.  Please see the Serpent manual for
        a complete description of contents.

    """
    if isinstance(detfile, basestring):
        with open(detfile, 'r') as mfile:
            f = mfile.read()
    else:
        f = detfile.read()

    # Keep comments around
    f = _replace_comments(f)

    # Replace matlab Arrays
    f = _replace_arrays(f)

    # Find detector variable names 
    det_names = re.findall(_detector_pattern, f)
    det_names = np.unique(det_names)

    # Append detector reshaping
    f += "\n\n# Reshape detectors\n"
    for dn in det_names:
        if dn + "E" in det_names:
            f += "{name}.shape = ({name}_VALS, 13)\n".format(name=dn)
        else:
            f += "{name}.shape = ({name_min_E}_EBINS, 3)\n".format(name=dn, name_min_E=dn[:-1])

    # Add imports to header
    header = "import numpy as np\n\n"

    # Add header to file
    f = header + f

    # Remove semicolons
    f = _replace_semicolons(f)

    # Write the file out
    if write_py:
        if isinstance(detfile, basestring):
            new_filename = detfile.rpartition('.')[0] + '.py'
        else:
            new_filename = detfile.name.rpartition('.')[0] + '.py'
        with open(new_filename, 'w') as pyfile:
            pyfile.write(f)

    # Execute the adjusted file
    det = {}
    exec(f, {}, det)

    return det
