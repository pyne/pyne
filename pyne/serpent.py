import re
import sys
from warnings import warn
from pyne.utils import QAWarning
 
import numpy as np

if sys.version_info[0] > 2:
    basestring = str
    
warn(__name__ + " is not yet QA compliant.", QAWarning)

_if_idx_str_serpent1 = (
    'if (exist("idx", "var"));\n'
    '  idx = idx + 1;\n'
    'else;\n'
    '  idx = 1;\n'
    'end;'
    )

_if_idx_str_serpent2 = (
    "if (exist('idx', 'var'));\n"
    '  idx = idx + 1;\n'
    'else;\n'
    '  idx = 1;\n'
    'end;'
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

_detector_pattern_all = r"(DET\w+)\s*=\s*"

_imaterial_line_pattern = r"(i[a-zA-Z]\w+)\s*=\s*\d*;"

def _delete_imaterial(s):
    """"Remove imaterial information from the top of Serpent2 *_dep.m file started from 'i' with nothing."""
    s = re.sub(_imaterial_line_pattern,
               '', s)
    return s

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
    s = re.sub(_numpy_array_pattern,
               lambda mo: 'np.array(' + mo.group(0) + ')', s)

    return s


def parse_res(resfile, write_py=False):
    """Converts a serpent results ``*_res.m`` output file to a dictionary (and
    optionally to a ``*_res.py`` file).

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
    if_idx_str = _if_idx_str_serpent1
    IDX = f.count(if_idx_str)
    if IDX == 0:
        if_idx_str = _if_idx_str_serpent2
        IDX = f.count(if_idx_str)

    # Replace if statements with something more meaningful
    fpart = f.partition(if_idx_str)
    f = fpart[0] + "idx = 0" + fpart[2]
    f = f.replace(if_idx_str, 'idx += 1')

    # Replace matlab Arrays
    f = _replace_arrays(f)

    # Add imports to header
    header = "import numpy as np\n\n"

    # Find all variables and shape
    vars_shape = np.array(list(set(re.findall(_lhs_variable_pattern, f))))
    vars_dtype = dict(re.findall(_rhs_variable_pattern, f))
    # Initialize variables to zero
    header = header + "# Initialize variables\n"
    for vs in vars_shape:
        # Determine shape
        s = re.search(r'\[.*?:(.*?)\]', vs[1])
        if s is None:
            vs_shape = ""
        else:
            vs_shape = s.group(1)
            vs_shape = vs_shape.split()
            vs_shape = ", ".join(vs_shape)

        # Determine Data type
        rhs = vars_dtype[vs[0]]
        if ("\'" in rhs) or ("\"" in rhs):
            dt = "'S{0}'".format(int(s.group(1)))
            vs_shape = ""
        elif ('.' in rhs) or ('E' in rhs) or ('e' in rhs):
            dt = "float"
        else:
            dt = "int"

        zero_line = "{0} = np.zeros([{1}, {2}], dtype={3})\n".format(vs[0],
                     IDX, vs_shape, dt)
        header = header + zero_line

    # Add IDx to file
    if 0 < IDX:
        header = header + "\n\n# Maximum Index\n\nIDX = {0}\n\n".format(IDX)

    # Add header to file
    f = header + f

    # Replace variable overrides
    vars = np.array(list(set(re.findall("(" + _lhs_variable_pattern + ")", f))))
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
    exec(f, res, res)
    if '__builtins__' in res:
        del res['__builtins__']
    return res


def parse_dep(depfile, write_py=False, make_mats=True):
    """Converts a serpent depletion ``*_dep.m`` output file to a dictionary (and
    optionally to a ``*_dep.py`` file).

    Parameters
    ----------
    depfile : str or file-like object
        Path to depletion file or a dep file handle.
    write_py : bool, optional
        Flag for whether to write the dep file to an analogous python file.
    make_mats : bool, optional
        Flag for whether or not to build Materials out of mass data and add
        these to the return dictionary.  Materials so added have names which
        end in '_MATERIAL'.

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

    # Remove imaterial information from the top of Serpent2 *_dep.m file
    f = _delete_imaterial(f)
    
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
        mat_gen_line = "{name}MATERIAL = [{name}VOLUME * Material(dict(zip(zai[:-2], {name}MDENS[:-2, col]))) for col in cols]\n"
        footer += ('\n\n# Construct materials\n'
                   'zai = list(map(int, ZAI))\n'
                   'cols = list(range(len(DAYS)))\n')
        base_names = re.findall('(MAT_\w*_)MDENS = ', f)
        for base_name in base_names:
            footer += mat_gen_line.format(name=base_name)
        footer += "TOT_MATERIAL = [Material(dict(zip(zai[:-2], TOT_MASS[:-2, col]))) for col in cols]\n"
        footer += "del zai, cols\n"

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
    exec(f, dep, dep)
    if '__builtins__' in dep:
        del dep['__builtins__']
    return dep


def parse_det(detfile, write_py=False):
    """Converts a serpent detector ``*_det.m`` output file to a dictionary (and
    optionally to a ``*_det.py`` file).

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
    all_det_names = re.findall(_detector_pattern_all, f)
    all_det_names = np.unique(all_det_names)

    is_serpent_1 = any([(dn.endswith('_VALS') and dn[:-5] in det_names) or
                        (dn.endswith('_EBINS') and dn[:-6] in det_names)
                        for dn in all_det_names])

    # Append detector reshaping
    f += '\n\n# Reshape detectors\n'
    for dn in det_names:
        if is_serpent_1:
            if dn + 'E' in det_names:
                f += '{name}.shape = ({name}_VALS, 13)\n'.format(name=dn)
            else:
                f += '{name}.shape = ({name_min_E}_EBINS, 3)\n'.format(name=dn,
                                                            name_min_E=dn[:-1])
        else:
            if (dn + 'T' in det_names):
                f += '{name}.shape = (len({name})//13, 13)\n'.format(name=dn)
            elif (dn + 'E' in det_names):
                f += '{name}.shape = (len({name})//12, 12)\n'.format(name=dn)
            else:
                f += '{name}.shape = (len({name})//3, 3)\n'.format(name=dn)

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

def parse_coe(coefile, grabB1=True):
    """ Used to parse output from serpent 2's ability to branch problem variables
    in order to generate group constants at several states. 
    
    See the serpent wiki at:
    http://serpent.vtt.fi/mediawiki/index.php/Automated_burnup_sequence

    Example:
        say that you apply this to a serpent output with fuel temperature changes
        named fuel0 through fuel14. Then, control rod positions were altered with
        branches rod0 through rod4. In total, this would make 75 sets of different
        group constants since all combinations are covered. To access INF_KINF with
        (ie k_infinity with no B1 corrections) in universe 2
        fuel0, rod0 at day 0 of burnup, you'd use:

        coe[1]         ["2"]     ["fuel0"]["rod0"]["INF_KINF"]
        #   ^ dep step  ^ uni ID  ^branches       ^ parameter name
        
        all parameters serpent gives are available by default

    Parameters
    ----------
    coefile : str or file-like object
        Path to *.coe file or a *.coe file handle.

    grabB1 : bool
        whether to include results from the B1 corrected spectrum
        if group constants were generated on full-core, set this to false
        to save memory.

    Returns
    -------
    coe : dict
        Nested dictionary of the parsed coe file. First dimension key is burnup index (1-indexed). Next
        is the universe group constants were generated in. Then, branch names represent
        the keys, in the order that serpent printed them. On the final dimension, the parameters returned by serpent
        make up the remaining keys. 

        note: B1 is used for leakage-corrected group constants generated on an infinite
        lattice. Don't use these values otherwise.

    """
    if isinstance(coefile, basestring):
        f = open(coefile, 'r')
    else:
        f = coefile

    coe = {}
    universes = [] # universes GCs were made in
    layers = None

    while 1:

        # this line describes branch and coe indices
        try:
            l = next(f).split()
        except StopIteration:
            break
        _,numbranch,coeIndex,nTot,nUni = tuple([int(item) for item in l])

        # this line describes what branch you're on
        l = next(f).split()
        # init list to hold branch names if it is still None
        if layers == None:
            layers = [ [] for i in range(int(l[0]))]
        for i,branch in enumerate(l[1:]):
            # append if not seen before
            if branch not in layers[i]:
                layers[i].append(branch)
        theseBranches=l[1:]

        next(f)
        l = next(f).split()
        nBU = int(l[2]) # number of burnup steps
        BUi = int(l[1]) # index of BU step
        BU  = float(l[0]) # burnup value
        if BUi not in coe.keys():
            coe[BUi] = {}

        #coe[BUi]['burnDays'] = BU
        # the above may be nice to have, but it goes against
        # the paradigm here of having only nested dicts
        # until the bottom level with data.

        # ok, should be at a totally new set of data now,
        # so, let's gather all params.
        for iUni in range(nUni):

            l = next(f).split()
            uni = l[0]
            if uni not in coe[BUi].keys():
                coe[BUi][uni] = {}
            nPara = int(l[1]) # num of params to follow

            thisCoe = coe[BUi][uni]
            # delve through the branches, reassigning the thisCoe reference
            # at each successive level
            for i,branch in enumerate(theseBranches):
                if branch not in thisCoe.keys():
                    thisCoe[branch] = {}
                if i==len(theseBranches)-1:
                    break
                thisCoe = thisCoe[branch]

            names = []
            values = [] # should all be numeric in nature
            for i in range(nPara):
                # build from lists to avoid frequently rehashing dicts
                l = next(f).split()
                name = l[0]
                # check if B1
                # also, yeah 'and' could be used, but bool checking
                # is faster than string comparison, so that's why to nest
                if not grabB1:
                    if name[0:2] == 'B1':
                        continue
                value = l[2:]
                names.append(name)
                if len(values) ==1:
                    values.append(float(value[0]))
                else:
                    values.append([float(item) for item in value])

            thisCoe[branch]=dict(zip(names,values))

    return coe

def print_coe(coedict,lvl=0,short=True):
    """ Pretty prints a dictionary returned
    by parse_coe. Use this to visualize its
    nested structure.
    
    Parameters:
    -----------
    coedict: dict
        particularly, one returned by parse_coe,
        although it will work on any nested dict.
    short:   bool
        short print, only prints first few params
        found. Set to false for a full STDOUT.
        
    Returns:
    --------
    None
    """
    # blah blah, recursive functions aren't fast
    # print header:
    #if lvl==0:
    #    print("BU i| uni\# | branches ...")
    keys = list(coedict.keys())
    if isinstance(coedict[keys[0]],dict):
        for k in keys:
            print(str('    '*lvl) + str(k))
            print_coe(coedict[k],lvl=lvl+1)
    else:
        if short:
            keys = keys[:5]
        for k in keys:
            print(str('    '*lvl) + str(k) + ' : ' + str(coedict[k]))
    return None
