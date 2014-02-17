from __future__ import print_function, division
from urllib import urlopen
import numpy as np
import numpy.lib.recfunctions
from pyne import nucname
import tables as tb
from pyne.dbgen.api import BASIC_FILTERS


def readtable(i, spdat):
    """
    Reads in a set of 5 html tables and returns corresponding yield data
    """
    parent = getdata(i, spdat)[0]
    pfinal = (parent.split('<strong>')[1]).split('</strong>')[0]
    pid = conv_to_id(pfinal)
    fpdata = getdata(i + 1, spdat)
    dt = np.dtype([('from_nuc', int), ('to_nuc', int),
                   ('thermal_yield', float), ('thermal_yield_err', float),
                   ('fast_yield', float), ('fast_yield_err', float),
                   ('_14MeV_yield', float), ('_14MeV_yield_err', float)
                  ])
    dfinal = np.zeros((len(fpdata),), dtype=dt)
    for index, item in enumerate(fpdata):
        dfinal[index]['from_nuc'] = pid
        dfinal[index]['to_nuc'] = conv_to_id(item)
    thermaldata = getdata(i + 2, spdat)
    for index, item in enumerate(thermaldata):
        dat, err = conv_to_num(item)
        dfinal[index]['thermal_yield'] = dat
        dfinal[index]['thermal_yield_err'] = err
    fastdata = getdata(i + 3, spdat)
    for index, item in enumerate(fastdata):
        dat, err = conv_to_num(item)
        dfinal[index]['fast_yield'] = dat
        dfinal[index]['fast_yield_err'] = err
    dtdata = getdata(i + 4, spdat)
    for index, item in enumerate(dtdata):
        dat, err = conv_to_num(item)
        dfinal[index]['_14MeV_yield'] = dat
        dfinal[index]['_14MeV_yield_err'] = err
    return dfinal


def conv_to_id(nuc):
    """
    Converts html nuclide names to nuclide ids
    """
    parts = nuc.split('-')
    return nucname.id(parts[1] + parts[2])


def conv_to_num(dstring):
    """
    Converts html number and error to floats
    """
    if dstring == '-':
        return 0, 0
    dat, err = dstring.split('&plusmn;')
    if '<sup>' in dat:
        dat = parse_num(dat)
    else:
        dat = float(dat)
    if '<sup>' in err:
        err = parse_num(err)
    else:
        err = float(err)
    return dat, err


def parse_num(dst):
    """
    Converts html numbers with exponents to floats
    """
    nums = dst.split('x')
    base = float(nums[0])
    exp = (nums[1].split('<sup>')[1]).split('</sup>')[0]
    return base * 10 ** float(exp)


def getpoint(line):
    """
    Gets data entries from html lines
    """
    spline = line.split('<tr><td class="xl28b">&nbsp;&nbsp;')
    if len(spline) > 1:
        data = spline[1].split('</td></tr>')[0]
    else:
        data = None
    return data


def getdata(i, spdat):
    """
    Gets the data from the nds html table
    """
    lines = spdat[i].splitlines()
    dlist = []
    for line in lines:
        d = getpoint(line)
        if d is None:
            continue
        dlist.append(d)
    return dlist


def make_fpy_table(nuc_data, build_dir=""):
    """Adds the NDS fission yields to the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    """
    build_filename = os.path.join(build_dir, 'wimsd-fpyield.html')
    with open(build_filename, 'r') as f:
        raw_data = f.read()
    spdat = raw_data.split("<table>")
    alldata = []
    for i in range(1, 31, 5):
        alldata.append(readtable(i, spdat))
    alldata = numpy.lib.recfunctions.stack_arrays(alldata)
    db = tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS)
    if not hasattr(db.root, 'neutron'):
        neutron_group = db.createGroup('/', 'neutron', 'Neutron Data')
    fpy_table = db.createTable('/neutron/', 'nds_fission_products', alldata,
                               'WIMSD Fission Product Yields, percent [unitless]')
    fpy_table.flush()
    db.close()

def grab_fpy(build_dir="", file_out='nds-fpyield.html'):
    """Grabs the NDS fission product yields from the IAEA website
    """
    build_filename = os.path.join(build_dir, file_out)
    local_filename = os.path.join(os.path.dirname(__file__), file_out)

    if os.path.exists(local_filename):
        shutil.copy(local_filename, build_filename)
        return

    nist = urllib2.urlopen('https://www-nds.iaea.org/sgnucdat/c2.htm')
    with open(build_filename, 'w') as f:
        f.write(nist.read())


def make_fpy(args):
    """Controller function for NDS fission products."""
    nuc_data, build_dir = args.nuc_data, args.build_dir
    # Check that the table exists
    with tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS) as f:
        if hasattr(f.root, 'neutron') and hasattr(f.root.neutron,
                                                  'nds_fission_products'):
            print('skipping WIMSD fission product yield table creation; '
                  'already exists.')
            return
    print("Grabbing NDS fission product yield data.")
    grab_fpy(build_dir)

    print('Making NDS fission product yield table.')
    make_fpy_table(nuc_data, build_dir)
