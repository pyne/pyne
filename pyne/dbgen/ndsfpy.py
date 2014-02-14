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
    dt = np.dtype([('parent', int), ('fission_product', int),
                   ('thermal_yield', float), ('fast_yield', float),
                   ('_14MeV_yield', float)
    ])
    dfinal = np.zeros((len(fpdata),), dtype=dt)
    for index, item in enumerate(fpdata):
        dfinal[index]['parent'] = pid
        dfinal[index]['fission_product'] = conv_to_id(item)
    thermaldata = getdata(i + 2, spdat)
    for index, item in enumerate(thermaldata):
        dfinal[index]['thermal_yield'] = conv_to_num(item)[0]
    fastdata = getdata(i + 3, spdat)
    for index, item in enumerate(fastdata):
        dfinal[index]['fast_yield'] = conv_to_num(item)[0]
    dtdata = getdata(i + 4, spdat)
    for index, item in enumerate(dtdata):
        dfinal[index]['_14MeV_yield'] = conv_to_num(item)[0]
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


def getndsfpdata(nuc_data):
    """Adds the NDS fission yields to the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    """
    doc = urlopen("https://www-nds.iaea.org/sgnucdat/c2.htm")
    dat = doc.read()
    spdat = dat.split("<table>")
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
    print('Making NDS fission product yield table.')
    getndsfpdata(nuc_data)
