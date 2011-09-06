import os
import re
import urllib2

from pyne import nucname


def grab_kaeri_nuclide(nuc, build_dir="", n=None):
    """Grabs a nuclide file from KAERI from the web and places 
    it in nuc.html in the build directory.

    Parameters
    ----------
    nuc : str, int
        nuclide, preferably in name form.
    build_dir : str, optional
        Directory to place html files in.
    n : None or int
        Optional flag on data to grab.  None = basic data, 
        2 = cross section summary, 3 = cross section graphs.
    """
    if not isinstance(nuc, basestring):
        nuc = nucname.name(nuc)

    if n is None:
        filename = os.path.join(build_dir, nuc + '.html')
        kaeri_url = 'http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc={0}'.format(nuc)
    else:
        filename = os.path.join(build_dir, '{nuc}_{n}.html'.format(nuc=nuc, n=n))
        kaeri_url = 'http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc={0}&n={n}'.format(nuc, n=n)
    print "    getting {0} and placing in {1}".format(nuc, filename)

    # Get the url 
    req = urllib2.Request(kaeri_url, headers={'User-Agent': 'Mozilla/5.0'})
    hdl = urllib2.urlopen(req)
    with open(filename, 'w') as f:
        f.write(hdl.read())

nat_iso_regex = re.compile('.*?/cgi-bin/nuclide[?]nuc=([A-Za-z]{1,2}\d{1,3}).*?[(].*?[)]')

def parse_for_natural_isotopes(htmlfile):
    """Parses an elemental html file, returning a set of naturally occuring isotopes."""
    nat_isos = set()
    with open(htmlfile, 'r') as f:
        for line in f:
            m = nat_iso_regex.search(line)
            if m is not None:
                nat_isos.add(nucname.zzaaam(m.group(1)))
    return nat_isos


all_iso_regex = re.compile('.*?/cgi-bin/nuclide[?]nuc=([A-Za-z]{1,2}\d{1,3})')

def parse_for_all_isotopes(htmlfile):
    """Parses an elemental html file, returning a set of all occuring isotopes."""
    isos = set()
    with open(htmlfile, 'r') as f:
        for line in f:
            m = all_iso_regex.search(line)
            if m is not None:
                isos.add(nucname.zzaaam(m.group(1)))
    return isos

