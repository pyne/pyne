from __future__ import print_function
import os
import re
import sys
from pyne.utils import QA_warn

try:
    import urllib.request as urllib2
    from urllib.error import URLError
except ImportError:
    import urllib2
    from urllib2 import URLError

from pyne import nucname

QA_warn(__name__)

if sys.version_info[0] > 2:
    basestring = str


def grab_kaeri_nuclide(nuc, build_dir="", n=None):
    """Grabs a nuclide file from KAERI from the web and places
    it a {nuc}.html file in the build directory.

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
        nuc = nucname.name(nuc).upper()

    if n is None:
        filename = os.path.join(build_dir, nuc + ".html")
        kaeri_url = "http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc={0}".format(nuc)
    else:
        filename = os.path.join(build_dir, "{nuc}_{n}.html".format(nuc=nuc, n=n))
        kaeri_url = "http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc={0}&n={n}".format(
            nuc, n=n
        )
    print("    getting {0} and placing in {1}".format(nuc, filename))

    # Get the url
    req = urllib2.Request(kaeri_url, headers={"User-Agent": "Mozilla/5.0"})
    hdl = urllib2.urlopen(req, timeout=30.0)
    i = 1

    # try reading in the data until it works or ten times
    read_in = False
    while (not read_in) and (i <= 10):
        try:
            kaeri_html = hdl.read()
            read_in = True
        except URLError:
            hdl.close()
            i += 1
            print(
                "    getting {0} and placing in {1}, attempt {2}".format(
                    nuc, filename, i
                )
            )
            hdl = urllib2.urlopen(req, timeout=30.0)

    # Write out to the file
    with open(filename, "w") as f:
        f.write(kaeri_html)


nat_iso_regex = re.compile(
    ".*?/cgi-bin/nuclide[?]nuc=([A-Za-z]{1,2}\d{1,3}).*?[(].*?[)]"
)


def parse_for_natural_isotopes(htmlfile):
    """Parses an elemental html file, returning a set of naturally occuring isotopes."""
    nat_isos = set()
    with open(htmlfile, "r") as f:
        for line in f:
            m = nat_iso_regex.search(line)
            if m is not None:
                nat_isos.add(nucname.id(m.group(1)))
    return nat_isos


all_iso_regex = re.compile(".*?/cgi-bin/nuclide[?]nuc=([A-Za-z]{1,2}\d{1,3})")


def parse_for_all_isotopes(htmlfile):
    """Parses an elemental html file, returning a set of all occuring isotopes."""
    isos = set()
    with open(htmlfile, "r") as f:
        for line in f:
            m = all_iso_regex.search(line)
            if m is not None:
                isos.add(nucname.id(m.group(1)))
    return isos
