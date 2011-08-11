import os
import urllib2

from pyne import nucname


def grab_kaeri_nuclide(nuc, build_dir=""):
    """Grabs a nuclide file from KAERI from the web and places 
    it in nuc.html in the build directory.

    Parameters
    ----------
    nuc : str, int
        nuclide, preferably in name form.
    build_dir : str
        Directory to place html files in.
    """
    if not isinstance(nuc, basestring):
        nuc = nucname.name(nuc)

    filename = os.path.join(build_dir, nuc + '.html')
    print "    getting {0} and placing in {1}".format(nuc, filename)

    kaeri = urllib2.urlopen('http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc={0}'.format(nuc))
    with open(filename, 'w') as f:
        f.write(kaeri.read())
