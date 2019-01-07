"""Some PyNE testing utilities."""
from __future__ import print_function
import os
import sys
import ssl
from hashlib import md5
try:
    import urllib.request as urllib
except ImportError:
    import urllib2 as urllib


def ssl_context():
    # this is compitble for both Python 2 & 3
    # on Python 3, you can do just ssl.SSLContext()
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def download_file(url, localfile, md5_hash):
    """Donwload a file and make sure its MD5 hash matches."""
    if os.path.isfile(localfile):
        with open(localfile, "rb") as f:
            html = f.read()
    else:
        msg = 'Downloading {0!r} to {1!r}'.format(url, localfile)
        print(msg, file=sys.stderr)
        req = urllib.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        f = urllib.urlopen(req, timeout=120.0, context=ssl_context())
        try:
            html = f.read()
        finally:
            f.close()
        with open(localfile, 'wb') as f:
            f.write(html)
    obs_hash = md5(html).hexdigest()
    if obs_hash != md5_hash:
        raise AssertionError("{} hash check failed; please try redownloading "
                             "the data file.".format(localfile))
