#! /usr/bin/env python
"""A simple script to run all of the IPython notebooks in this directory.
"""
import os
import sys
import traceback
from tempfile import NamedTemporaryFile
import subprocess

def execpy(filename, glb=None, loc=None):
    """A function equivalent to the Python 2.x execfile statement."""
    with io.open(filename, 'r') as f:
        src = f.read()
    exec(compile(src, filename, "exec"), glb, loc)

def execipynb(filename, glb=None, loc=None):
    """A function equivalent to the Python 2.x execfile statement but for IPython 
    notebooks."""
    out = NamedTemporaryFile()
    err = NamedTemporaryFile()
    rtn = subprocess.check_call(['ipython', 'nbconvert', '--to=python', '--stdout', 
                                 filename], stdout=out, stderr=err)
    out.seek(0)
    src = out.read()
    out.close()
    err.close()
    exec(compile(src, filename, "exec"), glb, loc)

def main():
    files = os.listdir('.')
    thisfile = os.path.split(__file__)[1]
    count = 0
    nsucc = 0
    nfail = 0
    summary = ""
    for f in files:
        base, ext = os.path.splitext(f)
        if f == thisfile:
            print("thisfile")
            continue
        elif ext == '.py':
            execer = execpy
        elif ext == '.ipynb':
            execer = execipynb
        else:
            continue
        count += 1
        try:
            execer(f)
            nsucc += 1
            sys.stderr.write('.')
            sys.stderr.flush()
        except KeyboardInterrupt:
            raise
        except:
            nfail += 1
            sys.stderr.write('F')
            sys.stderr.flush()
            summary += "\nFAILURE: " + f + ":\n----------\n"
            summary += traceback.format_exc()
    if count > 0:
        summary += "\n------\nTOTAL: {0} - ".format(count)
    if nsucc > 0:
        summary += " {0} successes".format(nsucc)
    if nfail > 0:
        summary += " {0} failures".format(nfail)
    print(summary)

if __name__ == '__main__':
    main()
