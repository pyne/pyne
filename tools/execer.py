#! /usr/bin/env python
"""A simple script to run all of the Jupyter notebooks in this directory.
"""
import io
import os
import sys
import traceback
from tempfile import NamedTemporaryFile
import subprocess


class get_ipython:
    """Mocks some ipython funniness, such as magic commands."""

    def magic(self, *args, **kwargs):
        pass


def execpy(filename, glb=None, loc=None):
    """A function equivalent to the Python 2.x execfile statement."""
    glb = {} if glb is None else glb
    with io.open(filename, "r") as f:
        src = f.read()
    exec(compile(src, filename, "exec"), glb, loc)


def execipynb(filename, glb=None, loc=None):
    """A function equivalent to the Python 2.x execfile statement but for IPython
    notebooks."""
    glb = {} if glb is None else glb
    glb["get_ipython"] = get_ipython
    out = NamedTemporaryFile()
    err = NamedTemporaryFile()
    env = os.environ.copy()
    env["TERM"] = "dumb"
    env["PYTHONIOENCODING"] = "utf-8"
    rtn = subprocess.check_call(
        ["jupyter", "nbconvert", "--to=python", "--stdout", filename],
        stdout=out,
        stderr=err,
        env=env,
    )
    out.seek(0)
    src = out.read()
    out.close()
    err.close()
    exec(compile(src, filename, "exec"), glb, loc)


def main():
    cwd = os.getcwd()
    execerdir = os.path.dirname(os.path.abspath(__file__))
    sys.path = [cwd if d == execerdir else d for d in sys.path]
    files = sorted(os.listdir(cwd))
    thisfile = os.path.split(__file__)[1]
    count = 0
    nsucc = 0
    nfail = 0
    summary = ""
    for f in files:
        base, ext = os.path.splitext(f)
        if f == thisfile:
            continue
        elif ext == ".py":
            execer = execpy
            cat = "cat " + f
        elif ext == ".ipynb":
            execer = execipynb
            cat = "jupyter nbconvert --to=python --stdout " + f
        else:
            continue
        count += 1
        try:
            execer(f)
            nsucc += 1
            sys.stderr.write(".")
            sys.stderr.flush()
        except KeyboardInterrupt:
            raise
        except:
            nfail += 1
            sys.stderr.write("F")
            sys.stderr.flush()
            summary += "\nFAILURE: " + f + ":\nCAT: " + cat + "\n----------\n"
            summary += traceback.format_exc()
    if count > 0:
        summary += "\n------\nTOTAL: {0} - ".format(count)
    if nsucc > 0:
        summary += " {0} successes".format(nsucc)
    if nfail > 0:
        summary += " {0} failures".format(nfail)
    print(summary)
    if nfail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
