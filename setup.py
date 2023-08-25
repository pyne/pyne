import os
import re
import sys
import ssl
import shutil
import platform
import subprocess
from glob import glob
from distutils import sysconfig
from contextlib import contextmanager
from skbuild import setup


# import src into pythonpath - needed to actually run decaygen/atomicgen
if "." not in sys.path:
    sys.path.append(os.getcwd() + "/src")


def absexpanduser(x):
    return os.path.abspath(os.path.expanduser(x))


IS_NT = os.name == "nt"
LOCALDIR = absexpanduser("~/.local")
CMAKE_BUILD_TYPES = {
    "none": "None",
    "debug": "Debug",
    "release": "Release",
    "relwithdebinfo": "RelWithDebInfo",
    "minsizerel": "MinSizeRel",
}
ON_DARWIN = platform.system() == "Darwin"
LIBEXT = ".dylib" if ON_DARWIN else ".so"

SKIP_OPTION = "SKIP"


@contextmanager
def indir(path):
    orig = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(orig)


@contextmanager
def cleanpypath(path):
    orig = sys.path
    sys.path = [p for p in sys.path if p != path]
    yield
    sys.path = orig


def ssl_context():
    # this is compitble for both Python 2 & 3
    # on Python 3, you can do just ssl.SSLContext()
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


DECAY_H = os.path.join("src", "decay.h")
DECAY_CPP = os.path.join("src", "decay.cpp")
DECAY_S = glob("src/decay*.s")
DECAY_H_REP = os.path.join("src", "_decay.h")
DECAY_CPP_REP = os.path.join("src", "_decay.cpp")
DECAY_URL = "http://raw.githubusercontent.com/pyne/data/master/decay.tar.gz"


CRAM_H = os.path.join("src", "cram.h")
CRAM_C = os.path.join("src", "cram.c")
CRAM_S = glob("src/cram*.s")
CRAM_URL = "http://raw.githubusercontent.com/pyne/data/master/cram.tar.gz"


local_ensdf_evaluators = [
    "alphad",
    "delta",
    "gtol",
    "bldhst",
    "hsicc",
    "hsmrg",
    "seqhst",
    "logft",
    "radd",
    "ruler",
]
local_ensdf_tools = [
    ["ensdf_processing/RADD/98AK04.in", "98AK04.in"],
    ["ensdf_processing/RADD/ELE.in", "ELE.in"],
]


def copy_ensdf_executables(exe_dest):
    print("Copying ENSDF Executables to install directory")
    # Hack for copying the executables the first time PyNE is installed, before
    # pyne has been added to the python path.
    if exe_dest[-4:] != "pyne":
        exe_dest = sysconfig.get_python_lib()
        for f in os.listdir(sysconfig.get_python_lib()):
            if re.match("pyne", f):
                exe_dest = exe_dest + "/" + f
        exe_dest = exe_dest + "/pyne"
    for tool in local_ensdf_evaluators:
        try:
            local_path = os.path.join("build", os.path.join("src", tool))
            dest_path = os.path.join(exe_dest, tool)
            shutil.copy(local_path, dest_path)
        except Exception:
            print(
                "Some ENSDF processing executables were unable to be copied to the \
                   install directory."
            )
    for tool in local_ensdf_tools:
        try:
            local_path = os.path.join("src", tool[0])
            dest_path = os.path.join(exe_dest, tool[1])
            shutil.copy(local_path, dest_path)
        except Exception:
            print(
                "Some ENSDF processing executables were unable to be copied to the \
                   install directory."
            )


ATOMIC_H = os.path.join("src", "atomic_data.h")
ATOMIC_CPP = os.path.join("src", "atomic_data.cpp")
ATOMIC_H_UNDER = os.path.join("src", "_atomic_data.h")
ATOMIC_CPP_UNDER = os.path.join("src", "_atomic_data.cpp")


def generate_atomic():
    with indir("src"):
        try:
            import atomicgen
        except ImportError:
            return False
        try:
            atomicgen.main()
        except Exception:
            return False
    return True


def ensure_atomic():
    # generate the data
    generated = generate_atomic()
    if generated:
        return
    # last resort - if generate atomic failed, use the backup
    if not os.path.isfile(ATOMIC_H) and not os.path.isfile(ATOMIC_CPP):
        print("!!! Could not generate atomic data, using backup.")
        shutil.copy(ATOMIC_H_UNDER, ATOMIC_H)
        shutil.copy(ATOMIC_CPP_UNDER, ATOMIC_CPP)


def ensure_nuc_data():
    import tempfile

    tdir = tempfile.gettempdir()
    with cleanpypath("."), cleanpypath(os.getcwd()), indir(tdir):
        from pyne.dbgen import nuc_data_make
        from pyne.dbgen.api import build_dir

        bdir = os.path.join(os.getcwd(), "build", build_dir)
        nuc_data_make.main(args=["-b", bdir])



def update_other_args(ns):
    if ns.hdf5 is not None:
        os.environ["HDF5_ROOT"] = ns.hdf5
    if ns.moab is not None:
        os.environ["MOAB_ROOT"] = ns.moab
    if ns.dagmc is not None:
        os.environ["DAGMC_ROOT"] = ns.dagmc

ensure_atomic()
setup(
    name="PyNE",
    version="1.0",
    description="Your project description",
    author="Your Name",
    author_email="your@email.com",
    # packages=["pyne"],  # Add your package(s) here
    cmake_args=[
        "-DPYTHON_EXECUTABLE:FILEPATH=" + sys.executable,
        "-DCMAKE_BUILD_TYPE:STRING=Release",
    ],
    cmake_install_dir=".",
)
