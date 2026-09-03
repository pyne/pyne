from xdress.utils import apiname
from xdress.doxygen import default_doxygen_config

package = "pyne"
packagedir = "pyne"

plugins = (
    "xdress.stlwrap",
    "xdress.autoall",
    "xdress.autodescribe",
    "xdress.doxygen",
    "xdress.cythongen",
)

doxygen_config = {
    "PROJECT_NAME": "PYNE",
    "EXTRACT_ALL": False,  # Note usage of python False
    "GENERATE_DOCBOOK": False,
    "GENERATE_LATEX": True,  # Could be 'YES' or True
}

doxyfile_name = "./build/tally_doxyfile"

includes = ["/mnt/data/opt/dagmc/hdf5/include"]

extra_types = "extra_types"  # non-default value

stlcontainers = [
    ("set", "str"),
    ("set", "int32"),
    ("map", "str", "str"),
    ("map", "str", "int32"),
    ("map", "int32", "str"),
    ("map", "str", "uint32"),
    ("map", "uint32", "str"),
    ("map", "str", "float64"),
    ("map", "uint32", "uint32"),
    ("map", "int32", "int32"),
    ("map", "int32", "float64"),
    ("map", "int32", "complex"),
    ("map", "uint32", "float64"),
    ("map", "str", ("vector", "float64")),
    ("map", "int32", ("vector", "float64")),
]

# stlcontainers_module = 'stlcontainers'

classes = [apiname("Tally", "src/tally.*", incfiles="tally.h")]
