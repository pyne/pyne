import os
from warnings import warn


if os.name == "nt":
    p = os.environ["PATH"].split(";")
    lib = os.path.join(os.path.split(__file__)[0], "lib")
    os.environ["PATH"] = ";".join([lib] + p)

try:
    from .pyne_version import PYNE_VERSION

    __version__ = PYNE_VERSION
    from .pyne_config import *

except ImportError:
    msg = (
        "Error importing PyNE: you should not try to import PyNE from "
        "its source directory; please exit the PyNE source tree, and relaunch "
        "your python interpreter from there."
    )
    warn(msg, Warning)
    raise
