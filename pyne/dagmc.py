"""
This module attempts to import the DAGMC module from the PyNE library.

DAGMC (Direct Accelerated Geometry Monte Carlo) is a set of interfaces 
developed by the University of Wisconsin-Madison which allows Monte Carlo 
radiation transport codes to use CAD-based geometries. It is used in PyNE 
to enable complex geometry definitions.

To learn more about DAGMC, visit: https://svalinn.github.io/DAGMC/

For more information on installing PyNE with DAGMC, visit: https://pyne.io/install/
"""

from warnings import warn

try:
    from pyne._dagmc import *
    HAVE_DAGMC = True
except ImportError:
    HAVE_DAGMC = False
    msg =(
    "\n\n\033[1m\033[91mDAGMC module is not available.\033[0m\n"
    "\033[1mPlease install PyNE with DAGMC enabled.\033[0m\n"
    "For more information, see: \033[94mhttps://pyne.io/install\033[0m\n"
    )
    warn(msg, Warning)
