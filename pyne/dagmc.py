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
    from pyne import _dagmc
    HAVE_DAGMC = True
except ImportError:
    warn(
        "DAGMC module is not available. Please install PyNE with DAGMC enabled.\n"
        "For more information, see: https://pyne.io/install/",
        ImportWarning
        )
    HAVE_DAGMC = False
