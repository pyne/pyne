import os
import warnings
import importlib.metadata

__version__ = importlib.metadata.version("pyne")

try:
    from .pyne_config import *
    from .paths import *
except ImportError:
    warnings.warn(
        "It seems that PyNE is being run from its source directory. "
        "This setup is not recommended as it may lead to unexpected behavior, "
        "such as conflicts between source and installed versions. "
        "Please run your script from outside the PyNE source tree.",
        RuntimeWarning,
    )
    raise

