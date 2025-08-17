__version__ = "0.2.5"

import sys

exec("from .FortranRecordReader import FortranRecordReader")
exec("from .FortranRecordWriter import FortranRecordWriter")
exec("from . import config")