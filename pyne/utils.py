from __future__ import division
import os

from distutils.dir_util import remove_tree

from pyne._utils import fromstring_split, fromstring_token, endftod,\
                        use_fast_endftod, fromendf_tok, toggle_warnings,\
                        use_warnings, fromendl_tok


class QAWarning(UserWarning):
    pass


barn_conv_dict = {
    'mb': 1E-3,
    'ub': 1E-6,
    'microbarn': 1E-6,
    'b': 1.0,
    'barn': 1.0,
    'barns': 1.0,
    'kb': 1E+3,
    'kilobarn': 1E+3,
    'cm2': 1E+24,
    'cm^2': 1E+24,
    }


def to_barns(xs, units):
    """Converts a cross section with units to barns.

    Parameters
    ----------
    xs :
        Cross section value in [units].
    units : str
        Units flag, eg 'b', 'microbarn'.

    Returns
    -------
    barn_xs :
        Cross section value in [barns].

    """
    return xs * barn_conv_dict[units.lower()]


def from_barns(xs, units):
    """Converts a cross section from barns to units.

    Parameters
    ----------
    xs :
        Cross section value in [barns].
    units : str
        Units flag, eg 'b', 'microbarn'.

    Returns
    -------
    unit_xs :
        Cross section value in [units].

    """
    return xs / barn_conv_dict[units.lower()]


#########################
### message functions ###
#########################

USE_COLOR = (os.name is 'posix')


def message(s):
    """Formats a message for printing.  If on a posix system the message will
    be in color.

    """
    head = "\033[1;32m" if USE_COLOR else "*** MESSAGE ***: "
    tail = "\033[0m" if USE_COLOR else ""
    msg = head + s + tail
    return msg


def failure(s):
    """Formats a fail message for printing.  If on a posix system the message
    will be in color.

    """
    head = "\033[1;31m" if USE_COLOR else "*** FAILURE ***: "
    tail = "\033[0m" if USE_COLOR else ""
    msg = head + s + tail
    return msg


def warning(s):
    """Formats a warning message for printing. If on a posix system the message
    will be in color.
    """
    head = "\033[1;33m" if USE_COLOR else "*** WARNING ***: "
    tail = "\033[0m" if USE_COLOR else ""
    msg = head + s + tail
    return msg

##################################
### Path manipulation routines ###
##################################


def remove(path):
    """Removes a path, or recursively a directory, or does nothing
    if path is neither a file nor a directory.

    """
    if os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        remove_tree(path, verbose=False)
    else:
        pass
