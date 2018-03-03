from __future__ import division
import os

from distutils.dir_util import remove_tree

from pyne._utils import fromstring_split, fromstring_token, endftod,\
                        use_fast_endftod, fromendf_tok, toggle_warnings,\
                        use_warnings, fromendl_tok


class QAWarning(UserWarning):
    pass

time_conv_dict = {'as': 1e-18,
                  'attosec': 1e-18,
                  'attosecond': 1e-18,
                  'attoseconds': 1e-18,
                  'fs': 1e-15,
                  'femtosec': 1e-15,
                  'femtosecond': 1e-15,
                  'femtoseconds': 1e-15,
                  'ps': 1e-12,
                  'picosec': 1e-12,
                  'picosecond': 1e-12,
                  'picoseconds': 1e-12,
                  'ns': 1e-9,
                  'nanosec': 1e-9,
                  'nanosecond': 1e-9,
                  'nanoseconds': 1e-9,
                  'us': 1e-6,
                  'microsec': 1e-6,
                  'microsecond': 1e-6,
                  'microseconds': 1e-6,
                  'ms': 1e-3,
                  'millisec': 1e-3,
                  'millisecond': 1e-3,
                  'milliseconds': 1e-3,
                  's': 1.0,
                  'sec': 1.0,
                  'second': 1.0,
                  'seconds': 1.0,
                  'm': 60.0,
                  'min': 60.0,
                  'minute': 60.0,
                  'minutes': 60.0,
                  'h': 3600.0,
                  'hour': 3600.0,
                  'hours': 3600.0,
                  'd': 86400.0,
                  'day': 86400.0,
                  'days': 86400.0,
                  'w': 86400.0*7.0,
                  'week': 86400.0*7.0,
                  'weeks': 86400.0*7.0,
                  'y': 86400.0*365.25,
                  'year': 86400.0*365.25,
                  'years': 86400.0*365.25,
                  'c': 86400.0*365.25*100,
                  'century': 86400.0*365.25*100,
                  'centuries': 86400.0*365.25*100,
                  }


def to_sec(input_time, units):
    """Converts a time with units to seconds.

    Parameters
    ----------
    input_time : number
        Time value in [units].
    units : str
        Units flag, eg 'min', 'ms', 'days'

    Returns
    -------
    sec_time : float
        Time value in [sec].

    """
    conv = time_conv_dict.get(units.lower(), None)
    if conv:
        sec_time = input_time * conv
        return sec_time
    else:
        raise ValueError('Invalid units: {0}'.format(units))

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
