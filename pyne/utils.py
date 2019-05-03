from __future__ import division
import os

from distutils.dir_util import remove_tree
import filecmp
from io import open

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


def str_to_unicode(s):
    """
    This function convert a str from binary or unicode to str (unicode).
    If it is a list of string, convert every element of the list.

    Parameters:
    -----------
    s : str or list of str

    Returns:
    --------
    s : text str or list of unicode str
    """
    if isinstance(s, str) or isinstance(s, bytes):
        # it is a str, convert to text str
        try:
            s = s.decode('utf-8')
        except:
            pass
        return s
    else:
        for i, item in enumerate(s):
            try:
                s[i] = item.decode('utf-8')
            except:
                pass
        return s    


def is_close(a, b, rel_tol=1e-9, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def is_float(s):
    """
    This function checks whether a string can be converted as a float number.
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def str_almost_same(s1, s2, rel_tol=1e-9):
    """
    This function is used to compare two string to check whether they are
    almost the same.
    Return True if two strings are exactly the same.
    Return True if two strings are almost the same with only slight difference
    of float decimals.
    Return False if two strings are different.
    """
    # if string can be converted to float number
    if is_float(s1) and is_float(s2):
        return is_close(float(s1), float(s2), rel_tol)
    else:
        # not a number
        return s1 == s2


def line_almost_same(l1, l2, rel_tol=1e-9):
    """
    This function is used to compare two lines (read from files). If they are
    the same, or almost the same (with only silight different on float
    numbers), return True. Ohterwise, return False.

    Parameters:
    -----------
    l1 : str
        Line 1
    l2 : str
        Line 2
    rel_tol : float
        Reletive tolerance for float comparison

    Returns:
    --------
    True, if two lines are the same. False, if they are different.
    """
    if l1 == l2:
        # exactly the same
        return True
    else:
        # There are differences
        ls1 = l1.strip().split()
        ls2 = l2.strip().split()
        if len(ls1) != len(ls2):
            return False
        else:
            # compare string elements of the line
            for i in range(len(ls1)):
                if str_almost_same(ls1[i], ls[2], rel_tol):
                    pass
                else:
                    return False
        return True


def file_almost_same(f1, f2, rel_tol=1e-9):
    """
    For some reasones, it's useful to compare two files that are almost the
    same. Two files, f1 and f2, the text contents are exactly the same, but
    there is a small difference in numbers. Such as the difference between
    'some text 9.5' and 'some text 9.500000000001'.
    For example, in PyNE test files, there are some expected file generated
    by python2, however, the the file generated by python3 may have difference
    in dicimals.

    Parameters:
    -----------
    f1 : str
        Path of file 1
    f2 : str
        Path of file 2
    rel_tol : float
        Reletive tolerance for float numbers

    Returns:
    True : bool
        If two file are exactly the same, or almost the same with only dicimal
        differences.
    False : bool
        If the difference the two files are more than dicimal differences.
    """
    if filecmp.cmp(f1, f2):
        # precheck
        return True
    else:
        # read line by line and then compare
        with open(f1, 'r') as f:
            lines1 = f.readlines()
            lines1 = str_to_unicode(lines1)
        with open(f2, 'r') as f:
            lines2 = f.readlines()
            lines2 = str_to_unicode(lines2)
        # check length of lines
        if len(lines1) != len(lines2):
            return False
        # check content line by line
        for i in range(len(lines1)):
            if line_almost_same(lines1[i], lines1[i], rel_tol):
                pass
            else:
                return False
    return True

