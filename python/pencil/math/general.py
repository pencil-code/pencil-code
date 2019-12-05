# general.py
#
# Various short mathematical helper functions.
#
# Authors: A. Schreiber (aschreib@isaac1.bc.rzg.mpg.de)
#          S. Candelaresi (iomsn1@gmail.com).
"""
Contains the definitions various short mathematicl function.
"""


def is_number(num_str):
    """
    Checks if string s is a number.

    call signature:

    is_number(num_str):

    Keyword arguments:

    *num_str*:
      String containing the number to be checked.
    """

    try:
        float(num_str)
        return True
    except ValueError:
        return False


def is_int(num_str):
    """
    Checks if string num_str is an int.

    call signature:

    is_int(num_str):

    Keyword arguments:

    *num_str*:
      String containing the number to be checked.
    """

    try:
        a = float(num_str)
        b = int(num_str)
    except ValueError:
        return False
    else:
        return a == b


def is_float(num_str):
    """
    Checks if string num_str is a float.

    call signature:

    is_float(num_str):

    Keyword arguments:

    *num_str*:
      String containing the number to be checked.
    """

    try:
        float(num_str)
    except ValueError:
        return False
    else:
        return True


def is_iterable(i):
    """
    Checks if i is an iterable.

    call signature:

    is_iterable(i):

    Keyword arguments:

    *i*:
      Object to be checked.
    """

    try:
        iter(i)
        return True
    except TypeError:
        return False
    return False


def log_range(a, b, include_ab=False, finer_factor=1):
    """
    Returns a log-space from a to b.
    	Example: log_range(4, 40) = [4, 5, 6, 7, 8, 9, 10, 20, 30, 40].

    call signature:

    log_range(a, b, include_ab=False, finer_factor=1):

    Keyword arguments:

    *a, b*:
      Start and end values.

    *include_ab*:
      Always include a and b in the interval.

    *finer_factor*:
      Factor by which the range is being refined.
    """

    import numpy as np

    ten_range = np.arange(1, 10, 1./finer_factor)
    low_exp = int(np.log10(a))
    high_exp = int(np.log10(b))
    exp_range = range(low_exp, high_exp+1)

    # Check consistency.
    if a == b:
        return a

    if include_ab:
        logrange = [a]
    else:
        logrange = []
    for i in exp_range:
        logrange = logrange + [10.0**i*j for j in ten_range]

    if include_ab:
        logrange = logrange + [b]

    # Clean up to requested range.
    seen = set()
    return np.array([c for c in logrange if c >= a and c <= b and c not in seen and not seen.add(c)])


def round_next_magnitude(num):
    """
    Rounds a number up to the next multiple of 10, e.g. 4.56 will be rounded to 10.

    call signature:

    round_next_magnitude(num):

    Keyword arguments:

    *num*:
      Number to be rounded.
    """

    import numpy as np

    return np.ceil(num/10.)*10.


def natural_sort(string_list, reverse=False):
    """
    Sort a list of float numbers in strings in a natural way.

    call signature:

    natural_sort(string_list, reverse=False):

    Keyword arguments:

    *string_list*:
      List of strings. Will be converted to float, then sorted.

    *reverse*:
      If true, list in reverse order.
    """

    import re

    convert = lambda text: float(text) if is_number(text) else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9.]+)', str(key))]

    return sorted(string_list, key=alphanum_key, reverse=reverse)
