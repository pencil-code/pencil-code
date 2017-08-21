
def natural_sort(l, reverse=False):
    """Method sorts list l of float numbers in strings in a natural way.

    Args:
        l:      list of strings, strings will be converted to float, then sorted
    """
    import re
    from .is_number import is_number

    convert = lambda text: float(text) if is_number(text) else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9.]+)', str(key)) ]

    return sorted(l, key = alphanum_key, reverse = reverse)
