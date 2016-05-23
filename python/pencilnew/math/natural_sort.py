
def natural_sort(l):
    """Method sorts list l of float numbers in strings in a natural way.

    Args:
        l:      list of strings, strings will be converted to float, then sorted
    """
    import re

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    
    return sorted(l, key = alphanum_key)
