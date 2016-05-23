

def is_int(x):
    """ Checks if x is an int. """
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b
