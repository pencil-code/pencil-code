
def is_float(x):
    """ Checks if x is a float. """
    try:
        a = float(x)
    except ValueError:
        return False
    else:
        return True
