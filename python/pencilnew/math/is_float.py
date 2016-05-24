
def is_float(s):
    """ Checks if string s is a float. """
    try:
        a = float(s)
    except ValueError:
        return False
    else:
        return True
