

def is_int(s):
    """ Checks if string s is an int. """
    try:
        a = float(s)
        b = int(s)
    except ValueError:
        return False
    else:
        return a == b
