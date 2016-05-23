def is_number(s):
    """ Checks if s is a number. """
    try:
      float(s)
      return True
    except ValueError:
      return False
