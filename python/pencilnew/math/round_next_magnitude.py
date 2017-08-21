def round_next_magnitude(num):
    """Rounds a number up to its next magnitude.
    e.g. 4,56 will be rounded to 10
    """
    import numpy as np

    return np.ceil(num/10.)*10.
