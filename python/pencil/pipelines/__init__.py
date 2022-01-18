"""
Read data and parameters from data directories.
"""

try:
    import f90nml
    from .clonesimulations import clone_sims
    from .makedict import parameter_table, make_sims_dict
except:
    print(
    "Warning: pipelines require f90nml be added to library with \
    'pip3 install f90nml' (Python 3) or \
    'pip install f90nml' (Python 2)."
    )

