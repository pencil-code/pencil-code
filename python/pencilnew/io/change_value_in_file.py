def change_value_in_file(filename, quantity, newValue, sim=False, filepath=False, DEBUG=False):
    """ Use to change a quantity in
        - *.in
        - *.local
        - submit*, i.e. submit.sh, submit.csh, files, only works if computer is readily specified in pencilnew.io.get_systemid

    Please add further functionallity by yourself!

    Args:
        filename:   can be "run.in", "start.in", "cparam.local"
        quantity:   variable to read in from file
        sim:        put simulation object here, file will be found by filename automatically
        filepath:   normally not needed, specify here where to find the file with filename, can be a list of paths if unshure
        DEBUG:      make dry run, tell me what you would do but dont change anything!
        silent:     suppress certain output by setting True

    Returns True if successful, else False
    """

    from . import get_value_from_file

    return_value = get_value_from_file(filename, quantity, change_quantity_to=newValue, sim=sim, filepath=filepath, DEBUG=DEBUG)

    if return_value == None:
        return False
    return True
