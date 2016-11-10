def change_value_in_file(filename, quantity, newValue, sim=False, filepath=False, DEBUG=False):
    from pencilnew.io import get_value_from_file

    return_value = get_value_from_file(filename, quantity, change_quantity_to=newValue, sim=sim, filepath=filepath, DEBUG=DEBUG):

    if return_value == False: return False
    else: return True
