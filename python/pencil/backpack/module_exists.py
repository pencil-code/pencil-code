def module_exists(MOD):
    """ Returns True if module MOD exists, else False. """

    import importlib

    found = importlib.util.find_spec(MOD) is not None
    return found
