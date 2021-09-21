def in_ipynb():
    """Returns True if executed in jupyter notebook, else False."""

    try:
        return "ZMQ" in str(get_ipython())
    except:
        return False

    # try:
    #     cfg = get_ipython().config
    #     if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
    #         return True
    #     else:
    #         return False
    # except NameError:
    #     return False
