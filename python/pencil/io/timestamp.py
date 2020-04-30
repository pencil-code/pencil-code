def timestamp():
    import datetime
    str = '%s' % datetime.datetime.now()
    return str.replace(' ', '_')
