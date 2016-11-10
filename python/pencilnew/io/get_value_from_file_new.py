def get_value_from_file(filename, quantity, change_quantity_to=False, sim=False, filepath=False, DEBUG=False):
    """ Use to read in a quantity from
        - *.in
        - *.local
        - submit*, i.e. submit.sh, submit.csh, files, only works if computer is readily specified in pencilnew.io.get_systemid

    Please add further functionallity by yourself!

    Args:
        file:       can be "run.in", "start.in", "cparam.local"
        quantity:   variable to read in from file
        sim:        put simulation object here, file will be found by filename automatically
        filepath:   normally not needed, specify here where to find the file with filename, can be a list of paths if unshure
        DEBUG:      make dry run, tell me what you would do but dont change anything!
    """

    import os
    from os.path import join, abspath, exists
    import pencilnew
    from pencilnew.math import is_number, is_float, is_int

    # prepare filename and quantity
    filename = filename.strip()                                             # get rid of whitespaces
    quantity = quantity.strip()

    # prepare search_path list to searth filename in
    if filepath == False:
        if sim == False: sim = pencilnew.get_sim()
        search_paths = [sim.path, join(sim.path, 'src')]                    # add other search paths here!!

    elif type(filepath) == type('string'):
        if filepath.endswith(filename): filepath = filepath[:-len(filename)]    # clean filepath if filename occures to be in there at the end
        search_paths = [abspath(filepath.strip())]                              # correct path format

    elif type(filepath) == type(['list']):
        search_paths = filepath

    else:
        print('! ERROR: Filename '+str(filename)+' could not be interprated or found!'); return False

    for search_path in search_paths:
        if filename in os.listdir(search_path):
            filepath = join(search_path, filename)
            break

    if DEBUG: print('~ DEBUG: Found file in '+filepath)

    # now having absolute filepath to file, lets check that file and find quantity inside!
    with open(filename, 'r') as f: data_raw = f.readlines()

    line_matches = []               # scan through file for differently for different files

    if filename.endswith('.in') or filename == 'cparam.local':
        SYM_COMMENT = '!'
        SYM_ASSIGN = '='
        SYM_SEPARATOR = ','

        for ii, line in enumerate(data_raw):
            if quantity in line.split(SYM_COMMENT)[0]: line_matches.append(ii)

    elif filename.startswith('submit') and filename.split('.')[-1] in ['csh', 'sh']:
        SYM_COMMENT = False
        SYM_ASSIGN = '='
        SYM_SEPARATOR = ','

        for ii, line in enumerate(data_raw):
            if line.startswith('#@') or line.startswith('# @'):
                if quantity in line: line_matches.append(ii)

    else:
        print('! ERROR: Filename unknown! No parsing possible! Please enhance this function to work with '+filename)

    if len(line_matches) > 1: print('! ERROR: Found more than one line with keyword "'+quantity+'" inside!'); return False
    if len(line_matches) == 0: print('! ERROR: Found no line with keyword "'+quantity+'" inside!'); return False
    line = data_raw[line_matches[0]].replace(' ','').replace('\n', '')       # get line with quantity inside

    qs = line.partition(quantity+SYM_ASSIGN)
    if SYM_ASSIGN in qs[-1]:
        qs = qs[:2]+qs[-1].partition(SYM_ASSIGN)
        qs = qs[:2]+qs[-1].partition(SYM_ASSIGN)
        qs = qs[:2]+qs[2].rpartition(',')+qs[3:]

    qs = list(qs)
    q = qs[2]

    # cleanup of q:
    if q.startswith("'") and q.endswith("'"): q = q[1:-1]       # quantity is string in 'STRING'
    elif q.startswith('"') and q.endswith('"'): q = q[1:-1]       # quantity is string in "STRING"
    elif not is_number(q[0]): q = q.strip()
    elif is_float:
        q = float(q)
        if is_int(q): q = int(q)

    if change_quantity_to != False:
        if DEBUG: print('~ DEBUG: Would change quantity '+quantity+' from '+str(q)+' to '+str(change_quantity_to))
        if DEBUG: print('~ DEBUG: old line:'+line)

        qs[2] = str(change_quantity_to)
        new_line = ''.join(qs).replace(SYM_SEPARATOR, SYM_SEPARATOR+' ')
        if DEBUG: print('~ DEBUG: new line:'+new_line)
        if not DEBUG: print('! NOT FINISHED: export changed line')

    return q



#location_of_ASSIGNMENT = line.find(quantity)+len(quantity)+1
#line[location_of_ASSIGNMENT:].split(SYM_ASSIGN)[0].rpartition(SYM_SEPARATOR)



#qs = line.split(SYM_SEPARATOR)         # split line to separate quantity names, but also separates values if its a tuple or a string containing SYM_SEPARATOR

#for ii, q in enumerate(qs):
#    if quantity in q: break




#if filename.endswith('.in'):
#    data = [q for q in data if not q.startswith('!')]
#    data = [q for q in data if not q.startswith('/')]
#    data = [q for q in data if not q.startswith('\n')]
#    data = [q for q in data if not q.startswith('&')]
