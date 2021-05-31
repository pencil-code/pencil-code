def get_value_from_file(filename, quantity, change_quantity_to=None, sim=False, filepath=False, DEBUG=False, silent=False):
    """ Use to read in a quantity from
        - *.in
        - *.local
        - submit*, i.e. submit.sh, submit.csh, files, only works if computer is readily specified in pc.io.get_systemid

    Please add further functionallity by yourself!

    Args:
        filename:   can be "run.in", "start.in", "cparam.local", path to that file is extraced from filepath or sim object
        quantity:   variable to read in from file
        sim:        put simulation object here, file will be found by filename automatically
        filepath:   normally not needed, specify here where to find the file with filename, can be a list of paths if unshure
        DEBUG:      make dry run, tell me what you would do but dont change anything!
        silent:     suppress certain output by setting True

    Return:
        Returns None if not successful
    """

    import os
    import numpy as np
    from os.path import join, abspath, exists, split, isfile
    from pencil import get_sim
    from pencil.math import is_number, is_float, is_int
    from pencil.io import timestamp, debug_breakpoint, mkdir
    import re
    import copy

    def string_to_tuple(s):
        q = s.split(',')

        if is_number(q[0]):
            q = np.array([float(t) for t in q])
            q_type = 'TUPLE_FLOAT'
            return q, q_type

        if q[0] == 'T' or q[0] == 'F':
            q = np.array([bool(t=='T') for t in q])
            q_type = 'TUPLE_BOOL'
            return q, q_type

        if type(q[0]) == type('string'):
            q = [t.replace('"','').replace("'", '') for t in q]
            q_type = 'TUPLE_STRING'
            return q, q_type

        print('! ERROR: Could not parse string '+s+' into a tuple!')
        print('! DEBUG_BREAKPOINT AKTIVATED - check out the following variables: string s, tuple q, first entry in tuple q[0]')
        debug_breakpoint(); return None, None

    def tuple_to_string(t, q_type):
        return ','.join([str(a) for a in t])


    ######## prepare filename and quantity
    filename = filename.strip()                                             # get rid of whitespaces
    quantity = quantity.strip()
    q_type = False                                                          # q_type will store the type of the quantity value once found and identified

    split_filename = split(filename)
    if sim == False and split_filename[0] != '' and filepath == False:
        filepath = split_filename[0]
        filename = split_filename[1]

    ######## find correct file
    # prepare search_path list to search filename in
    if filepath == False:
        if sim == False:
            sim = get_sim()
        else:
            filepath = sim.path
        search_paths = [sim.path, join(sim.path, 'src')]                    # add other search paths here!!

    elif type(filepath) == type('string'):
        if filepath.endswith(filename): filepath = filepath[:-len(filename)]    # clean filepath if filename occures to be in there at the end
        search_paths = [abspath(filepath.strip())]                              # correct path format

    elif type(filepath) == type(['list']):
        search_paths = filepath

    else:
        print('! ERROR: Filename '+str(filename)+' could not be interprated or found!'); return None

    absolute_filepath = None
    for search_path in search_paths:
        tmp_path = join(search_path, filename)
        if os.path.isfile(tmp_path):
            absolute_filepath = tmp_path
            break

    # Traps the case of not being able to find the file
    if absolute_filepath is None:
        if DEBUG:
            print('~ DEBUG: File {0} not found in {1}!'.format(filename, search_paths))
        return None


    ######## open file
    # now having absolute filepath to file, lets check that file and find quantity inside!
    if DEBUG: print('~ DEBUG: Found file {0} in {1}'.format(filename,filepath))

    with open(absolute_filepath, 'r') as f:
        data_raw = f.readlines()


    ######## find line in file which quantity in
    line_matches = []
    # scan through file for differently for different files
    if filename.endswith('.in') or 'cparam.local' or 'Makefile.local' in filename :
        FILE_IS = 'IN_LOCAL'
        SYM_COMMENT = '!'
        SYM_ASSIGN = '='
        SYM_SEPARATOR = ','

        for ii, line in enumerate(data_raw):
            if line.strip().startswith('&'):
                continue       # filter out lines with &something, e.g. &density_run_pars
            quantity_match_tmp = re.search('[^0-9a-zA-Z_]*{0}[^0-9a-zA-Z_]'.format(quantity),
                                           line.split(SYM_COMMENT)[0])
            # Check if this substring occurs as a string.
            if quantity_match_tmp:
                quantity_match = quantity_match_tmp
                if str.count(line[0:quantity_match.start()], "'") % 2 == 0 and \
                   str.count(line[0:quantity_match.start()], '"') % 2 == 0:
                    if line_matches:
                        line_matches[0] = ii
                    else:
                        line_matches.append(ii)

    elif filename.startswith('submit') and filename.split('.')[-1] in ['csh', 'sh']:
        FILE_IS = 'SUBMIT'
        SYM_COMMENT = False
        SYM_ASSIGN = '='
        SYM_SEPARATOR = ','

        for ii, line in enumerate(data_raw):
            if line.replace(' ', '').startswith('#@') and quantity in line:
                quantity_match_tmp = re.search('[^0-9a-zA-Z_]*{0}[^0-9a-zA-Z_]'.format(quantity),
                                               line.split(SYM_COMMENT)[0])
                if quantity_match_tmp:
                    quantity_match = quantity_match_tmp
                    if line_matches:
                        line_matches[0] = ii
                    else:
                        line_matches.append(ii)
    else:
        print('! ERROR: Filename unknown! No parsing possible! Please enhance this function to work with '+filename)

    if len(line_matches) > 1: print('! ERROR: Found more than one line with keyword "'+quantity+'" inside!'); return None
    if len(line_matches) == 0:
        if silent == False: print('! ERROR: Found no line with keyword "'+quantity+'" inside '+join(filepath, filename)+'!')
        return None

    filename = os.path.basename(filename)


    ######## get line with quantity inside
    line = data_raw[line_matches[0]]

    ######## do separation of quantity from rest of line, i.e. get rid of comments and other quantities defined in this line
    comment = ''
    if SYM_COMMENT:
        tmp = line.partition(SYM_COMMENT)                                       # strip away comment
        line = tmp[0]
        if tmp[-1] != '': comment = SYM_COMMENT+tmp[-1]                         # and store for later

#    line = line.replace(' ','').replace('\n', '')                               # do cleanup in this line

    # Find the position where the quantity is stored.
    pos_equal_sign_left = quantity_match.end() + str.find(line[quantity_match.end()-1:], '=')
#    pos_equal_sign_right = pos_equal_sign_left + str.find(line[pos_equal_sign_left:], ', *[0-9a-zA-Z][0-9a-zA-Z]* *= *[0-9a-zA-Z][0-9a-zA-Z]*')
    pos_equal_sign_right = pos_equal_sign_left + str.find(line[pos_equal_sign_left:], '=')
    if pos_equal_sign_right < pos_equal_sign_left:
        pos_equal_sign_right = -1
        pos_right_comma = -1
    else:
        pos_right_comma = str.rfind(line[:pos_equal_sign_right], ',')

    # Change the quantity in the line string.
    q = copy.copy(line[pos_equal_sign_left:pos_right_comma])
#    qs = line.partition(quantity+SYM_ASSIGN)
#    if SYM_ASSIGN in qs[-1]:
#        qs = qs[:2]+qs[-1].partition(SYM_ASSIGN)
#        #qs = qs[:2]+qs[-1].partition(SYM_ASSIGN)
#        qs = qs[:2]+qs[2].rpartition(',')+qs[3:]
#
#    qs = list(qs)
#    q = qs[2]

    while q.endswith('\t'): q = q[:-1]; comment = '\t'+comment                  # take care of trailing tabulator
    while q.endswith(','): q = q[:-1]                                           # remove trailing ,


    ######## do a cleanup of quantity value q and convert into string, float, int or array, also remember data type of q
    if q.startswith("'") and q.endswith("'"):                                   # quantity q is string in form of 'STRING'
        q = q[1:-1]
        q_type = 'STRING'

    elif q.startswith('"') and q.endswith('"'):                                 # quantity q is string in form of "STRING"
        q = q[1:-1]
        q_type = 'STRING'

    elif not is_number(q[0]):                                                      # quantity q is string in form of not beeing a number
        q = q.strip().replace('"','').replace("'", '')
        q_type = 'STRING'

    try:
        float(q)
        q_type = 'FLOAT'
        if is_int(q):
            q = int(q)
            q_type = 'INT'
    except:
        if type(q) == type('string') and ',' in q:
            q, q_type = string_to_tuple(q)                                          # q is a TULPE_something
            print('q = {0}, q_type = {1}'.format(q, q_type))

        if type(q) == type('string') and q in ['F', 'f']:                            # q is BOOL
            q = False
            q_type = 'BOOL'

        if type(q)== type('string') and q in ['T', 't']:
            q = True
            q_type = 'BOOL'

        if type(q) == type('string'):
            if is_number(q[0]):
                q_type = 'STRING'

    if q_type == False:       # catch if type of q was not recognized
        print('! ERROR: Couldnt identify the data type of the quantity value: '+str(q))
        DEBUG = True
        debug_breakpoint()
    elif DEBUG:
        print('~ DEBUG: q_type = '+q_type)
    if q_type == 'FLOAT':
        q = float(q)
    elif q_type == 'INT':
        q = int(q)


    ######## if value of quantity has to be changed do:
    if change_quantity_to != None:


        ####### prepare change_quantity_to for string injection
        if q_type == 'STRING':
            if not FILE_IS=='SUBMIT':
                change_quantity_to = "'"+change_quantity_to+"'"

        elif q_type == 'BOOL':
            change_quantity_to = bool(change_quantity_to in ['T', 't', True])
            if change_quantity_to == True:
                change_quantity_to = 'T'
            elif change_quantity_to == False:
                change_quantity_to = 'F'
            else:
                print('! ERROR: There is something deeply wrong here! change_quantity_to should be bool...')
                debug_breakpoint(); return None

        elif q_type == 'FLOAT':
            change_quantity_to = '%e' % change_quantity_to

        elif q_type.startswith('TUPLE'):
            if q_type.endswith('BOOL'):
                if type(change_quantity_to) == type(['list', 'of', 'bool', 'or', 'strings']):
                    for ii, val in enumerate(change_quantity_to):
                        if val in  ['T', 't', True]:
                            change_quantity_to[ii] = 'T'
                        elif val in  ['F', 'f', False]:
                            change_quantity_to[ii] = 'F'
                        else:
                            print('! ERROR: There is something deeply wrong here! change_quantity_to['+str(ii)+'] should be bool or string representation, but it is '+str(change_quantity_to[ii]))
                            debug_breakpoint(); return None
                change_quantity_to = ','.join([str(t) for t in change_quantity_to])
            if q_type.endswith('FLOAT'):
                change_quantity_to = str(list(change_quantity_to))[1:-1]
            if q_type.endswith('STRING'):
                change_quantity_to = str(list(change_quantity_to))[1:-1]

        if DEBUG: print('~ DEBUG: Would change quantity '+quantity+' from '+str(q)+' to '+str(change_quantity_to))
        q = str(change_quantity_to)

        ######## further formatting
        new_line = line[:pos_equal_sign_left] + q + line[pos_right_comma:] + '\t'+comment    # create new line and add comment stripped away before
#        new_line = ''.join(qs).replace(SYM_SEPARATOR, SYM_SEPARATOR+' ')+'\t'+comment    # create new line and add comment stripped away before
        new_line = new_line.rstrip()    # clean empty spaces on the right, no one needs that...
        if new_line[-1] != '\n': new_line = new_line+'\n'
        if FILE_IS=='SUBMIT': new_line = new_line.replace('#@', '#@ ').replace('=', ' = ')    # optimizing format of submit script

        if DEBUG:
            print('~ DEBUG: old line: '+str(data_raw[line_matches[0]])[:-1])
            print('~ DEBUG: new line: '+str(new_line)[:-1])

        if not DEBUG:
            ####### do backup of file before changing it
            from shutil import copyfile
            target = join(sim.path, 'pc/backups/'+timestamp())
            mkdir(target); target = join(target, filename)
            copyfile(absolute_filepath, target)

            # replace line in raw data
            data_raw[line_matches[0]] = new_line

            # save on drive
            f.close()
            with open(absolute_filepath, 'w') as f:
                for l in data_raw: f.write(l)

    ######## DONE!
    return q
