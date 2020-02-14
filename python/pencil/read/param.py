# param.py
#
# Read the parameters for the simulation.
# Requires: nl2python perl script (based on Wolfgang Dobler's nl2idl script).
# 14/02/20: default now load both init/run pars into Param(object), with
#           asdict and nest_dict also default. Only where parameter is present
#           in more than one module is it nested. Run parameter values replace
#           conflicting init values.
# 
# Authors:
# J. Oishi (joishi@amnh.org).
# S. Candelaresi (iomsn1@gmail.com)
# F. A. Gent (fred.gent.ncl@gmail.com)
"""
Contains the parameters of the simulation.
"""


def param(*args, **kwargs):
    """
    Read Pencil Code simulation parameters.
    Requires: nl2python perl script (based on Wolfgang Dobler's nl2idl script).

    call signature:

    read(datadir='data', param2=False, quiet=True, asdict=True, nest_dict=False)

    Keyword arguments:

    *datadir*:
      Directory where the data is stored.

    *param1*:
      Selects the set of parameters only from start.

    *param2*:
      Selects the set of parameters only from run.

    *quiet*
      Flag for switching of output.

    *asdict*
      Reads parameters as dictionary.

    *nest_dict*
      Reads parameters as nested dictionary.
    """

    param_tmp = Param()
    param_tmp.read(*args, **kwargs)
    return param_tmp


class Param(object):
    """
    Param -- holds the simulation parameters.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.keys = []


    def read(self, datadir='data', param1=False, param2=False, quiet=True,
             asdict=True, nest_dict=True):
        """
        Read Pencil Code simulation parameters.
        Requires: nl2python perl script (based on Wolfgang Dobler's nl2idl script).

        call signature:

        read(datadir='data', param2=False, quiet=True, asdict=True, nest_dict=False)

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *param1*:
          Selects the set of parameters only from start.

        *param2*:
          Selects the set of parameters only from run.

        *quiet*
          Flag for switching of output.

        *asdict*
          Reads parameters as dictionary.

        *nest_dict*
          Reads parameters as nested dictionary.
        """

        import os

        datadir = os.path.expanduser(datadir)

        #Collect files to be read into a list
        files = []
        if param1:
            files.append(os.path.join(datadir, 'param.nml'))
        elif param2:
            files.append(os.path.join(datadir, 'param2.nml'))
        else:
            files.append(os.path.join(datadir, 'param.nml'))
            files.append(os.path.join(datadir, 'param2.nml'))

        # Verify path of files in list.
        for filen in files:
            if not os.path.exists(filen):
                print("read.param: no such file {0}.".format(filen))
                raise ValueError

        # Read the parameters into a dictionary.
        param_list = dict()
        # Construct object from dictionary with Python
        if asdict:
            # Nesting parameters with same name under module attributes
            if nest_dict:
                for filen in files:
                    param_list, param_conflicts, pars_list = \
                                     self.__read_nml(param_list, filen, nest=True)
            # Parameters with same name will be written by last value
            else:
                for filen in files:
                    param_list, param_conflicts, pars_list = \
                                     self.__read_nml(param_list, filen)
            # If nesting occurs report conflicts and record nests to retain
            if len(param_conflicts) > 0:
                keep_list = []
                for key in param_conflicts.keys():
                    print(key,'as',param_conflicts[key][0],'in',
                          param_conflicts[key][1],'conflicts with',
                          param_conflicts[key][2],'in',
                          param_conflicts[key][3])
                    keep_list.append(param_conflicts[key][1])
                    keep_list.append(param_conflicts[key][3])
                # Create object container for nested contents
                class Foo(object):
                    pass
            for par in pars_list:
                if par in keep_list:
                    key_list = []
                    # Strip out nested objects that will not be required
                    for key in param_list[par].keys():
                        if not key in param_conflicts.keys():
                            key_list.append(key)
                    for key in key_list: 
                        param_list[par].__delitem__(key)
            if not quiet:
                print(param_list)
            # Construct class Params object attributes
            key_list = param_list.keys()
            for key in key_list:
                if key in pars_list:
                    # Nest only parameters with name conflicts
                    if key in keep_list:
                        ext_object = Foo()
                        for subkey in param_list[key].keys():
                            print(subkey, 'is nested under', key)
                            setattr(ext_object, subkey, param_list[key][subkey])
                        setattr(self, key, ext_object)
                else:
                    # Unique parameters saved unnested
                    setattr(self, key, param_list[key])
        # Construct object by calling external perl script
        else:
            # Execute output of nl2python script.
            for filen in files:
                cmd = 'nl2python ' + filen
                script = os.popen(cmd).read()
                if not quiet:
                    print(script)
                if script:
                    # This import is needed to execute the script.
                    import numpy
                    # TODO: Do this without calling a shell command.
                    # Done: via asdict
                    exec(script.replace("\n    ", "\nself.")[198:])
                    del(numpy)
                else:
                    print("Param.read: nl2python returned nothing!"+
                          " Is $PENCIL_HOME/bin in the path?")
                    return -1

        return 0


    def __param_formatter(self, string_part):
        """
        Formats the parameters from the files.

        call signature:

        __param_formatter(self, string_part)

        Keyword arguments:

        *string_part*:
          Part of the string to be formatted.
        """

        import re

        string_part = re.sub(" ", "", string_part)
        if string_part == "T":
            return True
        if string_part == "F":
            return False
        try:
            if "." in string_part:
                return float(string_part)
            return int(string_part)
        except:
            return re.sub("'", "", string_part)


    def __tuple_catch(self, string):
        """
        Catches name - value tuples in a string.

        call signature:

        __tuple_catch(self, string)

        Keyword arguments:

        *string*:
          String containing the name - value tuple.
        """

        if "(" in string:
            string = string.replace("(", "").replace(")", "").split(",")
            for j in range(len(string)):
                string[j] = self.__param_formatter(string[j])
            return tuple(string)
        return self.__param_formatter(string)


    def __read_nml(self, params, file_name, nest=False):
        """
        Reads in F90 namelist as dictionary object

        call signature:

        __read_nml(self, file_name, nest=False)

        Keyword arguments:

        *file_name*:
          Name of the file.

        *nest*
          For nested dictionaries.
        """

        import re

        r = re.compile(r'(?:[^,(]|\([^)]*\))+')

        # Contain the nested parameters to be retained
        param_conflicts = dict()
        # Contain the nest names for each parameter set
        super_name_list = []
        for rawline in open(file_name):
            if ' \n' in rawline and not '=' in rawline:
                continue
            line = rawline.rstrip('\n')
            if line == ' ':
                continue
            if " '," in line and not '=' in line:
                continue
            if line[1] == "&" or line[0] == "&":
                super_name = line[2:].lower().rsplit('_pars'
                                            )[0].rsplit('_init'
                                            )[0].rsplit('_run')[0]
                if nest:
                    params[super_name] = dict()
                    super_name_list.append(super_name)
            else:
                line = re.sub("^ ", "", line)
                if line != "/":
                    split = re.split("=", line)
                    name = re.sub(" ", "", split[0].lower())
                    value = split[1]
                    parts = r.findall(value)
                    value = []
                    for i in range(len(parts)):
                        if "*" in parts[i]:
                            s = parts[i].split("*")
                            for j in range(int(s[0])):
                                value += [self.__tuple_catch(s[1])]
                        else:
                            value += [self.__tuple_catch(parts[i])]
                    if len(value) == 1:
                        value = value[0]
                    if nest:
                        # Save all parameters nested and unnested
                        if name not in params.keys():
                            params[name] = value
                            params[super_name][name] = value
                        # If name conflict exists remove unnested copies
                        else:
                            if params[name] == value:
                                if not super_name in ('run','init'):
                                    params[super_name][name] = value
                            else:
                                if super_name in ('run','init'):
                                    params[name] = value
                                else:
                                    if not name in params[super_name].keys():
                                        # Record name conflict details 
                                        for alt_name in super_name_list:
                                            if name in params[alt_name].keys():
                                                param_conflicts[name] = (
                                                        value, super_name,
                                                        params[alt_name][name],
                                                        alt_name        )
                                                params.__delitem__(name)
                                            else:
                                                params[name] = value
                                        params[super_name][name] = value
                                    else:
                                        params[super_name][name] = value
                                        params[name] = value
                    else:
                        params[name] = value
        if len(super_name_list) > 0:
            if 'run' in super_name_list: super_name_list.remove('run')
            if 'init' in super_name_list: super_name_list.remove('init')
        return params, param_conflicts, super_name_list
