# param.py
#
# Read the parameters for the simulation.
# Requires: nl2python perl script (based on Wolfgang Dobler's nl2idl script).
#
# Authors:
# J. Oishi (joishi@amnh.org).
# S. Candelaresi (iomsn1@gmail.com)
"""
Contains the parameters of the simulation.
"""


def param(*args, **kwargs):
    """
    Read Pencil Code simulation parameters.
    Requires: nl2python perl script (based on Wolfgang Dobler's nl2idl script).

    call signature:

    read(datadir='data', param2=False, quiet=True, asdict=False, nest_dict=False)

    Keyword arguments:

    *datadir*:
      Directory where the data is stored.

    *param2*:
      Selects the set of parameters.

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


    def read(self, datadir='data', param2=False, quiet=True,
             asdict=False, nest_dict=False):
        """
        Read Pencil Code simulation parameters.
        Requires: nl2python perl script (based on Wolfgang Dobler's nl2idl script).

        call signature:

        read(datadir='data', param2=False, quiet=True, asdict=False, nest_dict=False)

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *param2*:
          Selects the set of parameters.

        *quiet*
          Flag for switching of output.

        *asdict*
          Reads parameters as dictionary.

        *nest_dict*
          Reads parameters as nested dictionary.
        """

        import os

        datadir = os.path.expanduser(datadir)

        if param2:
            filen = os.path.join(datadir, 'param2.nml')
        else:
            filen = os.path.join(datadir, 'param.nml')

        # Execute output of nl2python script.
        if not os.path.exists(filen):
            print("read.param: no such file {0}.".format(filen))
            raise ValueError

        # Read the parameters into a dictionary.
        if asdict:
            if nest_dict:
                param_list = self.__read_nml(filen, nest=True)
            else:
                param_list = self.__read_nml(filen)
            if not quiet:
                print(param_list)
            key_list = param_list.keys()
            for key in key_list:
                setattr(self, key, param_list[key])
        # Read the parameters as attributes to class Params.
        else:
            cmd = 'nl2python ' + filen
            script = os.popen(cmd).read()
            if not quiet:
                print(script)
            if script:
                # This import is needed to execute the script.
                import numpy
                # TODO: Do this without calling a shell command.
                exec(script.replace("\n    ", "\nself.")[198:])
                del(numpy)
            else:
                print("Param.read: nl2python returned nothing! Is $PENCIL_HOME/bin in the path?")
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


    def __read_nml(self, file_name, nest=False):
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

        params = dict()
        for rawline in open(file_name):
            line = rawline.rstrip('\n')
            if line == ' ':
                continue
            print(line)
            if line[0] == "&":
                super_name = line[1:].lower()
                if nest:
                    params[super_name] = dict()
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
                        params[super_name][name] = value
                    else:
                        params[name] = value
        return params
