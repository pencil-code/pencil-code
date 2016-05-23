# param.py
#
# Read the parameters for the simulation.
# Requires: nl2python perl script (based on Wolfgang Dobler's nl2idl script).
#
# Author: S. Candelaresi (iomsn1@gmail.com), J. Oishi (joishi@amnh.org).
"""
Contains the parameters of the simulation.
"""


def param(*args, **kwargs):
    """
    Read Pencil Code simulation parameters.
    Requires: nl2python perl script (based on Wolfgang Dobler's nl2idl script).

    call signature:

    read(data_dir='data/', param2=False, quiet=False)

    Keyword arguments:

    *data_dir*:
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


    def read(self, data_dir='data/', param2=False, quiet=False,
             asdict=False, nest_dict=False):
        """
        Read Pencil Code simulation parameters.
        Requires: nl2python perl script (based on Wolfgang Dobler's nl2idl script).

        call signature:

        read(self, data_dir='data/', param2=False, quiet=False)

        Keyword arguments:

        *data_dir*:
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

        data_dir = os.path.expanduser(data_dir)

        if param2:
            filen = os.path.join(data_dir, 'param2.nml')
        else:
            filen = os.path.join(data_dir, 'param.nml')

        # Execute output of nl2python script.
        if not os.path.exists(filen):
            print("Param.read: no such file {0}.".format(filen))
            raise ValueError

        if asdict:
            if nest_dict:
                param_list = self.__read_nml(filen, nest=True)
            else:
                param_list = self.__read_nml(filen)
            if not quiet:
                print(param_list)
        else:
            cmd = 'nl2python '+filen
            script = os.popen(cmd).read()
            if not quiet:
                print(script)
            if script:
                class Params:
                    pass
                exec(script.replace("\n    ", "\nParams.")[198:])
            else:
                print("Param.read: nl2python returned nothing! Is $PENCIL_HOME/bin in the path?")
                return -1
            param_list = Params()

        key_list = dir(param_list)
        for key in key_list:
            setattr(self, key, getattr(Params, key))


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
            else:
                return int(string_part)
        except:
            return re.sub("'", "", string_part)


    def __tuplecatch(self, string):
        """
        Catches name - value tuples in a string.

        call signature:

        __tuplecatch(self, string)

        Keyword arguments:

        *string*:
          String containing the name - value tuple.
        """

        if "(" in string:
            string = string.replace("(", "").replace(")", "").split(",")
            for j in range(len(string)):
                string[j] = self.__param_formatter(string[j])
            return tuple(string)
        else:
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
                            for i in range(int(s[0])):
                                value += [self.__tuplecatch(s[1])]
                        else:
                            value += [self.__tuplecatch(parts[i])]
                    if len(value) == 1:
                        value = value[0]
                    if nest:
                        params[super_name][name] = value
                    else:
                        params[name] = value
        return params
