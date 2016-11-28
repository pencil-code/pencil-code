# averages.py
#
# Read the average files.
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Contains the classes and methods to read average files.
"""


def aver(*args, **kwargs):
    """
    Read Pencil Code average data.

    call signature:

    read(self, plane_list=['xy', 'xz', 'yz'], data_dir='data'):

    Keyword arguments:

    *plane_list*:
      A list of the planes over which the averages were taken.

    *data_dir*:
      Directory where the data is stored.
    """

    averages_tmp = Averages()
    averages_tmp.read(*args, **kwargs)
    return averages_tmp


class Averages(object):
    """
    Averages -- holds Pencil Code averages data and methods.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        import numpy as np

        self.t = np.array([])


    def read(self, plane_list=['xy', 'xz', 'yz'], data_dir='data'):
        """
        Read Pencil Code average data.

        call signature:

        read(self, plane_list=['xy', 'xz', 'yz'], data_dir='data'):

        Keyword arguments:

        *plane_list*:
          A list of the planes over which the averages were taken.

        *data_dir*:
          Directory where the data is stored.
        """

        from pencilnew import read
        import numpy as np
        import os

        # Initialize the planes list.
        if plane_list:
            if isinstance(plane_list, list):
                plane_list = plane_list
            else:
                plane_list = [plane_list]
        else:
            plane_list = ['xy', 'xz', 'yz']

        # Determine which average files to read.
        in_file_name_list = []
        aver_file_name_list = []
        if plane_list.count('xy') > 0:
            in_file_name_list.append('xyaver.in')
            aver_file_name_list.append('xyaverages.dat')
        if plane_list.count('xz') > 0:
            in_file_name_list.append('xzaver.in')
            aver_file_name_list.append('xzaverages.dat')
        if plane_list.count('yz') > 0:
            in_file_name_list.append('yzaver.in')
            aver_file_name_list.append('yzaverages.dat')
        else:
            print("error: invalid plane name")
            return -1

        class Foo():
            pass

        for plane, in_file_name, aver_file_name in zip(plane_list, in_file_name_list, aver_file_name_list):
            # This one will store the data.
            ext_object = Foo()

            # Get the averaged quantities.
            file_id = open(os.path.join(os.path.dirname(data_dir), in_file_name))
            variables = file_id.readlines()
            file_id.close()
            for i in range(sum(list(map(self._equal_newline, variables)))):
                variables.remove('\n')
            n_vars = len(variables)
    
            # Determine the structure of the xy/xz/yz averages.
            nw_string = 'n' + 'xyz'[np.where(np.array(map(plane.find, 'xyz')) == -1)[0][0]]
            nw = getattr(read.dim(), nw_string)
            file_id = open(os.path.join(data_dir, aver_file_name))
            aver_lines = file_id.readlines()
            file_id.close()
            entry_length = int(nw*n_vars/8)
            n_times = len(aver_lines)/(1 + entry_length)
            
            # Prepare the data arrays.
            t = np.zeros(n_times, dtype=np.float32)
            data_raw = np.zeros([n_times, n_vars*nw])
            
            # Read the data
            line_idx = 0
            t_idx = -1
            for current_line in aver_lines:
                if line_idx % (entry_length+1) == 0:
                    t_idx += 1
                    t[t_idx] = np.float32(current_line)
                    raw_idx = 0
                else:
                    data_raw[t_idx, raw_idx*8:(raw_idx*8+8)] = list(map(np.float32, current_line.split()))
                    raw_idx += 1
                line_idx += 1
            
            # Restructure the raw data and add it to the Averages object.
            data_raw = np.reshape(data_raw, [n_times, n_vars, nw])
            var_idx = 0
            for var in variables:
                setattr(ext_object, var.strip(), data_raw[:, var_idx, :])
                var_idx += 1
        
            self.t = t
            setattr(self, plane, ext_object)
    
        del(data_raw)
        del(ext_object)

    
    def _equal_newline(self, line):
        """
        Determine if string is equal new line.
        """
        
        if line == '\n':
            return True
        else:
            return False
