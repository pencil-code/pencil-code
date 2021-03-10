# powers.py
#
# Read the power spectra of the simulation.
#
# Authors:
# B. Dintrans (boris.dintrans@irap.omp.eu)
# S. Candelaresi (iomsn1@gmail.com)
"""
Contains the classes and methods to read the power spectra.
"""


def power(*args, **kwargs):
    """
    Read the power spectra.

    call signature:

    power(datadir='data', quiet=False)

    Keyword arguments:

    *datadir*:
      Directory where the data is stored.

    *quiet*
      Flag for switching off output.
    """

    power_tmp = Power()
    power_tmp.read(*args, **kwargs)
    return power_tmp


class Power(object):
    """
    Power -- holds power spectrum data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.t = []


    def read(self, datadir='data', quiet=False):
        """
        Read the power spectra.

        call signature:

        power(datadir='data',quiet=False)

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *quiet*
          Flag for switching off output.
        """

        import os
        import numpy as np
        from .. import read

        # Find the exsisting power files.
        power_list = []
        file_list = []
        for file_name in os.listdir(datadir):
            if file_name[:5] == 'power' and file_name[-4:] == '.dat':
                if file_name[:6] == 'power_':
                    power_list.append(file_name.split('.')[0][6:])
                else:
                    power_list.append(file_name.split('.')[0][5:])
                file_list.append(file_name)

        # Determine the file and data structure.
        dim = read.dim(datadir=datadir)
        block_size = np.ceil(int(dim.nxgrid/2)/8.) + 1

        # Read the power spectra.
        for power_idx, file_name in enumerate(file_list):
            # For the moment, exclude some incompatible files.
            if file_name == 'powero.dat' or file_name == 'poweru.dat' or \
                file_name == 'powerb.dat' or file_name == 'powera.dat':
                continue
            if not quiet:
                print(file_name)

            # Read the raw file.
            infile = open(os.path.join(datadir, file_name), 'r')
            line_list = infile.readlines()
            infile.close()

            # Extract the numbers from the file strings.
            n_blocks = int(len(line_list)/block_size)
            if not file_name == 'power_krms.dat':
                time = []
                power_array = []
                for line_idx, line in enumerate(line_list):
                    if np.mod(line_idx, block_size) == 0:
                        time.append(float(line.strip()))
                    else:
                        for value_string in line.strip().split():
                            power_array.append(float(value_string))

                # Reformat into arrays.
                time = np.array(time)
                power_array = np.array(power_array).reshape([n_blocks, int(dim.nxgrid/2)]).astype(np.float32)

                self.t = time.astype(np.float32)
                setattr(self, power_list[power_idx], power_array)
            else:
                power_array = []
                for line_idx, line in enumerate(line_list):
                    if line_idx < block_size-1:
                        for value_string in line.strip().split():
                            power_array.append(float(value_string))
                power_array = np.array(power_array).reshape([int(dim.nxgrid/2)]).astype(np.float32)
                setattr(self, power_list[power_idx], power_array)
