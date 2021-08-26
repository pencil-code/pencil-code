# power.py
#
# Read the power spectra of the simulation.
#
"""
Contains the classes and methods to read the power spectra.
"""


def power(*args, **kwargs):
    """
    Read the power spectra.

    Signature:

    power(datadir='data', fn='', quiet=False)

    Parameters
    ----------
    *datadir*:  Directory where the data is stored.

    *fn*:  Filename to read.
          If a filename is given, only that power spectrum is read.
          By default it reads all the power spectrum files.

    *quiet*:    Flag for switching off output.

    Returns
    -------
    Class containing the different power spectrum as attributes.

    Notes
    -----
    Use the attribute keys to get a list of attributes

    Examples
    --------
    >>> pw = pc.read.power()
    >>> pw.keys()
    t
    kin
    krms
    hel_kin
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

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(self, datadir="data", fn="", quiet=False):
        """
        Read the power spectra.

        call signature:

        power(datadir='data', fn='', quiet=False)

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *fn*:
          Filename to read.
          If a filename is given, only that power spectrum is read.
          By default it reads all the power spectrum files.

        *quiet*
          Flag for switching off output.
        """

        import os
        import os.path as op
        import numpy as np
        from pencil import read
        from pencil.util import ffloat
        import sys
        import matplotlib as plt

        power_list = []
        file_list = []

        if fn:
            print("Reading only ", fn)
            try:
                if op.isfile(op.join(datadir, fn)):
                    # print("read one file")
                    if fn[:5] == "power" and fn[-4:] == ".dat":
                        if fn[:6] == "power_":
                            power_list.append(fn.split(".")[0][6:])
                            print("appending", fn.split(".")[0][6:])
                        else:
                            power_list.append(fn.split(".")[0][5:])
                            print("appending", fn.split(".")[0][5:])
                        file_list.append(fn)
                else:
                    print("File does not exist, exiting")
            except IOError:
                print("File does not exist, exiting")
                return

        else:

            # Find the existing power files.

            # power_list = []
            # file_list = []
            for fn in os.listdir(datadir):
                if fn[:5] == "power" and fn[-4:] == ".dat":
                    if fn[:6] == "power_":
                        power_list.append(fn.split(".")[0][6:])
                    else:
                        power_list.append(fn.split(".")[0][5:])
                    file_list.append(fn)

        # Determine the file and data structure.
        dim = read.dim(datadir=datadir)
        block_size = np.ceil(int(dim.nxgrid / 2) / 8.0) + 1

        # Read the power spectra.
        for power_idx, fn in enumerate(file_list):
            # Read the raw file.
            infile = open(os.path.join(datadir, fn), "r")
            line_list = infile.readlines()
            infile.close()

            # Extract the numbers from the file strings.
            n_blocks = int(len(line_list) / block_size)

            if not quiet:
                print(fn)

            # For the moment, exclude some incompatible files.
            # if fn == 'powero.dat' or fn == 'poweru.dat' or \
            if fn == "powero.dat" or fn == "powerb.dat" or fn == "powera.dat":
                continue
            elif (
                fn == "powerux_xy.dat"
                or fn == "poweruy_xy.dat"
                or fn == "poweruz_xy.dat"
            ):
                # This file has a different number of k

                # This files has the k vector, and irrational numbers
                # Get k vectors:
                nk = 0
                if "k_x" in line_list[1]:
                    nkx = int(
                        line_list[1]
                        .split()[line_list[1].split().index("k_x") + 1]
                        .split(")")[0][1:]
                    )
                    ini = 2
                    kx = []
                    for i in range(ini, int(np.ceil(nkx / 8)) + ini):
                        kx.append([float(j) for j in line_list[i].split()])
                    kx = np.array(list(plt.cbook.flatten(kx)))
                    setattr(self, "kx", kx)
                    ini = i + 1
                    nk = max(nk, nkx)

                if "k_y" in line_list[1]:
                    nky = int(
                        line_list[1]
                        .split()[line_list[1].split().index("k_y") + 1]
                        .split(")")[0][1:]
                    )
                    ky = []
                    for i in range(ini, int(np.ceil(nky / 8)) + ini):
                        ky.append([float(j) for j in line_list[i].split()])
                    ky = np.array(list(plt.cbook.flatten(ky)))
                    setattr(self, "ky", ky)
                    ini = i + 1
                    nk = max(nk, nky)

                if "k_z" in line_list[1]:
                    nkz = int(
                        line_list[1]
                        .split()[line_list[1].split().index("k_z") + 1]
                        .split(")")[0][1:]
                    )
                    kz = []
                    for i in range(ini, int(np.ceil(nkz / 8)) + ini):
                        kz.append([float(j) for j in line_list[i].split()])
                    kz = np.array(list(plt.cbook.flatten(ky)))
                    setattr(self, "kz", kz)
                    ini = i + 1
                    nk = max(nk, nkz)
                # Now read the rest of the file
                # print('ini', ini)
                line_list = line_list[ini:]
                if line_list[0].strip() == "-Infinity":
                    line_list = line_list[1:]
                if line_list[0][0] == "z":
                    line_list = line_list[2:]
                time = []
                power_array = []
                # print('nk', nk)
                block_size = np.ceil(int(nk * 2) / 16.0) + 1
                n_blocks = int(len(line_list) / block_size)
                for line_idx, line in enumerate(line_list):
                    if np.mod(line_idx, block_size) == 0:
                        # print(float(line.strip()))
                        time.append(float(line.strip()))
                        # print("line_idx", line_idx)
                    else:
                        maxi = len(line.strip().split())
                        for j in range(0, maxi, 2):
                            a = line.strip().split()[j]

                            b = line.strip().split()[j + 1]

                            power_array.append(complex(real=ffloat(a), imag=ffloat(b)))
                time = np.array(time)
                power_array = (
                    np.array(power_array)
                    .reshape([n_blocks, int(nk)])
                    .astype(np.complex)
                )

                self.t = time.astype(np.float32)
                setattr(self, power_list[power_idx], power_array)

            elif (
                fn == "poweruz_x.dat" or fn == "powerux_x.dat" or fn == "poweruy_x.dat"
            ):
                # this has irrational numbers

                time = []
                # print('complex reading of file ', fn)
                power_array = []
                for line_idx, line in enumerate(line_list):
                    if np.mod(line_idx, block_size) == 0:
                        # print(float(line.strip()))
                        time.append(float(line.strip()))
                    else:
                        for value_string in line.strip().split("( ")[1:]:
                            value_string = (
                                value_string.replace(")", "j")
                                .strip()
                                .replace(", ", "")
                                .replace(" ", "+")
                            )
                            power_array.append(complex(value_string))

                time = np.array(time)
                power_array = (
                    np.array(power_array)
                    .reshape([n_blocks, int(dim.nxgrid / 2)])
                    .astype(np.complex)
                )
                self.t = time.astype(np.float32)
                setattr(self, power_list[power_idx], power_array)

            elif fn == "power_krms.dat":
                power_array = []
                for line_idx, line in enumerate(line_list):
                    if line_idx < block_size - 1:
                        for value_string in line.strip().split():
                            power_array.append(float(value_string))
                power_array = (
                    np.array(power_array)
                    .reshape([int(dim.nxgrid / 2)])
                    .astype(np.float32)
                )
                setattr(self, power_list[power_idx], power_array)
            else:
                time = []
                power_array = []
                for line_idx, line in enumerate(line_list):
                    if np.mod(line_idx, block_size) == 0:
                        time.append(float(line.strip()))
                    else:
                        for value_string in line.strip().split():
                            power_array.append(ffloat(value_string))

                # Reformat into arrays.
                time = np.array(time)
                power_array = (
                    np.array(power_array)
                    .reshape([n_blocks, int(dim.nxgrid / 2)])
                    .astype(np.float32)
                )
                self.t = time.astype(np.float32)
                setattr(self, power_list[power_idx], power_array)
