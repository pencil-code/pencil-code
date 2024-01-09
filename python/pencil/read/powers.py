# power.py
#
# Read the power spectra of the simulation.
#
"""
Contains the classes and methods to read the power spectra.
"""

import os
import numpy as np
from pencil import read
from pencil.util import ffloat
import re

def power(*args, **kwargs):
    """
    power(datadir='data', file_name='', quiet=False)

    Read the power spectra.

    Parameters
    ----------
    datadir : string
        Directory where the data is stored.

    file_name : string
        Filename to read.
        If a filename is given, only that power spectrum is read.
        By default it reads all the power spectrum files.

    quiet : bool
        Flag for switching off output.

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

    def read(self, datadir="data", file_name="", quiet=False):
        """
        read(datadir='data', file_name='', quiet=False)
    
        Read the power spectra.
    
        Parameters
        ----------
        datadir : string
            Directory where the data is stored.
    
        file_name : string
            Filename to read.
            If a filename is given, only that power spectrum is read.
            By default it reads all the power spectrum files.
    
        quiet : bool
            Flag for switching off output.
    
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
        import os.path as op

        power_list = []
        file_list = []

        if file_name:
            print("Reading only ", file_name)
            try:
                if op.isfile(op.join(datadir, file_name)):
                    # print("read one file")
                    if file_name[:5] == "power" and file_name[-4:] == ".dat":
                        if file_name[:6] == "power_":
                            power_list.append(file_name.split(".")[0][6:])
                            print("appending", file_name.split(".")[0][6:])
                        else:
                            power_list.append(file_name.split(".")[0][5:])
                            print("appending", file_name.split(".")[0][5:])
                        file_list.append(file_name)
                else:
                    print("File does not exist, exiting")
            except IOError:
                print("File does not exist, exiting")
                return

        else:

            # Find the existing power files.

            # power_list = []
            # file_list = []
            for file_name in os.listdir(datadir):
                if file_name[:5] == "power" and file_name[-4:] == ".dat":
                    if file_name[:6] == "power_":
                        power_list.append(file_name.split(".")[0][6:])
                    else:
                        power_list.append(file_name.split(".")[0][5:])
                    file_list.append(file_name)

        # Read the power spectra.
        for power_name, file_name in zip(power_list, file_list):
            if not quiet:
                print(file_name)

            if (
                file_name == "powero.dat"
                or file_name == "powerb.dat"
                or file_name == "powera.dat"
                ):
                # Exclude some incompatible files.
                pass
            elif re.match("power.*_xy.dat", file_name):
                self._read_power2d(power_name, file_name, datadir)
            elif (
                file_name == "poweruz_x.dat"
                or file_name == "powerux_x.dat"
                or file_name == "poweruy_x.dat"
                ):
                self._read_power_1d(power_name, file_name, datadir)
            elif file_name == "power_krms.dat":
                self._read_power_krms(power_name, file_name, datadir)
            else:
                self._read_power(power_name, file_name, datadir)

    def _read_power2d(self, power_name, file_name, datadir):
        """
        Handles output of power_xy subroutine.
        """
        dim = read.dim(datadir=datadir)
        param = read.param(datadir=datadir)

        infile = open(os.path.join(datadir, file_name), "r")
        line_list = infile.readlines()
        infile.close()

        # Get k vectors:
        if param.lintegrate_shell:
            nk = int(
                line_list[1]
                .split()[line_list[1].split().index("k") + 1]
                .split(")")[0][1:]
                )
            ini = 2
            k = []
            for i in range(ini, int(np.ceil(nk / 8)) + ini):
                k.extend([float(j) for j in line_list[i].split()])
            k = np.array(k)
            self.k = k
            ini = i + 1
        else:
            nkx = int(
                line_list[1]
                .split()[line_list[1].split().index("k_x") + 1]
                .split(")")[0][1:]
                )
            ini = 2
            kx = []
            for i in range(ini, int(np.ceil(nkx / 8)) + ini):
                kx.extend([float(j) for j in line_list[i].split()])
            kx = np.array(kx)
            self.kx= kx
            ini = i + 1

            nky = int(
                line_list[1]
                .split()[line_list[1].split().index("k_y") + 1]
                .split(")")[0][1:]
                )
            ky = []
            for i in range(ini, int(np.ceil(nky / 8)) + ini):
                ky.extend([float(j) for j in line_list[i].split()])
            ky = np.array(ky)
            self.ky = ky
            ini = i + 1

            nk = nkx * nky

        # Now read z-positions, if any
        if param.lintegrate_z:
            nzpos = 1
        else:
            if "z-pos" in line_list[ini]:
                print("More than 1 z-pos")
                nzpos = int(re.search(r"\((\d+)\)", line_list[ini])[1])
                ini += 1
                zpos = np.array([float(j) for j in line_list[ini].split()])
                ini += 1
            else:
                nzpos = dim.nzgrid
                grid = read.grid(datadir=datadir, trim=True, quiet=True)
                zpos = grid.z
            self.zpos = zpos
        self.nzpos = nzpos

        # Now read the rest of the file
        line_list = line_list[ini:]
        time = []
        power_array = []
        linelen = len(line_list[1].strip().split())

        if param.lintegrate_shell:
            block_size = np.ceil(nk / 8) * nzpos + 1
        else:
            block_size = np.ceil(int(nk * nzpos) / 8) + 1

        for line_idx, line in enumerate(line_list):
            if np.mod(line_idx, block_size) == 0:
                time.append(float(line.strip()))
            else:
                if linelen == 8:
                    # real power spectrum
                    for value_string in line.strip().split():
                        power_array.append(ffloat(value_string))

                elif linelen == 16:
                    # complex power spectrum
                    real = line.strip().split()[0::2]
                    imag = line.strip().split()[1::2]
                    for a, b in zip(real, imag):
                        power_array.append(ffloat(a) + 1j * ffloat(b))

        time = np.array(time)

        if linelen == 8:
            power_array = np.array(power_array, dtype=np.float32)
        elif linelen == 16:
            power_array = np.array(power_array, dtype=complex)

        if param.lintegrate_shell or (dim.nxgrid == 1 or dim.nygrid == 1):
            power_array = power_array.reshape([len(time), nzpos, nk])
        else:
            power_array = power_array.reshape([len(time), nzpos, nky, nkx])

        self.t = time.astype(np.float32)
        setattr(self, power_name, power_array)

    def _read_power_1d(self, power_name, file_name, datadir):
        """
        Handle output of subroutine power_1d
        """
        dim = read.dim(datadir=datadir)

        block_size = np.ceil(int(dim.nxgrid / 2) / 8.0) + 1

        time = []
        power_array = []
        with open(os.path.join(datadir, file_name), "r") as f:
            for line_idx, line in enumerate(f):
                if np.mod(line_idx, block_size) == 0:
                    time.append(float(line.strip()))
                elif line.find(",") == -1:
                    # if the line does not contain ',', assume it represents a series of real numbers.
                    for value_string in line.strip().split():
                        power_array.append(float(value_string))
                else:
                    # Assume we have complex numbers.
                    for value_string in line.strip().split("( ")[1:]:
                        value_string = (
                            value_string.replace(")", "j")
                            .strip()
                            .replace(", ", "")
                            .replace(" ", "+")
                        )
                        power_array.append(complex(value_string))

        time = np.array(time)
        power_array = np.array(power_array).reshape([len(time), int(dim.nxgrid / 2)])
        self.t = time
        setattr(self, power_name, power_array)

    def _read_power_krms(self, power_name, file_name, datadir):
        """
        Read power_krms.dat.
        """
        dim = read.dim(datadir=datadir)

        block_size = np.ceil(int(dim.nxgrid / 2) / 8.0) + 1

        power_array = []
        with open(os.path.join(datadir, file_name), "r") as f:
            for line_idx, line in enumerate(f):
                # KG: Is this file expected to contain extra lines? If not, the if below can be removed.
                if line_idx < block_size - 1:
                    for value_string in line.strip().split():
                        power_array.append(float(value_string))
        power_array = (
            np.array(power_array).reshape([int(dim.nxgrid / 2)]).astype(np.float32)
        )
        setattr(self, power_name, power_array)

    def _read_power(self, power_name, file_name, datadir):
        """
        Handles output of power subroutine.
        """
        dim = read.dim(datadir=datadir)

        block_size = np.ceil(int(dim.nxgrid / 2) / 8.0) + 1

        time = []
        power_array = []
        with open(os.path.join(datadir, file_name), "r") as f:
            for line_idx, line in enumerate(f):
                if np.mod(line_idx, block_size) == 0:
                    time.append(float(line.strip()))
                else:
                    for value_string in line.strip().split():
                        power_array.append(ffloat(value_string))

        # Reformat into arrays.
        time = np.array(time)
        power_array = (
            np.array(power_array)
            .reshape([len(time), int(dim.nxgrid / 2)])
            .astype(np.float32)
        )
        self.t = time.astype(np.float32)
        setattr(self, power_name, power_array)
