# particles_to_vtk.py
#
# Convert the particle data into vtk (geometry) format.
#
# Authors:
# S. Candelaresi (iomsn1@gmail.com).
"""
Contains the class with the vtk data.
"""


def particles_to_vtk(
    var_file="pvar.dat",
    datadir="data",
    proc=-1,
    destination="particles.vtk",
    ti=-1,
    tf=-1,
    binary=True,
):
    """
    Read the pVAR files from Pencil Code and write them as vtk files.

    call signature:

    particles_to_vtk(var_file='pvar.dat', datadir='data', proc=-1,
                     destination='particles.vtk', ti=-1, tf=-1, binary=True)

    Keyword arguments:

    *var_file*:
      Name of the PVAR file.

    *datadir*:
      Directory where the data is stored.

    *proc*:
      Processor to be read. If -1 read all and assemble to one array.

    *destination*:
      Destination file (only if ti, tf > 0).

    *ti*:
      Initial time index for animation.

    *tf*:
      Final time index for animation.

    *binary*:
       Determines if binary or clear text data for the vtk files.
    """

    from pencil import read

    # Read the particle data and save it in a list.
    pvar_list = []
    if (ti >= 0) and (tf >= 0):
        for tidx in range(ti, tf):
            var_file = "PVAR{0}".format(tidx)
            pvar = read.pvar(varfile=var_file, datadir=datadir, proc=proc)
            pvar_list.append(pvar)
    else:
        pvar = read.pvar(varfile=var_file, datadir=datadir, proc=proc)
        pvar_list.append(pvar)

    # Convert the data into vtk.
    particles_vtk_tmp = ParticlesVtk()
    particles_vtk_tmp.convert_to_vtk(pvar_list)

    # Define some constants.
    particles_vtk_tmp.ti = ti
    particles_vtk_tmp.tf = tf
    particles_vtk_tmp.binary = binary
    particles_vtk_tmp.destination = destination

    # Write the data inot vtk files.
    particles_vtk_tmp.write_to_vtk()

    return particles_vtk_tmp


class ParticlesVtk(object):
    """
    ParticlesVtk -- holds the vtk particle data and methods.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.ipars = 0
        self.vtk_grid_data = []
        self.ti = -1
        self.tf = -1
        self.binary = True
        self.destination = "work.vtk"

    def convert_to_vtk(self, pvar_list):
        """
        Convert particle data into vtk format and return the vtk objects.

        call signature:

        convert_to_vtk(self, pvar_list)

        Keyword arguments:

        *pvar_list*:
          List of particle objects.
        """

        import numpy as np
        import vtk

        for pvar in pvar_list:
            points = np.vstack([pvar.xp, pvar.yp, pvar.zp])

            vtk_grid_data = vtk.vtkUnstructuredGrid()
            vtk_points = vtk.vtkPoints()

            # Add the data to the vtk points.
            for point_idx in range(points.shape[1]):
                vtk_points.InsertNextPoint(points[:, point_idx])
            # Add the vtk points to the vtk grid.
            vtk_grid_data.SetPoints(vtk_points)

            self.vtk_grid_data.append(vtk_grid_data)

    def write_to_vtk(self):
        """
        Write the grid data into vtk files.

        call signature:

            write_to_vtk(self)
        """

        import vtk

        if (self.ti >= 0) and (self.tf >= 0):
            for tidx in range(self.ti, self.tf):
                destination = "{0}{1}.vtk".format(self.destination, tidx)
                writer = vtk.vtkUnstructuredGridWriter()
                if self.binary:
                    writer.SetFileTypeToBinary()
                else:
                    writer.SetFileTypeToASCII()

                writer.SetFileName(destination)
                # Ensure compatability between vtk 5 and 6.
                try:
                    writer.SetInputData(self.vtk_grid_data[tidx - self.ti])
                except:
                    writer.SetInput(self.vtk_grid_data[tidx - self.ti])
                writer.Write()
        else:
            writer = vtk.vtkUnstructuredGridWriter()
            if self.binary:
                writer.SetFileTypeToBinary()
            else:
                writer.SetFileTypeToASCII()

            writer.SetFileName(self.destination)
            # Ensure compatability between vtk 5 and 6.
            try:
                writer.SetInputData(self.vtk_grid_data[0])
            except:
                writer.SetInput(self.vtk_grid_data[0])
            writer.Write()
