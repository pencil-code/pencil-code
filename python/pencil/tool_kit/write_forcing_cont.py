import numpy as np


def write_forcing_cont(a, outfile="forcing_cont.dat"):
    """
    Writes the file forcing_cont.dat that can be used to specify the form of the continuous forcing in forcing.f90

    Parameters
    ----------
    a : numpy array specifying the continuous forcing. Shape is expected to be (3,nz,ny,nx). Note that the order of the spatial indices is the same as in the rest of the Pencil Python module.
    
    outfile : file into which the array should be written

    Example usage
    -------------
    >>> import pencil as pc
    >>> import numpy as np
    >>> dim = pc.read.dim()
    >>> a = np.ones((3,dim.nx,dim.ny,dim.nz))
    >>> pc.tool_kit.write_forcing_cont(a)
    """
    a = np.transpose(a, axes=[0,3,2,1]) #flip order of spatial axes to what the Fortran code uses
    a_ = np.reshape(a, np.shape(a), order="F")

    with open(outfile, "w") as f:
        """
        Documention for list-directed IO: https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc5/index.html

        Note that tabs apparently cause problems in some Fortran compilers.
        """
        for elem, i in zip(np.nditer(a_, order="F"), range(0, np.size(a_))):
            f.write("{}\n".format(elem))
