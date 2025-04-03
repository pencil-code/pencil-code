import scipy
from scipy.io import FortranFile
import numpy as np

class FortranFileExt(FortranFile):

    def read_record(self, *dtypes, **kwargs):
        """
        Reads a record of a given type from the file.

        Parameters
        ----------
        *dtypes : dtypes, optional
            Data type(s) specifying the size and endianness of the data.

        Returns
        -------
        data : ndarray
            A 1-D array object.

        Raises
        ------
        FortranEOFError
            To signal that no further records are available
        FortranFormattingError
            To signal that the end of the file was encountered
            part-way through a record

        Notes
        -----
        If the record contains a multidimensional array, you can specify
        the size in the dtype. For example::

            INTEGER var(5,4)

        can be read with::

            read_record('(4,5)i4').T

        Note that this function does **not** assume the file data is in Fortran
        column major order, so you need to (i) swap the order of dimensions
        when reading and (ii) transpose the resulting array.

        Alternatively, you can read the data as a 1-D array and handle the
        ordering yourself. For example::

            read_record('i4').reshape(5, 4, order='F')

        For records that contain several variables or mixed types (as opposed
        to single scalar or array types), give them as separate arguments::

            double precision :: a
            integer :: b
            write(1) a, b

            record = f.read_record('<f4', '<i4')
            a = record[0]  # first number
            b = record[1]  # second number

        and if any of the variables are arrays, the shape can be specified as
        the third item in the relevant dtype::

            double precision :: a
            integer :: b(3,4)
            write(1) a, b

            record = f.read_record('<f4', np.dtype(('<i4', (4, 3))))
            a = record[0]
            b = record[1].T

        NumPy also supports a short syntax for this kind of type::

            record = f.read_record('<f4', '(3,3)<i4')

        See Also
        --------
        read_reals
        read_ints

        """
        dtype = kwargs.pop('dtype', None)
        if kwargs:
            raise ValueError(f"Unknown keyword arguments {tuple(kwargs.keys())}")

        if dtype is not None:
            dtypes = dtypes + (dtype,)
        elif not dtypes:
            raise ValueError('Must specify at least one dtype')

        dtypes = tuple(np.dtype(dtype) for dtype in dtypes)
        block_size = sum(dtype.itemsize for dtype in dtypes)   # size of one block of data defined by dtypes

        total_size = 0
        rec_sizes = np.array([],dtype=int)
        start_pos = self._fp.tell()
        num_subrecords = 0

        while True:
            first_size = self._read_size(eof_ok=True)
            size = first_size
            if size < 0:
                size = -size

            rec_sizes=np.append(rec_sizes,size)
            total_size += size
            num_subrecords += 1

            maxoffset = 2**31-1
            num_reps, remainder = divmod(size, maxoffset)
            for i in range(num_reps):
                self._fp.seek(maxoffset,1)

            self._fp.seek(remainder,1)
            second_size = abs(self._read_size())
            #print("FIRST-SECOND",first_size,second_size)

            if size != second_size:
                raise ValueError('Sizes do not agree in the header and footer for '
                                 'this record - check header dtype')
            if first_size > 0:
                break

        end_pos = self._fp.tell()

        num_blocks, remainder = divmod(total_size, block_size)  # num_blocks: number of blocks in logical record
        if remainder != 0:
            raise ValueError(f'Size obtained ({total_size}) is not a multiple of the '
                             f'dtypes given ({block_size}).')

        if len(dtypes) != 1 and total_size != block_size:       # if more than one type and num_blocks>1  -  yet unclear, see Fortran array of structures
            # Fortran does not write mixed type array items in interleaved order,
            # and it's not possible to guess the sizes of the arrays that were written.
            # The user must specify the exact sizes of each of the arrays.
            raise ValueError(f'Size obtained ({total_size}) does not match with the '
                             f'expected size ({block_size}) of multi-item record')
        #print('num_subrecords,total_size,rec_sizes,start_pos=', num_subrecords,total_size,rec_sizes,start_pos)

        self._fp.seek(start_pos)                       # reset file pointer to original record head
        if num_subrecords > 1:                         # combine all subrecords in tmparray
            tmparr = bytearray(total_size)
            istart = 0
            for i in range(num_subrecords):
                first_size = self._read_size(eof_ok=True)
                iend = istart+rec_sizes[i]
                r = self._fp.read(rec_sizes[i])
                if len(r) != rec_sizes[i]:
                    raise FortranFormattingError(
                        "End of file in the middle of a record")
                tmparr[istart:iend] = r
                second_size = self._read_size(eof_ok=True)
                istart=iend

        data = []
        for dtype in dtypes:
            if num_subrecords == 1:
                first_size = self._read_size(eof_ok=True)
                r = np.fromfile(self._fp, dtype=dtype, count=num_blocks)
                if len(r) != num_blocks:
                    raise FortranFormattingError(
                        "End of file in the middle of a record")
                if dtype.shape != ():
                    # Squeeze outmost block dimension for array items
                    if num_blocks == 1:
                        assert r.shape == (1,) + dtype.shape
                        r = r[0]

                data.append(r)
            else:
                #print('DTYPE,NUM_BLOCKS=',dtype,num_blocks)
                r = np.frombuffer(tmparr, dtype=dtype, count=num_blocks)
                data.append(r)

        self._fp.seek(end_pos)                       # reset file pointer behind original record end

        # Unpack result
        if len(dtypes) == 1:
            return data[0]
        else:
            return tuple(data)
