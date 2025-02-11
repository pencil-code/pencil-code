# Authors: Matthew Brett, Travis Oliphant

"""
Class for reading and writing numpy arrays from / to binary files
"""

import struct
import sys

import numpy as N

__all__ = ['sys_endian_code', 'npfile']

sys_endian_code = (sys.byteorder == 'little') and '<' or '>'

class npfile(object):
    ''' Class for reading and writing numpy arrays to/from files
    
    Inputs:
      file_name -- The complete path name to the file to open
                   or an open file-like object
      permission -- Open the file with given permissions: ('r', 'w', 'a')
                    for reading, writing, or appending.  This is the same
                    as the mode argument in the builtin open command.
      format -- The byte-ordering of the file:
                (['native', 'n'], ['ieee-le', 'l'], ['ieee-be', 'B']) for
                native, little-endian, or big-endian respectively.

    Attributes:
      endian   -- default endian code for reading / writing
      order    -- default order for reading writing ('C' or 'F')
      file     -- file object containing read / written data

    Methods:
      seek, tell, close  -- as for file objects
      rewind             -- set read position to beginning of file
      read_raw           -- read string data from file (read method of file)
      write_raw          -- write string data to file (write method of file)
      read_array         -- read numpy array from binary file data
      write_array        -- write numpy array contents to binary file
      
    Example use:
    >>> from StringIO import StringIO
    >>> import numpy as N
    >>> from scipy.io import npfile
    >>> arr = N.arange(10).reshape(5,2)
    >>> # Make file-like object (could also be file name)
    >>> my_file = StringIO()
    >>> npf = npfile(my_file)
    >>> npf.write_array(arr)
    >>> npf.rewind()
    >>> npf.read_array((5,2), arr.dtype)
    >>> npf.close()
    >>> # Or read write in Fortran order, Big endian
    >>> # and read back in C, system endian
    >>> my_file = StringIO()
    >>> npf = npfile(my_file, order='F', endian='>')
    >>> npf.write_array(arr)
    >>> npf.rewind()
    >>> npf.read_array((5,2), arr.dtype)
    '''

    def __init__(self, file_name,
                 permission='rb',
                 endian = 'dtype',
                 order = 'C'):
        if 'b' not in permission: permission += 'b'
        if isinstance(file_name, str):
            self.file = open(file_name, permission)
        else:
            try:
                closed = file_name.closed
            except AttributeError:
                raise TypeError('Need filename or file object as input')
            if closed:
                raise TypeError('File object should be open')
            self.file = file_name
        self.endian = endian
        self.order = order

    def get_endian(self):
        return self._endian
    def set_endian(self, endian_code):
        self._endian = self.parse_endian(endian_code)
    endian = property(get_endian, set_endian, None, 'get/set endian code')
                                     
    def parse_endian(self, endian_code):
        ''' Returns valid endian code from wider input options'''
        if endian_code in ['native', 'n', 'N','default', '=']:
            return sys_endian_code
        elif endian_code in ['swapped', 's', 'S']:
            return sys_endian_code == '<' and '>' or '<'
        elif endian_code in ['ieee-le','l','L','little-endian',
                             'little','le','<']:
            return '<'
        elif endian_code in ['ieee-be','B','b','big-endian',
                             'big','be', '>']:
            return '>'
        elif endian_code == 'dtype':
            return 'dtype'
        else:
            raise ValueError("Unrecognized endian code: " + endian_code)
        return

    def __del__(self):
        try:
            self.file.close()
        except:
            pass

    def close(self):
        self.file.close()

    def seek(self, *args):
        self.file.seek(*args)

    def tell(self):
        return self.file.tell()

    def rewind(self,howmany=None):
        """Rewind a file to its beginning or by a specified amount.
        """
        if howmany is None:
            self.seek(0)
        else:
            self.seek(-howmany,1)

    def read_raw(self, size=-1):
        """Read raw bytes from file as string."""
        return self.file.read(size)

    def write_raw(self, str):
        """Write string to file as raw bytes."""
        return self.file.write(str)

    def remaining_bytes(self):
        cur_pos = self.tell()
        self.seek(0, 2)
        end_pos = self.tell()
        self.seek(cur_pos)
        return end_pos - cur_pos

    def _endian_order(self, endian, order):
        ''' Housekeeping function to return endian, order from input args '''
        if endian is None:
            endian = self.endian
        else:
            endian = self.parse_endian(endian)
        if order is None:
            order = self.order
        return endian, order

    def _endian_from_dtype(self, dt):
        dt_endian = dt.byteorder
        if dt_endian == '=':
            dt_endian = sys_endian_code
        return dt_endian
    
    def write_array(self, data, endian=None, order=None):
        ''' Write to open file object the flattened numpy array data

        Inputs
        data      - numpy array or object convertable to array
        endian    - endianness of written data
                    (can be None, 'dtype', '<', '>')
                    (if None, get from self.endian)
        order     - order of array to write (C, F)
                    (if None from self.order)
        '''
        endian, order = self._endian_order(endian, order)
        data = N.asarray(data)
        dt_endian = self._endian_from_dtype(data.dtype)
        if not endian == 'dtype':
            if dt_endian != endian:
                data = data.byteswap()
        self.file.write(data.tostring(order=order))
        
    def read_array(self, dt, shape=-1, endian=None, order=None):
        '''Read data from file and return it in a numpy array.
        
        Inputs
        ------
        dt        - dtype of array to be read
        shape     - shape of output array, or number of elements
                    (-1 as number of elements or element in shape
                    means unknown dimension as in reshape; size
                    of array calculated from remaining bytes in file)
        endian    - endianness of data in file
                    (can be None, 'dtype', '<', '>')
                    (if None, get from self.endian)
        order     - order of array in file (C, F)
                    (if None get from self.order)

        Outputs
        arr       - array from file with given dtype (dt)
        '''
        endian, order = self._endian_order(endian, order)
        dt = N.dtype(dt)
        try:
            shape = list(shape)
        except TypeError:
            shape = [shape]
        minus_ones = shape.count(-1)
        if minus_ones == 0:
            pass
        elif minus_ones == 1:
            known_dimensions_size = -N.prod(shape,axis=0) * dt.itemsize
            unknown_dimension_size, illegal = divmod(self.remaining_bytes(),
                                                     known_dimensions_size)
            if illegal:
                raise ValueError("unknown dimension doesn't match filesize")
            shape[shape.index(-1)] = unknown_dimension_size
        else:
            raise ValueError(
                "illegal -1 count; can only specify one unknown dimension")
        sz = dt.itemsize * N.prod(shape)
        dt_endian = self._endian_from_dtype(dt)
        buf = self.file.read(sz)
        arr = N.ndarray(shape=shape,
                         dtype=dt,
                         buffer=buf,
                         order=order)
        if (not endian == 'dtype') and (dt_endian != endian):
            return arr.byteswap()
        return arr.copy()

    def fort_write(self,data,endian=None,order=None,head_size=4):
        """Write a Fortran binary record from a numpy array

        Inputs:

          fmt -- If a string then it represents the same format string as
                 used by struct.pack.  The remaining arguments are passed
                 to struct.pack.

                 If fmt is an array, then this array will be written as
                 a Fortran record using the output type args[0].

          *args -- Arguments representing data to write.
        """
        endian, order = self._endian_order(endian, order)
        if endian == '<':
            nfmt = "<"
        elif endian == '>':
            nfmt = ">"
        else:
            nfmt = ""
        if head_size == 4:
            nfmt+= 'i'
        elif head_size == 8:
            nfmt+='L'
        else:
            raise TypeError("Unknown head_size. Valid vaules are 4 & 8.")

        #outstr = struct.pack(data.dtype,data.tostring(order=order))
        outstr = data.tostring(order=order)
        strlen = struct.pack(nfmt,len(outstr))
        self.file.write(strlen)
        self.file.write(outstr)
        self.file.write(strlen)

    def fort_read(self, dt, shape=-1, endian=None, order=None, head_size=4):
        '''Read data from a fortran binary record and return it in a numpy array.

        note that fortran records give a 4-byte (or 8-byte if you use
        gfortran) header describing the number of bytes in a
        record. if shape does not agree with this, an error is
        raised. 

        Inputs
        ------
        dt        - dtype of array to be read
        shape     - shape of output array, or number of elements
                    (-1 as number of elements or element in shape
                    means unknown dimension as in reshape; size
                    of array calculated from remaining bytes in file)
        endian    - endianness of data in file
                    (can be None, 'dtype', '<', '>')
                    (if None, get from self.endian)
        order     - order of array in file (C, F)
                    (if None get from self.order)
        head_size - the size in bytes of the fortran record header.
                    (almost always 4, except gfortran uses 8 by default)

        Outputs
        arr       - array from file with given dtype (dt)
        '''

        endian, order = self._endian_order(endian, order)
        dt = N.dtype(dt)
        dt_endian = self._endian_from_dtype(dt)
        #first, read header
        buf = self.file.read(head_size)
        header = N.ndarray(shape=(1),dtype='i',buffer=buf)[0]
        if (not endian == 'dtype') and (dt_endian != endian):
            header = header.byteswap()
        try:
            shape = list(shape)
        except TypeError:
            shape = [shape]
        minus_ones = shape.count(-1)
        if minus_ones == 0:
            pass
        elif minus_ones == 1:
            known_dimensions_size = -N.prod(shape,axis=0) * dt.itemsize
            unknown_dimension_size, illegal = divmod(header,
                                                     known_dimensions_size)
            if illegal:
                raise ValueError("unknown dimension doesn't match record size")
            shape[shape.index(-1)] = unknown_dimension_size
        else:
            raise ValueError(
                "illegal -1 count; can only specify one unknown dimension")
        sz = dt.itemsize * N.prod(shape)

        buf = self.file.read(sz)
        arr = N.ndarray(shape=shape,
                         dtype=dt,
                         buffer=buf,
                         order=order)
        #fortran record ends with the header repeated. skip this.
        self.seek(head_size,1)
        if (not endian == 'dtype') and (dt_endian != endian):
            return arr.byteswap()
        return arr.copy()
