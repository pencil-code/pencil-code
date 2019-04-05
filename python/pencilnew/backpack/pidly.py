"""pIDLy 0.2.7: IDL within Python.

Control ITT's IDL (Interactive Data Language) from within Python.

https://github.com/anthonyjsmith/pIDLy
http://pypi.python.org/pypi/pIDLy/

Requirements:

* Pexpect
* NumPy

Usage:

>>> import pidly
>>> idl = pidly.IDL()

>>  print(idl.__doc__)

Consult the docstrings or README.txt in the source distribution for
further information.

Copyright (c) 2008-2017, Anthony Smith
anthonysmith80@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

"""
from __future__ import print_function
import sys
import re
import weakref
import atexit
import unittest
import tempfile
import os
import errno
from datetime import datetime

import numpy
import pexpect
from scipy.io import idl as idlio

now = datetime.now
__version__ = '0.2.6'
STR_DELIMITER = '!@#'  # To distinguish items in an array of strings
try:
    __IPYTHON__
    _ipython = True
except NameError:
    _ipython = False


# List of weak references to IDL sessions (primarily for IPython)
try:
    weakrefs_to_pidly_sessions  # Allow IDL to continue when module is reloaded
except NameError:
    weakrefs_to_pidly_sessions = []


def close_all():
    for pidly_weakref in weakrefs_to_pidly_sessions: #[:]:
        if pidly_weakref() is not None:
            pidly_weakref().close()


# Clean exit in IPython
atexit.register(close_all)


class IDLInputOverflowError(Exception):
    """Expression too long for IDL to receive."""


class IDL(pexpect.spawn):
    """pidly.IDL() : Launch IDL session within Python.

    The IDL class inherits from pexpect.spawn.  Consult pexpect
    documentation for details of further methods.

    Usage:

    Initiate:
    >>> import pidly
    >>> idl = pidly.IDL()

    Or:
        idl = pidly.IDL('/path/to/idl')

    Execute commands:
    >>> idl('x = total([1, 1], /int)')

    Retrieve values:
    >>> print(idl.ev('x'))
    2

    Or (slightly slower):
    >>> print(idl.x)
    2

    Evaluate expressions:
    >>> print(idl.ev('x ^ 2'))
    4

    Use cache (IDL save) to handle large arrays:
    >>> idl('x=[1,2,3,4,5,6]')
    >>> print(idl.ev('x', use_cache=True))
    [1 2 3 4 5 6]

    Transfer a list of IDL variables, using cache:
    >>> idl('y=[1,2,3,4,5,6]')
    >>> xy = idl.ev_list(['x','y'], use_cache=True)
    >>> print(sorted(xy.keys()))
    ['x', 'y']
    >>> print(xy['x'])
    [1 2 3 4 5 6]

    Assign value from Python expression:
    >>> idl.x = 2 + 2
    >>> print(idl.x)
    4

    Or:
    >>> idl('x', 2 + 2)
    >>> print(idl.x)
    4

    Perform IDL function on Python expression(s):
    >>> idl.func('reform', range(4), 2, 2)
    array([[0, 1],
           [2, 3]])

    Or (slightly slower):
    >>> idl.reform(range(4), 2, 2)
    array([[0, 1],
           [2, 3]])

    With keywords (/L64 -> L64=True or L64=1)
    >>> x = idl.histogram(range(4), binsize=3, L64=True)
    >>> print(x)
    [3 1]
    >>> print(x.dtype)
    int64

    IDL procedure with Python argument(s):
    >>> idl.pro('plot', range(10), range(10), xstyle=True, ystyle=True)

    Interactive mode:
    >>  idl.interact()
    IDL> print, x
         4
    IDL> ^D
    >>>

    Close:
    >>> idl.close()

    pIDLy supports the transfer of:
    * ints, longs, ...
    * floats, doubles, ...
    * strings
    * arrays of the above types, with arbitrary size and shape
    * dictionaries <-> structures & lists of dicts <-> arrays of structures
      but with certain limitations on transfer from Python to IDL

    [NB if getting Syntax Errors when passing large arrays to IDL, try using
    >>  idl = pidly.IDL(long_delay=0.05)
    default is 0.02.]

    """

    def __init__(self, *arguments, **kwargs):
        # NB set .logfile_send = sys.stdout to monitor commands sent

        # The speed Py -> IDL is limited by the IDL input buffer,
        # which can overload.  self.delaybeforesend is increased
        # before sending long arrays.
        self.short_delay = kwargs.pop('short_delay', 0)
        # As small as possible! (try 0.014)
        self.long_delay = kwargs.pop('long_delay', 0.02)
        # Temp dir for get_from_save
        self._cache_dir = kwargs.pop("cache_dir", None)
        self.use_cache = kwargs.pop("use_cache", False)

        # There are limits to how data may be passed to IDL:
        # max_sendline is the number of bytes that may be sent in one line
        # (Final, additional, byte is the end-of-line)
        self.max_sendline = kwargs.pop('max_sendline', 1023)

        # max_idl_code_area is the maximum number of bytes that may
        # be input as an IDL command. This limit may be reached by
        # splitting the line and sending as more than one send()
        # command
        # (Must be 1244 or greater, otherwise command in ev() fails)
        self.max_idl_code_area = kwargs.pop('max_idl_code_area', 2046) # 2048?

        # Number of array elements in IDL code area limited to 249 (don't ask)
        # (Was 331 on IDL 6.x and 7[?], but 249 on 7.1.1)
        self.max_n_elements_code_area = kwargs.pop('max_n_elements_code_area',
                                                   249)

        # Custom IDL prompt
        self.idl_prompt = kwargs.pop('idl_prompt', 'IDL> ')

        # Begin
        if len(arguments) == 0:
            arguments = ('idl',)
        if 'timeout' not in kwargs:
            kwargs['timeout'] = None
        pexpect.spawn.__init__(self, *arguments, **kwargs)
        self.delaybeforesend = self.short_delay  # Changed for long commands
        self._wait_for_prompt()

        # For clean exit in IPython and for close_all()
        weakrefs_to_pidly_sessions.append(weakref.ref(self))

        self.ready = True  # For __setattr__


    def close(self):
        """
        Close IDL session.

        Try to call IDL exit function - this way you may still be able
        to terminate IDL if it is a called indirectly through a script.
        """
        if self.isalive():
            try:
                self.ex('exit', print_output=False, ret=True)
            except OSError:
                pass
        super(IDL, self).close()


    def ex(self, expression, assignment_value=None,
           print_output=True, ret=False):
        """Execute a command in IDL.

        If assignment_value is set (to a Python expression), this value is
        assigned to the IDL variable named in expression.

        """

        # Assign value to expression?
        if assignment_value is not None:
            expression = self._python_to_idl_input(assignment_value,
                                                   expression)

        # List of commands to execute?
        if hasattr(expression, '__iter__') and not isinstance(expression, (str, bytes, bytearray)):
            # Long assignments are broken down into lists: iterate then return
            # Or can receive a list of commands directly
            output = []
            print("Sending", len(expression), "commands to IDL:", end=' ')
            sys.stdout.flush()
            for exp in expression:
                sys.stdout.write(".")
                sys.stdout.flush()
                out = self.ex(exp, print_output=print_output, ret=ret)
                if out:
                    output.append(out)
                self.delaybeforesend = self.long_delay
            self.delaybeforesend = self.short_delay
            print("done.")
            if ret:
                return ''.join(output)
            else:
                return

        # Send expression to IDL
        if self._send_expression_to_idl(expression):  # Any bytes sent?
            # Wait for command to be completed, and optionally print output
            self.readline()  # First line of output = IDL command repeated
            idl_output = self._wait_for_prompt(print_output=print_output)

            # Return IDL output
            if ret and idl_output:
                return idl_output


    def ev(self, expression, print_output=True, use_cache=None):
        """Return the value of an IDL expression as a numpy.ndarray."""

        # Evaluate expression and store as an IDL variable
        if expression != 'pidly_tmp':
            self.ex('pidly_tmp=' + expression, print_output=print_output)

        if use_cache is None:
            use_cache = self.use_cache
        if use_cache:
            return self._get_from_save(['pidly_tmp'])['pidly_tmp']
        else:
            return self._get_from_print()


    def ev_list(self, names, print_output=True, use_cache=None):
        """Return a dictionary containing values of IDL variables in list names."""
        if not isinstance(names, (list, tuple)):
            raise ValueError("input should be a list or tuple, not ", type(names))
        if use_cache is None:
            use_cache = self.use_cache
        if use_cache:
            return self._get_from_save(names)
        else:
            result = {}
            for eachname in names:
                result[eachname] = self.ev(eachname, print_output=print_output)
            return result


    def interact(self, show_prompt=True, **kwargs):
        """Interactive IDL shell. Press ^D to return to Python."""

        if show_prompt:
            print(self.idl_prompt, end=' ')
        sys.stdout.flush()
        if 'escape_character' not in kwargs:
            kwargs['escape_character'] = '\x04'
        pexpect.spawn.interact(self, **kwargs)
        if not _ipython:
            print("")
    interact.__doc__ += "\n\n        " + pexpect.spawn.interact.__doc__


    def variables(self):
        """Return list of names of defined IDL variables."""

        # Retrieve output from IDL help command ('help' without 'output=...'
        # prints one screen at a time, waiting for spacebar...)
        self.ex('help, output=pidly_tmp')
        help_output = self.ev('pidly_tmp', use_cache=False).tolist()  # String array
        variable_list = []
        for line in help_output[1:]:  # 1st line = "%..."
            if line.startswith('Compiled'):  # End of variable list
                break
            elif line and not line.startswith('%'):
                variable_list.append(line.split()[0])
        return variable_list


    def func(self, name, *args, **kwargs):
        """Evaluate IDL function."""

        try:
            string_bits = []
            for i, arg in enumerate(args):
                string_bits.append(self._python_to_idl_input(arg))
            for j, kwarg in enumerate(kwargs):
                string_bits.append(
                    kwarg + '=' + self._python_to_idl_input(kwargs[kwarg]))
            return self.ev(name + '(' + ','.join(string_bits) + ')', use_cache=False)
        except IDLInputOverflowError:
            string_bits = []
##             arg_vals = args
##             kwarg_vals = kwargs
##             for arg in args:
##                 arg_vals.append(arg)
##             for kwarg in kwargs:
##                 kwarg_vals.append(kwarg)
            for i, arg in enumerate(args):
                self.ex('pidly_tmp' + str(i), arg)
                string_bits.append('pidly_tmp' + str(i))
            for j, kwarg in enumerate(kwargs):
                self.ex('pidly_tmp' + str(len(args) + j), kwargs[kwarg])
                string_bits.append(
                    kwarg + '=' + 'pidly_tmp' + str(len(args) + j))
            return self.ev(name + '(' + ','.join(string_bits) + ')', use_cache=False)


    def pro(self, name, *args, **kwargs):
        """Execute IDL procedure."""
        try:
            string_bits = []
            for i, arg in enumerate(args):
                string_bits.append(self._python_to_idl_input(arg))
            for j, kwarg in enumerate(kwargs):
                string_bits.append(
                    kwarg + '=' + self._python_to_idl_input(kwargs[kwarg]))
            return self.ex(name + ',' + ','.join(string_bits))
        except IDLInputOverflowError:
            string_bits = []
            for i, arg in enumerate(args):
                self.ex('pidly_tmp' + str(i), arg)
                string_bits.append('pidly_tmp' + str(i))
            for j, kwarg in enumerate(kwargs):
                self.ex('pidly_tmp' + str(len(args) + j), kwargs[kwarg])
                string_bits.append(
                    kwarg + '=' + 'pidly_tmp' + str(len(args) + j))
            return self.ex(name + ',' + ','.join(string_bits))


    # Special methods


    # Calling the instance is the same as executing an IDL command.
    __call__ = ex


    def __getattr__(self, name):
        """Get IDL attribute.

        idl.x: return the value of 'x' from IDL.
        idl.f(x,y,...): return the value of IDL f() of Python variables x,y,...

        """

        # idl.x
        if name.upper() in self.variables():  # Takes time!
            return self.ev(name)

        # idl.f(x,y,...)
        elif (not name.startswith('_')
              and name not in ['trait_names', 'pidly_tmp']):  # for IPython
            def idl_function(*args, **kwargs):
                return self.func(name, *args, **kwargs)
            return idl_function


    def __setattr__(self, name, value):
        """Set IDL attribute: idl.x = ..."""

        if 'ready' in self.__dict__:
            # __init__ has finished
            if name in self.__dict__:
                pexpect.spawn.__setattr__(self, name, value)
            else:
                self.ex(name, value)
        else:
            # __init__ in progress
            pexpect.spawn.__setattr__(self, name, value)


    # "PRIVATE" METHODS


    # Sending to IDL


    def _send_expression_to_idl(self, expression):
        """Send a string to IDL and return the no. of bytes sent (or False)."""
        # Only method that sends anything to IDL

        if len(expression) > self.max_sendline:
            if len(expression) <= self.max_idl_code_area:
                # Long line: need to send it in chunks
                expression += '\n'
                for i in range((len(expression) - 1)
                               // (self.max_sendline + 1) + 1):
                    self.send(expression[(self.max_sendline + 1) * i
                                         : (self.max_sendline + 1) * (i + 1)])
                    self.delaybeforesend = self.long_delay
                self.delaybeforesend = self.short_delay
                return True
            else:
                raise IDLInputOverflowError("Expression too long for IDL to receive: cannot execute")
        else:
            if not self.isalive():
                raise OSError("IDL session is not alive.")
            return self.sendline(expression)


    def _python_to_idl_input(self, python_input, assign_to=None):
        """Take Python value and return string suitable for IDL assignment.

        For long input, returns a list of executable strings.
        """

        if assign_to is not None:
            assign_to_str = assign_to + "="
        else:
            assign_to_str = ''

        if isinstance(python_input, str):
            # Strings need additional quotes
            idl_input = assign_to_str + "\'" + python_input + "\'"

        else:
            # Convert to numpy array
            py_in = numpy.array(python_input)

            # Display warning if conversion has changed the array values
            if ((not isinstance(python_input, numpy.ndarray))
                and py_in.tolist() != python_input
                and py_in.dtype.name[0:3] == 'str'):
                print("(!) Conversion to numpy.array has changed input from:", file=sys.stderr)
                print(python_input, file=sys.stderr)
                print("to:", file=sys.stderr)
                print(py_in.tolist(), file=sys.stderr)

            # String format (must have commas between elements)
            if py_in.dtype.name == 'float64':
                str_py_in = ''.join([
                    '[',
                    ','.join(str('%.17e' % s)
                             for s in py_in.flatten().tolist()),
                    ']']).replace(' ', '').replace('e', 'd')
            elif py_in.dtype.name[0:3] == 'str':
                str_py_in = str(py_in.flatten().tolist()).replace("', ", "',")
            else:
                str_py_in = str(py_in.flatten().tolist()).replace(" ", "")
            if len(py_in.shape) > 1:
                # IDL can't handle list concatenations with > 3 dimensions
                str_py_in_shape = ("reform(" + str_py_in + ","
                                   + str(py_in.shape[::-1])[1:-1] + ")")
            elif len(py_in.shape) > 0:
                str_py_in_shape = str_py_in
            else:
                str_py_in_shape = str_py_in[1:-1]  # Remove '[' and ']'

            # Dictionary?  Convert to IDL structure
            if ((not hasattr(py_in.tolist(), 'keys')
                 and hasattr(py_in.tolist(), '__iter__')
                 and hasattr(py_in.tolist()[0], 'keys'))
                or hasattr(py_in.tolist(), 'keys')):
                return self._python_to_idl_structure(python_input, assign_to)
            else:
                # Cast as appropriate type
                idl_input = self._idl_cast_from_dtype(py_in.dtype,
                                                      str_py_in_shape)

            idl_input = ''.join([assign_to_str, idl_input])

            if (len(idl_input) > self.max_idl_code_area
                or len(py_in.flatten()) > self.max_n_elements_code_area):
                # String too long!  Need to create list of shorter commands
                if assign_to is None:
                    raise IDLInputOverflowError("Expression too long for IDL to receive")
                else:
                    idl_input = self._split_idl_assignment(py_in, str_py_in,
                                                           assign_to)

        return idl_input


    def _split_idl_assignment(self, py_in, str_py_in, assign_to):
        """Take a very long numpy array and return a list of commands
        to execute in order to assign this value to an IDL variable."""

        if assign_to is None:
            print("(!) No assign_to set.", file=sys.stderr)

        idl_input = []
        extend_string = ''

        # Each assignment string must be no longer than max_idl_code_area:
        #      assign_to=[assign_to,<max_length>
        max_length = self.max_idl_code_area - 2 * len(assign_to) - 3

        # In addition, code area limited by number of elements of array
        max_n_elements = self.max_n_elements_code_area

        # Loop until string has been split up into manageable chunks
        array_string_remaining = str_py_in[1:]  # Remove '['
        while len(array_string_remaining) > 1:
            # Take the first max_n_elements elements (separated by
            # commas) of the first max_length characters of the string
            max_string = re.match('([^,]*[,\]]){,' + str(max_n_elements) + '}',
                                  array_string_remaining[:max_length]).group()

            # Remove final character (',' or ']') from max_string
            idl_input.append(assign_to + "=[" + extend_string +
                             max_string[:-1] + "]")

            # Not for the first time round
            extend_string = ''.join([assign_to, ","])

            # What's left?
            array_string_remaining = array_string_remaining[len(max_string):]

        if len(py_in.shape) > 1:
            # Convert data type and shape
            idl_input.append(assign_to + "=" + "reform(" +
                             self._idl_cast_from_dtype(py_in.dtype, assign_to)
                             + ", " + str(py_in.shape[::-1])[1:-1] + ")")
        else:
            # Convert data type
            idl_input.append(assign_to + "=" +
                             self._idl_cast_from_dtype(py_in.dtype, assign_to))

        return idl_input


    def _idl_cast_from_dtype(self, dtype, idl_str):
        """Take a NumPy dtype and return an expression to cast an IDL
        expression as appropriate type."""

        if dtype.name[0:3] == 'str':
            return idl_str

        # NaN and Inf
        idl_str = idl_str.replace('inf', '!values.f_infinity')
        idl_str = idl_str.replace('nan', '!values.f_nan')

        if dtype.name == 'bool':
            return "byte(" + str(int(eval(idl_str))) + ")"
        elif dtype.name == 'uint8':
            return "byte(" + idl_str + ")"
        elif dtype.name == 'int16':
            return "fix(" + idl_str + ")"
        elif dtype.name == 'uint16':
            return "uint(" + idl_str + ")"
        elif dtype.name == 'int32':
            return "long(" + idl_str + ")"
        elif dtype.name == 'uint32':
            return "ulong(" + idl_str + ")"
        elif dtype.name == 'int64':
            return "long64(" + idl_str.replace('L', 'LL') + ")"
        elif dtype.name == 'uint64':
            return "ulong64(" + idl_str.replace('L', 'LL') + ")"
        elif dtype.name == 'float8':  # Not a NumPy type?
            print("Warning: converting 8-bit to 32-bit float.", file=sys.stderr)
            return "float(" + idl_str + ")"
        elif dtype.name == 'float16':  # Not a NumPy type?
            print("Warning: converting 16-bit to 32-bit float.", file=sys.stderr)
            return "float(" + idl_str + ")"
        elif dtype.name == 'float32':
            return "float(" + idl_str + ")"
        elif dtype.name == 'float64':
            return "double(" + idl_str + ")"
        else:
            print("(!) Could not convert NumPy dtype", \
                  dtype.name, "to IDL.", file=sys.stderr)
            return


    def _python_to_idl_structure(self, python_input, assign_to):
        """Given a Python dictionary, or a (simple, 1D) list of dictionaries,
        return list of command(s) to assign IDL structure to assign_to."""

        # Convert from numpy array if necessary
        py_in = numpy.array(python_input).tolist()

        try:
            return self._python_to_idl_structure_short(py_in, assign_to)
        except IDLInputOverflowError:
            if assign_to is not None:
                return self._python_to_idl_structure_long(py_in, assign_to)
            else:
                raise


    def _python_to_idl_structure_short(self, py_in, assign_to):
        """Given a Python dictionary, or a (simple, 1D) list of dictionaries,
        return single command to assign IDL structure to assign_to."""

        if hasattr(py_in, 'keys'):
            py_in_list = [py_in]
        else:
            py_in_list = py_in

        idl_input_list = []
        for row in py_in_list:
            struct_fields = []
            for key in row:
                struct_fields.append(''.join(
                    [key, ':', self._python_to_idl_input(row[key])]))
            idl_input_list.append('{' + ','.join(struct_fields) + '}')

        if hasattr(py_in, 'keys'):
            idl_input = idl_input_list[0]
        else:
            idl_input = '[' + ','.join(idl_input_list) + ']'

        if assign_to is not None:
            idl_input = assign_to + '=' + idl_input

        if len(idl_input) > self.max_idl_code_area:
            # String too long!  Need to create list of shorter commands
            raise IDLInputOverflowError("Expression too long for IDL to receive")
        else:
            return idl_input


    def _python_to_idl_structure_long(self, py_in, assign_to):
        """Given a Python dictionary, or a (simple, 1D) list of dictionaries,
        return a list of commands to assign IDL structure to assign_to."""

        if hasattr(py_in, 'keys'):
            n_rows = 1
            py_in_row = py_in
        else:  # List of dictionaries
            n_rows = len(py_in)
            py_in_row = py_in[0]

        # Make one row
        struct_fields = []
        for key in py_in_row:
            struct_fields.append(''.join(
                [key, ':', self._python_to_idl_input(py_in_row[key])]))

        # Commands to execute
        idl_input = [assign_to + '={' + ','.join(struct_fields) + '}']

        if n_rows > 1:
            # Fill rows with data
            for row in py_in[1:]:
                struct_fields = []
                for key in row:
                    struct_fields.append(''.join(
                        [key, ':', self._python_to_idl_input(row[key])]))
                idl_input.append(assign_to + '=[' + assign_to + ',{'
                                     + ','.join(struct_fields) + '}]')

        return idl_input


    # Receiving from IDL


    def _get_from_save(self, names):
        """Get values for a list of names.

        Use save command in IDL to save arrays/structure into file and use
        readsav to read it into python objects. This will save a lot of time
        when handling large arrays.
        """
        if self._cache_dir is None:
            self._cache_dir=tempfile.mkdtemp()
        else:
            if sys.version_info.major < 3:
                try:
                    os.makedirs(self._cache_dir)
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
            else:
                os.makedirs(self._cache_dir, exist_ok=True)
        if names != ['pidly_tmp']:
            valuelist = self.variables()
            for eachname in names:
                if eachname.upper() not in valuelist:
                    raise NameError("name {} not in idl variable list".format(eachname))
        savePath = os.path.join(self._cache_dir, 'pidly.sav')
        todoStr = "save," + ",".join(names) + ",file='{}'".format(savePath)
        self.ex(todoStr)
        return idlio.readsav(savePath)


    def _get_from_print(self):
        """Get IDL's string representation of expression."""
        # Special treatment of arrays of type string (type == 7)
        # Floats/doubles printed with appropriate (or too much) precision
        idl_output = self.ex(
            # Arrays of strings
            'if size(pidly_tmp,/type) eq 7 and n_elements(pidly_tmp) gt 1 '
            + 'then print,strjoin(reform(pidly_tmp,n_elements(pidly_tmp)),'
            + '\'' + STR_DELIMITER + '\') '
            # Structures
            + 'else if n_elements(pidly_tmp) gt 0 '
            + 'and size(pidly_tmp,/type) eq 8 then '
            + 'for i = 0, n_elements(pidly_tmp) - 1 do begin'
            + '  print, "{" & '
            + '  for j = 0, n_tags(pidly_tmp) - 1 do '
            + '    if size(pidly_tmp[i].(j),/type) eq 5 then '
            + '      print,reform(pidly_tmp[i].(j), '
            + '        n_elements(pidly_tmp[i].(j))), format="(E26.18)" '
            + '    else if size(pidly_tmp[i].(j),/type) eq 4 then '
            + '      print,reform(pidly_tmp[i].(j),'
            + '         n_elements(pidly_tmp[i].(j))), format="(E19.11)" '
            + '    else if size(pidly_tmp[i].(j),/type) eq 7 '
            + '      and n_elements(pidly_tmp[i].(j)) gt 1 then '
            + '      print,strjoin(reform(pidly_tmp[i].(j),'
            + '                    n_elements(pidly_tmp[i].(j))),'
            + '      \'' + STR_DELIMITER + '\') '
            + '    else print, reform(pidly_tmp[i].(j),'
            + '         n_elements(pidly_tmp[i].(j))) & '
            + '  print, "}" &'
            + '  endfor '
            # Doubles
            + 'else if n_elements(pidly_tmp) gt 0 '
            + 'and size(pidly_tmp,/type) eq 5 then print,reform(pidly_tmp,'
            + 'n_elements(pidly_tmp)), format="(E26.18)" '
            # Floats
            + 'else if n_elements(pidly_tmp) gt 0 '
            + 'and size(pidly_tmp,/type) eq 4 then print,reform(pidly_tmp,'
            + 'n_elements(pidly_tmp)), format="(E19.11)" '
            # Other
            + 'else if n_elements(pidly_tmp) gt 0 then print,reform(pidly_tmp,'
            + 'n_elements(pidly_tmp))',
            print_output=False, ret=True)

        # Parse this string into a python variable
        if idl_output:
            idl_type_dims = self.ex(
                'print,size(pidly_tmp,/type) & '
                + 'print,size(pidly_tmp,/dimensions)',
                print_output=False, ret=True).splitlines()
            idl_type = int(idl_type_dims[0])
            idl_dims = numpy.array(
                ''.join(idl_type_dims[1:]).split()).astype(int)
            if idl_type == 8:  # Structure
                self.ex('help,pidly_tmp,/struct,output=pidly_tmp_2',
                        print_output=False)
                idl_struct = self.ev('pidly_tmp_2', use_cache=False).tolist()
                #idl_output = ''.join(['{', idl_output, '}'])
                return self._idl_struct_to_python(idl_output, idl_struct)
            else:
                return self._idl_output_to_python(idl_output, idl_type,
                                                  idl_dims)


    def _wait_for_prompt(self, print_output=False):
        """Read IDL output until IDL prompt displayed and return."""
        index = 1
        output_lines = []
        stop = False
        halt = False
        while index == 1:
            try:
                # 0: waiting for input, 1: output received, 2: IDL exited
                index = self.expect([self.idl_prompt, '\n', pexpect.EOF])
            except KeyboardInterrupt:
                print("\npIDLy: KeyboardInterrupt")
                if not _ipython:
                    print(self.before)
                    sys.stdout.flush()
                self.interact(show_prompt=False)
                break
            if sys.version_info.major < 3:
                new_line = self.before.replace('\r', '')
            else:
                new_line = self.before.decode().replace('\r', '')
            if new_line.startswith('% Stop encountered:'):
                stop = True
            if new_line.startswith('% Execution halted at:'):
                halt = True
            if new_line:
                output_lines.append(new_line)
                if print_output:
                    print(new_line)  # Live output while waiting for prompt
        if halt or stop:
            if not print_output:  # print output anyway
                print('\n'.join(output_lines))
            self.interact()
        return '\n'.join(output_lines)


    def _idl_output_to_python(self, idl_output, idl_type, idl_dims):
        """Take output from IDL print statement and return value."""

        # Find Python dtype and shape
        dtype = self._dtype_from_idl_type(idl_type)  # = None for string
        shape = self._shape_from_idl_dims(idl_dims)

        # Split the output into separate items
        if idl_type == 7:
            value = idl_output.split(STR_DELIMITER)
        else:
            value = idl_output.split()

        if value:
            if idl_type == 7:  # String
                if shape == ():
                    # Concatenate string
                    value = ' '.join(value)

            # Convert to numpy.array of appropriate type
            if dtype is None:
                value = numpy.array(value)
            else:
                value = numpy.array(value).astype(dtype)

            # Reshape array
            if numpy.product(shape) != value.size:
                print("(!) Could not reshape array.", file=sys.stderr)
            else:
                value = value.reshape(shape)

            return value


    def _dtype_from_idl_type(self, idl_type):
        """Convert IDL type to numpy dtype."""

        if idl_type is not None:
            python_idl_types = [
                None, 'uint8', 'int16', 'int32', 'float32', 'float64',
                None, None, None, None, None,
                None, 'uint16', 'uint32', 'int64', 'uint64']
            dtype = python_idl_types[idl_type]
            if dtype is None and idl_type != 7:
                print("(!) Could not convert IDL type " \
                      + str(idl_type) + " to Python.", file=sys.stderr)
        else:
            dtype = None
        return dtype


    def _shape_from_idl_dims(self, idl_dims):
        """Convert IDL dimensions to numpy shape."""

        # Dimensions run the opposite way, Python vs IDL
        shape = numpy.array(idl_dims[::-1]).tolist()
        if shape == [0]:
            shape = []

        return tuple(shape)


    def _idl_struct_to_python(self, idl_output, idl_struct):
        """Take print output of IDL structure and return Python dictionary.

        No spaces allowed in field names.

        """
        # Create meta-dictionary
        dict_def = []
        j = 1  # Omit first line
        for i in range(1, int(idl_struct[0].split()[3]) + 1):  # Number of tags
            # For each field, store (name, dtype, shape) in the meta-dictionary
            idl_struct_split = idl_struct[j].replace(', ',',').split()
            name = idl_struct_split[0]
            try:
                idl_type = self._idl_type_from_name(idl_struct_split[1])
                idl_dims = self._dims_from_struct_help(idl_struct_split[2])
            except IndexError:  # Some descriptions span two lines
                j += 1
                idl_type = self._idl_type_from_name(idl_struct_split[0])
                idl_dims = self._dims_from_struct_help(idl_struct_split[1])
            dict_def.append((name, idl_type, idl_dims))
            j += 1

        dict_list = []
        idl_output = ''.join(idl_output)
        rows = idl_output[2:-1].split('\n}\n{\n')  # Remove {\n and } at ends
        for row in rows:  # Each structure in array of structures
            # Create a dictionary for each row
            items = row.splitlines()
            dict_row = {}
            for name, idl_type, idl_dims in dict_def:
                idl_out_list = []
                if idl_type == 7:  # String
                    idl_out = items.pop(0)
                else:
                    line_items = items.pop(0).split()
                    for i in range(max(numpy.product(idl_dims), 1)): # No. values
                        try:
                            idl_out_list.append(line_items.pop(0))
                        except IndexError:
                            line_items = items.pop(0).split()
                            idl_out_list.append(line_items.pop(0))
                    idl_out = ' '.join(idl_out_list)
                dict_row[name.lower()] = self._idl_output_to_python(
                    idl_out, idl_type, idl_dims)
            dict_list.append(dict_row)

        if len(dict_list) == 1:
            return dict_list[0]
        else:
            return dict_list


    def _idl_type_from_name(self, type):
        """Return IDL type code, given type name."""

        if type == 'BYTE':
            return 1
        elif type == 'INT':
            return 2
        elif type == 'LONG':
            return 3
        elif type == 'FLOAT':
            return 4
        elif type == 'DOUBLE':
            return 5
        elif type == 'COMPLEX':
            return 6
        elif type == 'STRING':
            return 7
        elif type == 'STRUCT':
            return 8
        elif type == 'DCOMPLEX':
            return 9
        elif type == 'POINTER':
            return 10
        elif type == 'OBJREF':
            return 11
        elif type == 'UINT':
            return 12
        elif type == 'ULONG':
            return 13
        elif type == 'LONG64':
            return 14
        elif type == 'ULONG64':
            return 15


    def _dims_from_struct_help(self, struct_help):
        """Return IDL dims given description from structure."""

        try:
            dims = re.search('(?<=Array\[).*\]', struct_help).group()[:-1]
            idl_dims = numpy.array(dims.split(',')).astype(int)
            return idl_dims
        except AttributeError:
            return [0]


class TestPidly(unittest.TestCase):
    """Unit tests for pIDLy."""

    def setUp(self):
        if len(sys.argv) > 1 and sys.argv[0].endswith('test_pidly.py'):
            self.idl = IDL(sys.argv[1])
        else:
            self.idl = IDL()
        self.start_time = None
        self.mid_time = None
        self.end_time = None

    def tearDown(self):
        def t(time_delta):
            return time_delta.seconds + time_delta.microseconds / 1000000.
        self.idl.close()
        try:
            print("%0.5ss/%0.5ss " % (t(self.mid_time - self.start_time),
                                      t(self.end_time - self.mid_time)), end=' ')
            sys.stdout.flush()
        except TypeError:
            try:
                print("%0.5ss " % t(self.end_time - self.start_time), end=' ')
                sys.stdout.flush()
            except TypeError:
                pass

    def sendAndReceive(self, x):
        self.start_time = now()
        self.idl.x = x
        self.mid_time = now()
        y = self.idl.ev('x')
        #y = self.idl.x  # Twice as slow!
        self.end_time = now()
        return y

    def test_idl_dead(self):
        self.idl.close()
        with self.assertRaises(OSError):
            self.idl('print, 1')

    def test_longest_line(self):
        s = ["x='"]
        for i in range(self.idl.max_sendline - 4):
            s.append('a')
        s.append("'")
        x = ''.join(s)
        self.start_time = now()
        self.idl.sendline(x)
        self.idl.expect('IDL> ')
        self.mid_time = now()
        y = self.idl.x
        self.end_time = now()
        self.assertEqual(y, x[3:-1])

    def test_longest_string(self):
        n = self.idl.max_idl_code_area - 10
        x = ''.join(["x='"] + ["a" for i in range(n)] + ["'"])
        self.start_time = now()
        self.idl(x)
        self.mid_time = now()
        y = self.idl.x
        self.end_time = now()
        self.assertEqual(y, x[3:-1])

    def test_longest_string_overflow(self):
        s = ["x='"]
        for i in range(self.idl.max_idl_code_area - 3):
            s.append('a')
        s.append("'")
        x = ''.join(s)
        self.assertRaises(IDLInputOverflowError, self.idl, x)

    def test_20_function_calls(self):
        x = numpy.random.random(20)
        y = numpy.zeros(20)
        self.start_time = now()
        for i in range(20):
            y[i] = self.idl.sin(x[i])
        self.end_time = now()
        self.assertEqual(y.tolist(), numpy.sin(x).tolist())

    def test_20_function_calls_explicit(self):
        x = numpy.random.random(20)
        y = numpy.zeros(20)
        self.start_time = now()
        for i in range(20):
            y[i] = self.idl.func('sin', x[i])
        self.end_time = now()
        self.assertEqual(y.tolist(), numpy.sin(x).tolist())

    def test_20_function_calls_really_explicit(self):
        x = numpy.random.random()
        self.idl.x = x
        self.start_time = now()
        for i in range(20):
            y = self.idl.ev('sin(x)')
        self.end_time = now()
        self.assertEqual(y, numpy.sin(x))

    def test_long_function_call(self):
        x = numpy.random.random(1000)
        self.start_time = now()
        y = self.idl.total(x)
        self.end_time = now()
        self.assertEqual(y, sum(x))

    def test_long_function_dicts(self):
        x = [{'a':numpy.random.random()} for i in range(100)]
        self.start_time = now()
        y = self.idl.n_elements(x)
        self.end_time = now()
        self.assertEqual(y, len(x))

    def test_longish_function_call(self):
        # Total length raises IDLInputOverflowError, but each arg is short
        x = numpy.random.random(84)
        self.start_time = now()
        y = self.idl.total(x, double=True)
        self.end_time = now()
        self.assertEqual(y, sum(x))

    def test_single_integer(self):
        x = 2
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_single_element_list_int(self):
        x = [2]
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_single_float(self):
        x = 2.123
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_single_float32(self):
        x = numpy.array(2.123, dtype='float32')
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertEqual(numpy.array(x).shape, y.shape)
        self.assertEqual(numpy.array(x).dtype, y.dtype)

    def test_huge_double(self):
        x = -1e100
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertEqual(numpy.array(x).shape, y.shape)
        self.assertEqual(numpy.array(x).dtype, y.dtype)

    def test_infinity(self):
        x = numpy.inf
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertTrue(numpy.isinf(x))
        self.assertTrue(numpy.isinf(y))
        self.assertEqual(numpy.array(x).shape, y.shape)
        self.assertEqual(numpy.array(x).dtype, y.dtype)

    def test_infinity_neg(self):
        x = -numpy.inf
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertTrue(numpy.isneginf(x))
        self.assertTrue(numpy.isneginf(y))
        self.assertEqual(numpy.array(x).shape, y.shape)
        self.assertEqual(numpy.array(x).dtype, y.dtype)

    def test_nan(self):
        x = numpy.nan
        y = self.sendAndReceive(x)
        # NB NaN != NaN
        self.assertTrue(numpy.isnan(x))
        self.assertTrue(numpy.isnan(y))
        self.assertEqual(numpy.array(x).shape, y.shape)
        self.assertEqual(numpy.array(x).dtype, y.dtype)

    def test_inf_nan_array(self):
        x = [1.2, numpy.nan, numpy.inf, -numpy.inf, 3, numpy.nan, 4]
        y = self.sendAndReceive(x)
        # NB NaN != NaN
        self.assertTrue(all(numpy.isnan([x[1], x[5]])))
        self.assertTrue(all(numpy.isnan([y[1], y[5]])))
        self.assertEqual(x[0::2], y.tolist()[0::2])
        self.assertEqual(x[3], y[3])
        self.assertEqual(numpy.array(x).shape, y.shape)
        self.assertEqual(numpy.array(x).dtype, y.dtype)

    def test_single_string(self):
        x = 'hello'
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_multi_word_string(self):
        x = 'hello there'
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_list_of_strings(self):
        x = [' 5 2  5k 2', '4  33 55 1 ', ' 4 ', '2', '  3   2']
        y = self.sendAndReceive(x)

        self.assertEqual(y.tolist(), x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_3d_list_of_strings(self):
        x = [[['b', ' d   s    '], ['ff ', 's d  ewew']],
             [['', 'f'], ['gs', 'a']]]
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_long_integer(self):
        x = 25525252525525
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_int_array(self):
        x = [1,2,4,2555]
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_long_int_array(self):
        x = [1,2,3,25525252525525,23, 5, 6, 5, 2, 5, 7, 8, 3, 5]
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_2d_int_array(self):
        x = [[1,2,3],[3,4,5]]
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_3d_int_array(self):
        x = [[[ 1, 2, 3],[ 4, 5, 6]],[[ 7, 8, 9],[10,11,12]],
             [[13,14,15],[16,17,18]],[[19,20,21],[22,23,24]]]
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x)
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_mixed_array_warning(self):
        x = [22,23.,'24']
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), numpy.array(x).tolist())
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_simple_dictionary(self):
        x = dict(a=2)
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        self.assertEqual(numpy.array(x['a']).dtype, y['a'].dtype)
        self.assertEqual(numpy.array(x).shape, numpy.array(y).shape)

    def test_3_element_dict(self):
        x = {'a':'a h ', 'b':1e7, 'c':0.7}
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        for key in x:
            self.assertEqual(numpy.array(x[key]).dtype, y[key].dtype)
        self.assertEqual(numpy.array(x).shape, numpy.array(y).shape)

    def test_dict_string_arrays(self):
        x = {'a': numpy.arange(2*5*1).reshape(2,5,1),
             'b': [[['b', ' d   s    '], ['ff ', 's d  ewew']],
                   [['', 'f'], ['gs', 'a']]],
             'c': 5, 'd': [1, 5, 2.3]}
        y = self.sendAndReceive(x)
        self.assertEqual(sorted(x.keys()), sorted(y.keys()))
        self.assertEqual([y[key].tolist() for key in y],
                         [numpy.array(x[key]).tolist() for key in y])
        self.assertEqual([y[key].dtype for key in y],
                         [numpy.array(x[key]).dtype for key in y])
        self.assertEqual([y[key].shape for key in y],
                         [numpy.array(x[key]).shape for key in y])

    def test_3_3_element_dicts(self):
        x = [{'a':'ah', 'b':1000, 'c':0.7}, {'a':'be', 'b':8, 'c':-6.3},
             {'a':'x', 'b':0, 'c':81.}]
        y = self.sendAndReceive(x)
        self.assertEqual(y, x)
        for i, d in enumerate(x):
            for key in d:
                self.assertEqual(numpy.array(d[key]).dtype,
                                 y[i][key].dtype)
        self.assertEqual(numpy.array(x).shape, numpy.array(y).shape)

    def test_100_dicts_float32_double_string(self):
        x = [{'a':numpy.random.random(2).astype('float32'),
              'b':numpy.random.random(1),
              'c':(numpy.random.random(1)*100).astype('str')}
             for i in range(100)]
        y = self.sendAndReceive(x)
        for i, d in enumerate(x):
            self.assertEqual([y[i][key].tolist() for key in y[i]],
                             [d[key].tolist() for key in d])
        for i, d in enumerate(x):
            for key in d:
                self.assertEqual(d[key][0].dtype,
                                 y[i][key][0].dtype)
        self.assertEqual(numpy.array(x).shape, numpy.array(y).shape)


    def test_4d_int_array(self):
        x = numpy.arange(2*3*4*5).reshape(2,3,4,5)
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x.tolist())
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_6d_int_array_tests_max_n_elements_code_area(self):
        x = numpy.arange(2*3*4*5*6*7).reshape(2,3,4,5,6,7)
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x.tolist())
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_8d_int_array(self):
        x = numpy.arange(2*3*1*1*4*5*6*7).reshape(2,3,1,1,4,5,6,7)
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x.tolist())
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_50_doubles(self):
        x = numpy.random.random(50)
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x.tolist())
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_50_float32(self):
        x = numpy.random.random(50).astype('float32')
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x.tolist())
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_50_doubles_1e30(self):
        x = numpy.random.random(50) * 1e30
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x.tolist())
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_50_float32_1e30(self):
        x = numpy.random.random(50).astype('float32') * 1e30
        y = self.sendAndReceive(x)
        self.assertEqual(y.tolist(), x.tolist())
        self.assertEqual(numpy.array(x).dtype, y.dtype)
        self.assertEqual(numpy.array(x).shape, y.shape)

    def test_speed_20000_doubles(self):
        x = numpy.random.random(20000)
        y = self.sendAndReceive(x)
        # Don't print all 20000 floats if fails!
        self.assertTrue(y.tolist() == x.tolist())
        self.assertEqual(x.dtype, y.dtype)
        self.assertEqual(y.shape, x.shape)

    def test_speed_5000_doubles_one_by_one(self):
        n = 5000
        x = numpy.random.random(n)
        self.start_time = now()
        self.idl('x = dblarr(' + str(n) + ')')
        for i in range(n):
            self.idl('x[' + str(i) + ']', x[i])
            if (i + 1) % 100 == 0:
                print((i + 1), end=' ')
                sys.stdout.flush()
        self.mid_time = now()
        y = self.idl.x
        self.end_time = now()
        self.assertTrue(y.tolist() == x.tolist())
        self.assertEqual(y.dtype, x.dtype)
        self.assertEqual(y.shape, x.shape)

    def test_ev_with_cache(self):
        n = 5000
        seed = 1
        self.idl('y = randomu({}, {})'.format(seed, n))
        self.start_time = now()
        y_from_ev = self.idl.y
        self.mid_time = now()
        y_from_ev_cache = self.idl.ev('y', use_cache=True)
        self.end_time = now()
        self.assertEqual(y_from_ev.shape, (n,))
        self.assertEqual(y_from_ev_cache.shape, (n,))

    def test_ev_list(self):
        n = 5000
        seed = 1
        self.idl('x = randomu({}, {})'.format(seed, n))
        self.idl('y = randomu({}, {})'.format(seed, n))
        self.start_time = now()
        x_from_ev = self.idl.x
        y_from_ev = self.idl.y
        self.mid_time = now()
        xy_from_ev_list = self.idl.ev_list(['x', 'y'])
        self.end_time = now()
        # Don't print all values if fails!
        self.assertTrue(x_from_ev.tolist() == xy_from_ev_list['x'].tolist())
        self.assertTrue(y_from_ev.tolist() == xy_from_ev_list['y'].tolist())
        self.assertEqual(xy_from_ev_list['x'].shape, (n,))
        self.assertEqual(xy_from_ev_list['y'].shape, (n,))

    def test_ev_list_with_cache(self):
        n = 5000
        seed = 1
        self.idl('x = randomu({}, {})'.format(seed, n))
        self.idl('y = randomu({}, {})'.format(seed, n))
        self.start_time = now()
        xy_from_ev_list = self.idl.ev_list(['x','y'])
        self.mid_time = now()
        self.idl.use_cache = True
        xy_from_ev_list_cache = self.idl.ev_list(['x','y'])
        self.end_time = now()
        # Don't print all values if fails!
        self.assertTrue(xy_from_ev_list_cache['x'].tolist() == xy_from_ev_list['x'].tolist())
        self.assertTrue(xy_from_ev_list_cache['y'].tolist() == xy_from_ev_list['y'].tolist())
        self.assertEqual(xy_from_ev_list['x'].shape, (n,))
        self.assertEqual(xy_from_ev_list['y'].shape, (n,))

def test():
    """Run full tests on pIDLy."""
    # Use python pidly.py or python3 pidly.py.

    if len(sys.argv) > 1 and sys.argv[0].endswith('pidly.py'):
        idl = IDL(sys.argv[1])
    else:
        idl = IDL()

    print("pIDLy", __version__ + ": running full tests.")
    print("IDL", end=' ')
    idl('print,!VERSION')
    print("pexpect", pexpect.__version__)
    print("Showing (Python -> IDL time) / (IDL -> Python time).\n")

    import doctest
    import pidly

    suite = unittest.TestLoader().loadTestsFromTestCase(TestPidly)
    suite.addTest(doctest.DocTestSuite(pidly))
    unittest.TextTestRunner(verbosity=2).run(suite)


if __name__ == "__main__":
    test()
