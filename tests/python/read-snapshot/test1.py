#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Read time series."""

import os
import sys
import pencil as pc

input_dir = '../../input'


def main(args):
    var = read_var(os.path.join(input_dir, 'serial-1'))
    write_summary('test1.out', var)


def read_var(datadir, filename='var.dat'):
    return pc.read_var(
        datadir=datadir, varfile=filename, trimall=True, quiet=True
    )


def write_summary(filename, data):
    outputdir = os.path.dirname(sys.argv[0])
    output = open(os.path.join(outputdir, filename), 'w')
    for var in 'lnrho', 'ss', 'ux', 'uy', 'uz':
        values = getattr(data, var)
        output.write('min(%s): %g\n' % (var, values.min()))
        output.write('avg(%s): %g\n' % (var, values.mean()))
        output.write('max(%s): %g\n' % (var, values.max()))
        output.write('std(%s): %g\n' % (var, values.std()))
        output.write('\n')
    output.close()


if __name__ == '__main__':
    main(sys.argv[1:])

# 'dx', 'dy', 'dz'
# 'lnrho', 'ss', 'uu', 'ux', 'uy', 'uz',
# 'x', 'y', 'z'
# 't'
