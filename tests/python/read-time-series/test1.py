#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Read time series."""

import numpy as np
import os
import sys
import pencil as pc

input_dir = './input'


def main(args):
    ts = read_ts(os.path.join(input_dir, 'serial-1'))
    write_summary('test1.out', ts)


def read_ts(datadir, filename='time-series.dat'):
    return pc.read_ts(
        datadir=datadir, filename=filename, plot_data=False, quiet=True
    )


def write_summary(filename, ts):
    outputdir = os.path.dirname(sys.argv[0])
    output = open(os.path.join(outputdir, filename), 'w')
    for var in ts.keys:
        values = getattr(ts, var)
        output.write('min(%s): %g\n' % (var, values.min()))
        output.write('avg(%s): %g\n' % (var, values.mean()))
        output.write('max(%s): %g\n' % (var, values.max()))
        output.write('std(%s): %g\n' % (var, values.std()))
        output.write('\n')
    output.close()

if __name__ == '__main__':
    main(sys.argv[1:])
