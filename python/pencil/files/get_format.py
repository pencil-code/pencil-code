# $Id$
#
# check data format from data/var.general
#
# Check the endianness of the machine that wrote the data
#
# Author: C.P. McNally
# 

def get_format(datadir = 'data/'):
     
    if 'format = lsb' in open(datadir+'var.general').read():
      format = 'little'
    elif 'format = msb' in open(datadir+'var.general').read():
      format = 'big'
    else :
      format = 'little'

    return format
