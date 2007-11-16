# $Id: ts.py,v 1.1 2007-11-16 13:57:04 joishi Exp $
#
# read time_series.dat and return a TimeSeries class of 1D numpy
# arrrays
#
#
import sys
import os.path 
import re
import numpy as N

class TimeSeries:
    """
    TimeSeries -- holds pencil code time series data. each variable is
    represented by a data member of the class.
    """
    pass

def read_ts(filename='time_series.dat',datadir='data/',double=0,print_std=0,quiet=0):
    """
    read_ts -- reads Pencil Code time series data. modeled after idl function of same name.
    """
    datadir = os.path.expanduser(datadir)
    infile = open(datadir+filename,"r")
    lines = infile.readlines()
    infile.close()
    
    # need to handle cases where restart AND print.in changes, but not right away
    # idl version uses input_table function with a STOP_AT and FILEPOSITION keywords
    nlines_init=len(lines)
    keys=[]
    data = N.zeros((nlines_init,len(keys)))
    nlines=0
    for line in lines:
        if (re.search("^#--",line)):
            # read header and create keys for dictionary
            line = line.strip("#-\n")
            keys_new = re.split("-+",line)
            if (keys_new != keys):
                n_newrows = abs(len(keys_new) - len(keys))
                data = N.append(data,N.zeros((nlines_init,n_newrows)),axis=1)
                keys = keys_new
        else:
            row = N.array(map(float,re.split(" +",line.strip(" \n"))))
            data[nlines,:] = row
            nlines += 1
    
    #clean up data
    data = N.resize(data,(nlines,len(keys)))
    
    if (not quiet):
        print "Read",nlines,"lines."

    #assemble into a TimeSeries class
    obj = TimeSeries()
    for i in range(0,len(keys)):
        setattr(obj,keys[i],data[:,i])

    return obj
            
if __name__=='__main__':
    read_ts.__doc__
