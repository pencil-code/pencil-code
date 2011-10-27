# $Id$
#
# read time_series.dat and return a TimeSeries class of 1D numpy
# arrrays
#
#
import sys
import os.path 
import re
import numpy as N
import pylab as P

class read_ts:
    """
    read_ts -- holds pencil code time series data. each variable is
    represented by a data member of the class.
    """

    def __init__(self, filename='time_series.dat', datadir='data',
                 double=0, print_std=0, quiet=0, plot_data=True, comment_char='#'):
        """
        constructor:
        -----------
        __init__ -- reads Pencil Code time series data. 
        Modeled after idl function of same name.

        params:
        ______
         string: filename  ='time_series.dat'
         string: datadir   = 'data'
         logical: double    = 0
         logical: print_std = 0
         logical: quiet     = 0
        """
        try:
            datadir = os.path.expanduser(datadir)
            infile = open(datadir+'/'+filename, "r")
            lines = infile.readlines()
            infile.close()
        except IOError:
            return
    
        # need to handle cases where restart AND print.in changes, 
        # but not right away
        # idl version uses input_table function with a STOP_AT 
        # and FILEPOSITION keywords
        nlines_init = len(lines)
        keys = []
        data = N.zeros((nlines_init, len(keys)))
        nlines = 0
        for line in lines:
            if re.search("^%s--" % comment_char, line):
                # read header and create keys for dictionary
                line = line.strip("%s-\n" % comment_char)
                keys_new = re.split("-+", line)
                if keys_new != keys:
                    n_newrows = abs(len(keys_new) - len(keys))
                    data = N.append(data, N.zeros((nlines_init, n_newrows)),
                                 axis=1)
                    keys = keys_new
            else:
                try:
                    row = N.array(map(float, re.split(" +", line.strip(" \n"))))
                    data[nlines, :] = row
                    nlines += 1
                except ValueError:
                    print "Invalid data on line %i. Skipping." % nlines
        #clean up data
        data = N.resize(data,(nlines,len(keys)))
    
        if (not quiet):
            print "Read",nlines,"lines."

        #assemble into a TimeSeries class
        for i in range(0,len(keys)):
            setattr(self,keys[i],data[:,i])
       
        if hasattr(self,'t') and plot_data: self.plot()

            
    def plot(self):
        """
          plot:
          ----
            Do two plots:
            Try to plot urms(t), ruzm(t) and brms(t), if any of these
            three is not available or zero, fill the list with the first two 
            variables other than `it' and `dt*'
        """
# speed the graphics (in connection with an ending P.show())
        P.ioff()
        elim = re.compile(r'dt|it|__|plot')
        # every argument of the read_ts class is listed in listargs
        listargs = dir(self) 
        # to eliminate it, dt*, and __*__ names
        for item in dir(self):
            if re.match(elim,item):
                listargs.remove(item)
        cnt = 0
        if (hasattr(self, 'urms') and self.urms.max() != 0.):
            cnt += 1
            P.subplot(2, 1, cnt)
            P.semilogy(self.t, self.urms)
            P.xlabel('Time')
            P.ylabel('urms')
            listargs.remove('urms')
        if (hasattr(self, 'brms') and self.brms.max() != 0.):
            cnt += 1
            P.subplot(2, 1, cnt)
            P.semilogy(self.t, self.brms)
            P.xlabel('Time')
            P.ylabel('brms')
            listargs.remove('brms')
        if (hasattr(self, 'ruzm') and self.ruzm.max() != 0.):
            cnt += 1
            P.subplot(2, 1, cnt)
            P.plot(self.t, self.ruzm)
            P.xlabel('Time')
            P.ylabel('ruzm')
            listargs.remove('ruzm')
        else:
            listargs.remove('t')
            i = 0
            while cnt <= 1:
                cnt += 1
                P.subplot(2, 1, cnt)
                P.plot(self.t, getattr(self, listargs[i]))
                P.xlabel('Time')
                P.ylabel(listargs[i])
                i += 1 
        P.show()
        P.ion()
     
if __name__=='__main__':
    read_ts.__doc__
