# $Id$
"""
param.py -- reads parameters from the f90 namelists
Author: J. Oishi (joishi@amnh.org)

REQUIRES: nl2python perl script (based on Wolfgang Dobler's nl2idl script)

todo: think about single/double precision; make things into numpy arrays?
"""
import numpy
import os

def read_param(datadir='data/',param2=False,quiet=False,asdict=False,nestdict=False):


    datadir = os.path.expanduser(datadir)

    if (param2):
        filen = os.path.join(datadir,'param2.nml')
    else:
        filen = os.path.join(datadir,'param.nml')

    # will need to read dim.dat to deal with precision, should that be necessary
    #dim = read_dim(datadir)

    # execute output of nl2python script
    if (not os.path.exists(filen)):
        #print "read_param: no such file",filen # Python 2
        print("read_param: no such file " + filen)
        raise ValueError

    if asdict:
        import paramdict
        if nestdict: ret = paramdict.ReadNML(filen,nest=True)
        else: ret = paramdict.ReadNML(filen)
        #if not quiet: print ret # Python 2
        if not quiet: print(ret)
    else:
        cmd = 'nl2python '+filen
        script = os.popen(cmd).read()
        #if (not quiet): print script; # Python 2
        if (not quiet): print(script);
        if (len(script) != 0):
            class Params: pass
            exec(script.replace("\n    ","\nParams.")[198:])
        else:
            #print "read_param: nl2python returned nothing! is $PENCIL_HOME/bin in the path?" # Python 2
            print("read_param: nl2python returned nothing! is $PENCIL_HOME/bin in the path?")
            return -1
        
        ret = Params() # par is the name of the class

    return ret

if __name__ == '__main__':
   #print self.__doc__ # Python 2
   print(self.__doc__)
