# $Id: ts.py 13781 2011-05-04 06:25:24Z iomsn $
#
# Two routines to compute and read data computed from the extisting VAR files.
# Meant to be used after the simulation in case not all diagnostic variables
# were plotted.
#
# Author: Simon Candelaresi (iomsn@physto.se, iomsn1@googlemail.com).
# 
#

import pencil as pc
import numpy as np
import struct

def post_compute(variables = ['b2m'], datadir = 'data'):
    """
    Compute diagnostic variables from the VAR files.

    call signature::
    
      post_compute(variables = ['urms'], datadir = 'data')
    
    Read the VAR files and compute the diagnostic variables.
    Write the result in 'post_evaluation.dat'.
    
    Keyword arguments:
    
      *variables*:
        The diagnostic variables to be computed.
        
      *datadir*:
        Data directory.
        
    """

    # open the destination file for writing
    fd = open(datadir + '/post_evaluation.dat', 'w')
    
    # write the header
    fd.write('#--t-------------')
    for variable in variables:
        fd.write('--{0:15}'.format((variable+'--------------------')[:15]))
    fd.write('\n')
    
    # read the list of all VAR files
    var_list_file = open('data/proc0/varN.list')
    var_list = var_list_file.readlines()
    var_list_file.close()
    
    if (len(variables[0]) == 1):        
        variables = [variables]
        
    # check if jj needs to be computed
    bb_flag = False
    # array containing the variables which depend on bb
    b_dep = ['ab_int', 'jb_int', 'b1m', 'b2m', 'bm2', 'abm', 'abrms', 'jbm', 'brms']
    if (len(set(variables + b_dep)) < len(variables + b_dep)):
        bb_flag = True
        
    # check if jj needs to be computed
    jj_flag = False
    # array containing the variables which depend on jj
    j_dep = ['jb_int', 'j2m', 'jm2', 'jbm', 'jrms']
    if (len(set(variables + j_dep)) < len(variables + j_dep)):
        jj_flag = True
    
    for var_file in var_list[:3]:
        # read the var file
        var = pc.read_var(varfile = var_file[:-1], datadir = datadir, quiet = True)
        # the output string which will be written in the destination file
        out_string = '{0:1.9e}  '.format(np.float64(var.t))
        if bb_flag:
            bb = pc.curl(var.aa, var.dx, var.dy, var.dz)
        if jj_flag:
            jj = pc.curl2(var.aa, var.dx, var.dy, var.dz)
            
        for variable in variables:
            if variable == 'b2m':
                b2m = (pc.dot2(bb)).mean()
                out_string += '{0:1.9e}  '.format(b2m)
            elif variable == 'j2m':
                j2m = (pc.dot2(jj)).mean()
                out_string += '{0:1.9e}  '.format(j2m)
                break
            elif variable == 'abm':
                abm = (pc.dot(var.aa, bb)).mean
                out_string += '{0:1.9e}  '.format(b2m)
                    
        fd.write(out_string+'\n')
            
    fd.close()



#  Under Construction

#def read_post(datadir = 'data'):
    #"""
    #Read the post processed diagnostic variables.

    #call signature::
    
      #read_post(datadir = 'data')
    
    #Read the post processed diagnostic variables from 'data/post_evaluation.dat.
    #Return an object with the variables.
    
    #Keyword arguments:
        
      #*datadir*:
        #Data directory.
        
    #"""

    ## the class which contains the post processing data
    #class post:
        #def __init__(self):
    ## open the destination file for writing
    #fd = open(datadir + '/post_evaluation.dat', 'r')
    
    ## check 

    #fd.close()

