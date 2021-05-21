#!/usr/bin/env python
#coding: utf-8
# structure_function.py
# Written by James F. Hollins
# included by Frederick Gent (fred.gent.ncl@gmail.com)

import numpy as np
import scipy as sp

def structure_function(arr,y,x):
    n=int(np.size(y)/2)
    m=int(np.size(x)/2)
    D=np.zeros((n,m))
    for i in range(0, n, 1):
        yshift=np.roll(arr,i,axis=0)
        for j in range(0, m, 1):
            xshift=np.roll(yshift,j,axis=1)
            shift_diff=xshift-arr
            shift_squared=np.power(shift_diff,2)
            shift_squared_cut=shift_squared[:,int(m/2):m]
            shift_squared_cut=\
                shift_squared_cut[np.logical_not(np.isnan(shift_squared_cut))]
            if np.size(shift_squared_cut)==0:
                D[i,j]=float('nan')
            else:
                D[i,j]=np.mean(shift_squared_cut)
    return D

def structure_function_shift_range(arr,nmin,nmax,x):
    n=nmax-nmin
    m=int(np.size(x)/2)
    D=np.zeros((n,m))
    for i in range(nmin, nmax, 1):
        yshift=np.roll(arr,i,axis=0)
        for j in range(0, m, 1):
            xshift=np.roll(yshift,j,axis=1)
            shift_diff=xshift-arr
            shift_squared=np.power(shift_diff,2)
            shift_squared_cut=shift_squared[:,int(m/2):m]
            shift_squared_cut=\
                shift_squared_cut[np.logical_not(np.isnan(shift_squared_cut))]
            if np.size(shift_squared_cut)==0:
                D[(i-nmin),j]=float('nan')
            else:
                D[(i-nmin),j]=np.mean(shift_squared_cut)
    return D
