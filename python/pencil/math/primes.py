#!/usr/bin/env python
#coding: utf-8
# primes.py
# Written by Fred Gent (fred.gent.ncl@gmail.com)
"""
Find the prime factorization of the simulation grid, useful for checking 
permissible remesh parameters, parallelization options and fft compatibility.
"""
import numpy as np
from . import natural_sort
def prime_factors(n):
    i = 2
    factors = list()
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors

def divisors(factors, nmax, nmin=8):
    divisors = list()
    nn = len(factors)
    n = factors[0]
    divisors.append(n)
    for i in range(0,nn-1):
        n = factors[i-1] * n
        if n <= nmax/nmin and not n in divisors and np.mod(nmax,n)==0:
            divisors.append(n)
        for j in range(i,nn):
            m = factors[j]*factors[i]
            if m <= nmax/nmin and not m in divisors and np.mod(nmax,m)==0:
                divisors.append(m)
            for k in range(j+1,nn):
                p = factors[k]*factors[i-1]*factors[j]
                if p <= nmax/nmin and not p in divisors and np.mod(nmax,p)==0:
                    divisors.append(p)
    return natural_sort(divisors)

def common_factors(nx,ny,nz,nmin=8):
    factors = list()
    xfactors = prime_factors(nx)
    yfactors = prime_factors(ny)
    zfactors = prime_factors(nz)
    xdivisors = divisors(xfactors,nx,nmin=nmin)
    ydivisors = divisors(yfactors,ny,nmin=nmin)
    zdivisors = divisors(zfactors,nz,nmin=nmin)
    print('divisors x:',xdivisors,'y:',ydivisors,'z:',zdivisors)

    return [xdivisors, ydivisors, zdivisors]
           #[natural_sort(xyproc), natural_sort(yzproc)],\
           
