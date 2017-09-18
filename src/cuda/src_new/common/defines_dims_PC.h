/*
*   Contains the defines that are extracted from Pencil Code
*   and read to the config files in config.cc
*/
#pragma once

//Astaroth requires nghost to be know at compile-time:
//Otherwise 
//  a) It seriously limits the optimizations we could do
//  b) Calling device functions with dynamically allocated memory is wonky (CUDA 6)
//  c) The implementations will break with bound size other than 3 (2017-07-20) 
//#define PC_BOUND_SIZE nghost

#define PC_NX nx
#define PC_NY ny
#define PC_NZ nz

#define PC_MX mx 
#define PC_MY my
#define PC_MZ mz

#define PC_W_GRID_SIZE nw
#define PC_GRID_SIZE   mw

#define PC_NX_MIN (l1-1)//Parentheses for safety, f.ex. l1-1*n != (l1-1)*n
#define PC_NY_MIN (m1-1)
#define PC_NZ_MIN (n1-1)

#define PC_NX_MAX (l2-1)
#define PC_NY_MAX (m2-1)
#define PC_NZ_MAX (n2-1)

