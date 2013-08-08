#include "sz.h"

void kern(const double *f, double *df)
{
  int i, j;
  for(i = 0; i < N; ++i)
    for(j = 0; j < N; ++j)
      df[i * S + j] += (45./60.) * (F(i+1,j) - F(i-1,j)) +
                       (-9./60.) * (F(i+2,j) - F(i-2,j)) +
                       ( 1./60.) * (F(i+3,j) - F(i-3,j)) +
                       (45./60.) * (F(i,j+1) - F(i,j-1)) +
                       (-9./60.) * (F(i,j+2) - F(i,j-2)) +
                       ( 1./60.) * (F(i,j+3) - F(i,j-3)) ;
}
