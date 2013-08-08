#include "sz.h"

#define P(j) p[j + H]

void kern(const double *f, double *df)
{
  int i, j;
  for(i = 0; i < N; ++i) {
    double p[S];

    for(j = 0; j < S; ++j)
      p[j] = F(i,j-H);

    for(j = 0; j < N; ++j)
      df[i * S + j] += (45./60.) * (F(i+1,j) - F(i-1,j)) +
                       (-9./60.) * (F(i+2,j) - F(i-2,j)) +
                       ( 1./60.) * (F(i+3,j) - F(i-3,j)) +
                       (45./60.) * (P(  j+1) - P(  j-1)) +
                       (-9./60.) * (P(  j+2) - P(  j-2)) +
                       ( 1./60.) * (P(  j+3) - P(  j-3)) ;
  }
}
