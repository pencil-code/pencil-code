#include "sz.h"

#define P(i,j) p[i+H][j+H]

void kern(const double *f, double *df)
{
  int i, j;

  double p[2 * H + 1][S];
  for(j = 0; j < S; ++j) p[1][j] = F(-3,j-H);
  for(j = 0; j < S; ++j) p[2][j] = F(-2,j-H);
  for(j = 0; j < S; ++j) p[3][j] = F(-1,j-H);
  for(j = 0; j < S; ++j) p[4][j] = F( 0,j-H);
  for(j = 0; j < S; ++j) p[5][j] = F( 1,j-H);
  for(j = 0; j < S; ++j) p[6][j] = F( 2,j-H);

  for(i = 0; i < N; ++i) {
    for(j = 0; j < S; ++j) p[0][j] = p[1][j];
    for(j = 0; j < S; ++j) p[1][j] = p[2][j];
    for(j = 0; j < S; ++j) p[2][j] = p[3][j];
    for(j = 0; j < S; ++j) p[3][j] = p[4][j];
    for(j = 0; j < S; ++j) p[4][j] = p[5][j];
    for(j = 0; j < S; ++j) p[5][j] = p[6][j];
    for(j = 0; j < S; ++j) p[6][j] = F(i+3,j-H);
    for(j = 0; j < N; ++j)
      df[i * S + j] += (45./60.) * (P(  1,j) - P( -1,j)) +
                       (-9./60.) * (P(  2,j) - P( -2,j)) +
                       ( 1./60.) * (P(  3,j) - P( -3,j)) +
                       (45./60.) * (P(0,j+1) - P(0,j-1)) +
                       (-9./60.) * (P(0,j+2) - P(0,j-2)) +
                       ( 1./60.) * (P(0,j+3) - P(0,j-3)) ;
  }
}
