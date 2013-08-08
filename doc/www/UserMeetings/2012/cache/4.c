#include "sz.h"

#define W 512
#define X (W + 2 * H)

#define P(i,j) p[i+H][j+H]

void kern(const double *f, double *df)
{
  int h, i, j;

  for(h = 0; h < N / W; ++h) {
    const int off = h * W;
    double p[2 * H + 1][X];
    for(j = 0; j < X; ++j) p[1][j] = F(-3,j+off-H);
    for(j = 0; j < X; ++j) p[2][j] = F(-2,j+off-H);
    for(j = 0; j < X; ++j) p[3][j] = F(-1,j+off-H);
    for(j = 0; j < X; ++j) p[4][j] = F( 0,j+off-H);
    for(j = 0; j < X; ++j) p[5][j] = F( 1,j+off-H);
    for(j = 0; j < X; ++j) p[6][j] = F( 2,j+off-H);

    for(i = 0; i < N; ++i) {
      for(j = 0; j < X; ++j) p[0][j] = p[1][j];
      for(j = 0; j < X; ++j) p[1][j] = p[2][j];
      for(j = 0; j < X; ++j) p[2][j] = p[3][j];
      for(j = 0; j < X; ++j) p[3][j] = p[4][j];
      for(j = 0; j < X; ++j) p[4][j] = p[5][j];
      for(j = 0; j < X; ++j) p[5][j] = p[6][j];
      for(j = 0; j < X; ++j) p[6][j] = F(i+3,j+off-H);
      for(j = 0; j < W; ++j)
        df[i * S + j + off] += (45./60.) * (P(  1,j) - P( -1,j)) +
                               (-9./60.) * (P(  2,j) - P( -2,j)) +
                               ( 1./60.) * (P(  3,j) - P( -3,j)) +
                               (45./60.) * (P(0,j+1) - P(0,j-1)) +
                               (-9./60.) * (P(0,j+2) - P(0,j-2)) +
                               ( 1./60.) * (P(0,j+3) - P(0,j-3)) ;
    }
  }
}
