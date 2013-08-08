#include <stdio.h>
#include <stdlib.h>
#include "cycle.h"
#include "sz.h"

int main()
{
  ticks t0, t1, s0, s1;
  double ds, ms = 1e100;

  FILE *file;

  int h;
  double *f, *df;

  f  = ((double *)malloc(sizeof(double) * S * S));
  df = ((double *)malloc(sizeof(double) * S * S));
  for(h = 0; h < S * S; ++h) f[h] = df[h] = h;
  f  += H * (S + 1);
  df += H * (S + 1);

  t0 = getticks();
  for(h = 0; h < 64; ++h) {
    s0 = getticks();
    kern(f, df);
    s1 = getticks();
    ds = elapsed(s1, s0);
    if(ms > ds) ms = ds;
  }
  t1 = getticks();
  ds = elapsed(t1, t0) / 64;

  df -= H * (S + 1);
  f  -= H * (S + 1);
  file = fopen("out", "w");
  fwrite(df, sizeof(double), S * S, file);
  fclose(file);
  free(df);
  free(f );

  printf("min(ticks) = %g; average(ticks) = %g\n", ms, ds);

  return 0;
}
