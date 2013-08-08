#define N 1024
#define H 3
#define S (N + 2 * H)

#define F(i,j) f[(i) * S + (j)]

void kern(const double *f, double *df);
